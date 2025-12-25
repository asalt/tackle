from __future__ import annotations

import glob
import gzip
import hashlib
import json
import logging
import os
import pickle
import time
from dataclasses import dataclass
from typing import Any, ClassVar, Iterable, Optional, Tuple


def _env_truthy(value: Any) -> bool:
    return str(value).strip().lower() in ("1", "true", "yes", "on")


def _stable_json_dumps(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


@dataclass(frozen=True)
class CacheConfig:
    enabled: bool
    cache_dir: str
    max_entries: int
    max_bytes: int
    bust: Optional[str]
    schema_version: int = 1

    ENV_ENABLE: ClassVar[str] = "TACKLE_CACHE"
    ENV_DIR: ClassVar[str] = "TACKLE_CACHE_DIR"
    ENV_MAX_ENTRIES: ClassVar[str] = "TACKLE_CACHE_MAX_ENTRIES"
    ENV_MAX_BYTES: ClassVar[str] = "TACKLE_CACHE_MAX_BYTES"
    ENV_BUST: ClassVar[str] = "TACKLE_CACHE_BUST"
    DEFAULT_MAX_ENTRIES: ClassVar[int] = 5
    DEFAULT_MAX_BYTES: ClassVar[int] = 2 * 1024 * 1024 * 1024
    DEFAULT_DIRNAME: ClassVar[str] = ".tackle_cache"

    @staticmethod
    def default_cache_dir(base_dir: Optional[str]) -> str:
        base_dir = base_dir or "."
        abs_base = os.path.abspath(base_dir)
        parent = os.path.dirname(abs_base)
        return os.path.join(parent, CacheConfig.DEFAULT_DIRNAME)

    @classmethod
    def from_env(cls, base_dir: Optional[str]) -> CacheConfig:
        enabled = _env_truthy(os.getenv(cls.ENV_ENABLE, ""))

        override_dir = os.getenv(cls.ENV_DIR)
        cache_dir = override_dir or cls.default_cache_dir(base_dir)

        try:
            max_entries = int(os.getenv(cls.ENV_MAX_ENTRIES, cls.DEFAULT_MAX_ENTRIES))
        except ValueError:
            max_entries = cls.DEFAULT_MAX_ENTRIES

        try:
            max_bytes = int(os.getenv(cls.ENV_MAX_BYTES, cls.DEFAULT_MAX_BYTES))
        except ValueError:
            max_bytes = cls.DEFAULT_MAX_BYTES

        bust = os.getenv(cls.ENV_BUST) or None

        return cls(
            enabled=enabled,
            cache_dir=cache_dir,
            max_entries=max_entries,
            max_bytes=max_bytes,
            bust=bust,
            schema_version=cls.schema_version,
        )


class SnapshotCache:
    def __init__(self, config: CacheConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self._logger = logger

    @staticmethod
    def hash_file(path: Optional[str]) -> Optional[str]:
        if not path:
            return None
        if not os.path.exists(path):
            return None
        hasher = hashlib.sha1()
        try:
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    hasher.update(chunk)
        except OSError:
            return None
        return hasher.hexdigest()

    @staticmethod
    def hash_iterable(values: Optional[Iterable[Any]]) -> Tuple[Optional[str], int]:
        if not values:
            return None, 0
        hasher = hashlib.sha1()
        items = sorted(str(v) for v in values)
        for item in items:
            hasher.update(item.encode("utf-8"))
            hasher.update(b"\0")
        return hasher.hexdigest(), len(items)

    @staticmethod
    def normalize_value(value: Any) -> Any:
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        if isinstance(value, os.PathLike):
            return os.fspath(value)
        if isinstance(value, dict):
            return {
                str(k): SnapshotCache.normalize_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        if isinstance(value, (list, tuple, set)):
            items = [SnapshotCache.normalize_value(v) for v in value]
            if isinstance(value, set):
                items = sorted(items, key=lambda x: str(x))
            return items
        return str(value)

    @classmethod
    def hash_key(cls, key: Any) -> Tuple[str, Any]:
        normalized = cls.normalize_value(key)
        payload = _stable_json_dumps(normalized)
        return hashlib.sha1(payload.encode("utf-8")).hexdigest(), normalized

    @staticmethod
    def read_meta(path: str) -> Optional[dict[str, Any]]:
        try:
            with open(path, "r") as handle:
                return json.load(handle)
        except (OSError, json.JSONDecodeError):
            return None

    @staticmethod
    def write_meta(path: str, meta: dict[str, Any]) -> bool:
        try:
            with open(path, "w") as handle:
                json.dump(meta, handle, sort_keys=True, indent=2)
        except OSError:
            return False
        return True

    @staticmethod
    def evict_entries(cache_dir: str, max_entries: int, max_bytes: int) -> None:
        entries = []
        meta_files = glob.glob(os.path.join(cache_dir, "cache_*.json"))
        for meta_path in meta_files:
            meta = SnapshotCache.read_meta(meta_path)
            if not meta:
                continue
            key_hash = meta.get("key_hash")
            if not key_hash:
                continue
            data_path = os.path.join(cache_dir, f"cache_{key_hash}.pkl.gz")
            if not os.path.exists(data_path):
                try:
                    os.remove(meta_path)
                except OSError:
                    pass
                continue
            size_bytes = os.path.getsize(data_path) + os.path.getsize(meta_path)
            last_access = meta.get("last_access", meta.get("created_at", 0))
            entries.append(
                {
                    "key_hash": key_hash,
                    "meta_path": meta_path,
                    "data_path": data_path,
                    "size_bytes": size_bytes,
                    "last_access": last_access,
                }
            )

        entries.sort(key=lambda entry: entry["last_access"])
        total_bytes = sum(entry["size_bytes"] for entry in entries)
        while entries and (len(entries) > max_entries or total_bytes > max_bytes):
            entry = entries.pop(0)
            for path in (entry["data_path"], entry["meta_path"]):
                try:
                    os.remove(path)
                except OSError:
                    pass
            total_bytes -= entry["size_bytes"]

    def _paths(self, key_hash: str) -> Tuple[str, str]:
        base = os.path.join(self.config.cache_dir, f"cache_{key_hash}")
        return base + ".pkl.gz", base + ".json"

    def load(self, key: Any) -> Optional[Any]:
        if not self.config.enabled:
            return None
        key_hash, _ = self.hash_key(key)
        data_path, meta_path = self._paths(key_hash)
        if not os.path.exists(data_path) or not os.path.exists(meta_path):
            return None
        meta = self.read_meta(meta_path)
        if not meta:
            return None
        if meta.get("schema_version") != self.config.schema_version:
            return None
        if meta.get("key_hash") != key_hash:
            return None
        try:
            with gzip.open(data_path, "rb") as handle:
                payload = pickle.load(handle)
        except (OSError, EOFError, pickle.UnpicklingError):
            return None

        meta["last_access"] = time.time()
        meta["size_bytes"] = os.path.getsize(data_path)
        self.write_meta(meta_path, meta)
        if self._logger is not None:
            self._logger.info("Cache hit: %s", data_path)
        return payload

    def save(self, key: Any, payload: Any) -> bool:
        if not self.config.enabled:
            return False
        os.makedirs(self.config.cache_dir, exist_ok=True)
        key_hash, normalized_key = self.hash_key(key)
        data_path, meta_path = self._paths(key_hash)
        try:
            with gzip.open(data_path, "wb") as handle:
                pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
        except OSError:
            return False

        now = time.time()
        meta = {
            "schema_version": self.config.schema_version,
            "key_hash": key_hash,
            "created_at": now,
            "last_access": now,
            "size_bytes": os.path.getsize(data_path),
            "key": normalized_key,
        }
        self.write_meta(meta_path, meta)
        self.evict_entries(
            self.config.cache_dir, self.config.max_entries, self.config.max_bytes
        )
        if self._logger is not None:
            self._logger.info("Cache write: %s", data_path)
        return True


@dataclass(frozen=True)
class DataSnapshot:
    data: Any
    df_filtered: Any
    col_metadata: Any
    taxon_ratios: Any
    gid_funcat_mapping: Any
    metric_values: Any
    areas: Any
    areas_log: Any
    mask: Any
    zeros: Any
    minval: Any
    batch_applied: Any
    normed: Any

    @classmethod
    def from_data(cls, data_obj: Any) -> "DataSnapshot":
        metric_values = getattr(data_obj, "_metric_values", None)
        normalized_metric_values = None
        if metric_values is not None:
            normalized_metric_values = {k: dict(v) for k, v in metric_values.items()}

        return cls(
            data=getattr(data_obj, "data", None),
            df_filtered=getattr(data_obj, "df_filtered", None),
            col_metadata=getattr(data_obj, "col_metadata", None),
            taxon_ratios=getattr(data_obj, "taxon_ratios", None),
            gid_funcat_mapping=getattr(data_obj, "gid_funcat_mapping", None),
            metric_values=normalized_metric_values,
            areas=getattr(data_obj, "_areas", None),
            areas_log=getattr(data_obj, "_areas_log", None),
            mask=getattr(data_obj, "_mask", None),
            zeros=getattr(data_obj, "_zeros", None),
            minval=getattr(data_obj, "minval", None),
            batch_applied=getattr(data_obj, "batch_applied", None),
            normed=getattr(data_obj, "normed", False),
        )

    def to_payload(self) -> dict[str, Any]:
        return {
            "data": self.data,
            "df_filtered": self.df_filtered,
            "col_metadata": self.col_metadata,
            "taxon_ratios": self.taxon_ratios,
            "gid_funcat_mapping": self.gid_funcat_mapping,
            "metric_values": self.metric_values,
            "areas": self.areas,
            "areas_log": self.areas_log,
            "mask": self.mask,
            "zeros": self.zeros,
            "minval": self.minval,
            "batch_applied": self.batch_applied,
            "normed": self.normed,
        }

    @classmethod
    def from_payload(cls, payload: Any) -> "DataSnapshot":
        if isinstance(payload, DataSnapshot):
            return payload

        if not isinstance(payload, dict):
            raise TypeError("Cache payload is not a dict")

        return cls(
            data=payload["data"],
            df_filtered=payload["df_filtered"],
            col_metadata=payload["col_metadata"],
            taxon_ratios=payload["taxon_ratios"],
            gid_funcat_mapping=payload["gid_funcat_mapping"],
            metric_values=payload.get("metric_values"),
            areas=payload["areas"],
            areas_log=payload["areas_log"],
            mask=payload["mask"],
            zeros=payload["zeros"],
            minval=payload.get("minval"),
            batch_applied=payload.get("batch_applied"),
            normed=payload.get("normed", False),
        )

    def apply_to(self, data_obj: Any) -> None:
        data_obj.data = self.data
        data_obj.df_filtered = self.df_filtered
        data_obj.col_metadata = self.col_metadata
        data_obj.taxon_ratios = self.taxon_ratios
        data_obj.gid_funcat_mapping = self.gid_funcat_mapping
        data_obj._metric_values = self.metric_values
        data_obj._areas = self.areas
        data_obj._areas_log = self.areas_log
        data_obj._mask = self.mask
        data_obj._zeros = self.zeros
        data_obj.minval = self.minval
        data_obj.batch_applied = self.batch_applied
        data_obj.normed = self.normed

        if hasattr(data_obj, "_gid_symbol"):
            data_obj._gid_symbol = None
        if hasattr(data_obj, "panel"):
            data_obj.panel = None
        if hasattr(data_obj, "exps"):
            data_obj.exps = None


class DataSnapshotCache:
    def __init__(
        self,
        config: CacheConfig,
        *,
        only_local: bool,
        logger: Optional[logging.Logger] = None,
    ):
        self.config = config
        self.only_local = bool(only_local)
        self._logger = logger
        self._cache = SnapshotCache(config, logger=logger)
        self._key: Optional[dict[str, Any]] = None
        self.hit = False
        self.saved = False

    @classmethod
    def from_env(
        cls,
        base_dir: Optional[str],
        *,
        only_local: bool,
        logger: Optional[logging.Logger] = None,
    ) -> "DataSnapshotCache":
        return cls(CacheConfig.from_env(base_dir), only_local=only_local, logger=logger)

    def try_restore(self, data_obj: Any) -> bool:
        if not self.config.enabled:
            return False
        self._key = self.build_key(data_obj)
        payload = self._cache.load(self._key)
        if payload is None:
            return False
        try:
            DataSnapshot.from_payload(payload).apply_to(data_obj)
        except Exception:
            return False
        self.hit = True
        return True

    def maybe_save(self, data_obj: Any) -> bool:
        if not self.config.enabled or self.hit or self.saved:
            return False
        if self._key is None:
            self._key = self.build_key(data_obj)
        try:
            payload = DataSnapshot.from_data(data_obj).to_payload()
        except Exception:
            return False
        self.saved = self._cache.save(self._key, payload)
        return self.saved

    def build_key(self, data_obj: Any) -> dict[str, Any]:
        geneid_hash, geneid_count = SnapshotCache.hash_iterable(
            getattr(data_obj, "geneid_subset", None)
        )
        ignore_hash, ignore_count = SnapshotCache.hash_iterable(
            getattr(data_obj, "ignore_geneid_subset", None)
        )

        genefile_norm = getattr(data_obj, "genefile_norm", None)
        genefile_norm_hash = SnapshotCache.hash_file(genefile_norm) if genefile_norm else None

        additional_info = getattr(data_obj, "additional_info", None)
        additional_info_hash = (
            SnapshotCache.hash_file(additional_info) if additional_info else None
        )

        funcats = getattr(data_obj, "funcats", None)
        if funcats and not isinstance(funcats, str):
            funcats = sorted(funcats)
        funcats_inverse = getattr(data_obj, "funcats_inverse", None)
        if funcats_inverse and not isinstance(funcats_inverse, str):
            funcats_inverse = sorted(funcats_inverse)
        cluster_annotate_cols = getattr(data_obj, "cluster_annotate_cols", None)
        if cluster_annotate_cols and not isinstance(cluster_annotate_cols, str):
            cluster_annotate_cols = sorted(cluster_annotate_cols)

        experiment_file = getattr(data_obj, "experiment_file")
        return {
            "schema_version": self.config.schema_version,
            "config_path": os.path.abspath(experiment_file),
            "config_hash": SnapshotCache.hash_file(experiment_file),
            "additional_info_hash": additional_info_hash,
            "data_dir": os.path.abspath(getattr(data_obj, "data_dir")),
            "only_local": self.only_local,
            "taxon": getattr(data_obj, "taxon"),
            "funcats": funcats,
            "funcats_inverse": funcats_inverse,
            "geneid_subset_hash": geneid_hash,
            "geneid_subset_count": geneid_count,
            "ignore_geneid_subset_hash": ignore_hash,
            "ignore_geneid_subset_count": ignore_count,
            "unique_pepts": getattr(data_obj, "unique_pepts"),
            "non_zeros": getattr(data_obj, "non_zeros"),
            "nonzero_subgroup": SnapshotCache.normalize_value(
                getattr(data_obj, "nonzero_subgroup")
            ),
            "SRA": getattr(data_obj, "SRA"),
            "number_sra": getattr(data_obj, "number_sra"),
            "normalize_across_species": getattr(data_obj, "normalize_across_species"),
            "export_all": getattr(data_obj, "export_all"),
            "cluster_annotate_cols": cluster_annotate_cols,
            "metrics": getattr(data_obj, "metrics"),
            "metrics_after_filter": getattr(data_obj, "metrics_after_filter"),
            "metrics_unnormed_area": getattr(data_obj, "metrics_unnormed_area"),
            "ifot": getattr(data_obj, "ifot"),
            "ifot_ki": getattr(data_obj, "ifot_ki"),
            "ifot_tf": getattr(data_obj, "ifot_tf"),
            "median": getattr(data_obj, "median"),
            "quantile75": getattr(data_obj, "quantile75"),
            "quantile90": getattr(data_obj, "quantile90"),
            "genefile_norm": os.path.abspath(genefile_norm) if genefile_norm else None,
            "genefile_norm_hash": genefile_norm_hash,
            "impute_missing_values": getattr(data_obj, "impute_missing_values"),
            "fill_na_zero": getattr(data_obj, "fill_na_zero"),
            "imputation_backend": getattr(data_obj, "imputation_backend")
            if getattr(data_obj, "impute_missing_values")
            else None,
            "lupine_mode": getattr(data_obj, "lupine_mode")
            if getattr(data_obj, "impute_missing_values")
            else None,
            "batch": getattr(data_obj, "batch"),
            "batch_nonparametric": getattr(data_obj, "batch_nonparametric"),
            "batch_noimputation": getattr(data_obj, "batch_noimputation"),
            "covariate": getattr(data_obj, "covariate"),
            "norm_info": SnapshotCache.normalize_value(getattr(data_obj, "norm_info")),
            "cache_bust": self.config.bust,
        }
