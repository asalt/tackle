from __future__ import annotations

import json
import os
import shutil
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from .hashutil import sha256_file


SCHEMA_VERSION = 1


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def derive_analysis_dir(*, result_dir: str | Path, conf_path: str | Path, name: str) -> Path:
    base = Path(result_dir)
    conf_stem = Path(conf_path).stem
    return base / conf_stem / name


def _tackle_version() -> str:
    try:
        import importlib.metadata

        return importlib.metadata.version("tackle")
    except Exception:
        return "unknown"


@dataclass(frozen=True)
class BuildSpec:
    taxon: str = "all"
    non_zeros: int = 1
    unique_pepts: int = 0
    normalize_across_species: bool = False
    fill_na_zero: bool = True
    impute_missing_values: bool = False
    imputation_backend: str = "gaussian"
    lupine_mode: str = "local"
    norm_method: str = "none"  # one of: none, median, ifot, ifot_ki, ifot_tf, quantile75, quantile90, genefile
    genefile_norm: Optional[str] = None
    ref_group_cols: Optional[list[str]] = None
    ref_label_col: str = "label"
    ref_control_value: Optional[str] = None


@dataclass(frozen=True)
class AnalysisConfig:
    schema_version: int
    created_at: str
    tackle_version: str
    analysis_dir: str
    result_dir: str
    data_dir: str
    analysis_name: str
    run_name: str
    conf_original: str
    conf_copied: str
    conf_sha256: str
    build_spec: BuildSpec

    def to_json(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["build_spec"] = asdict(self.build_spec)
        return payload

    @classmethod
    def from_json(cls, payload: dict[str, Any]) -> "AnalysisConfig":
        build_spec = BuildSpec(**payload["build_spec"])
        return cls(
            schema_version=int(payload["schema_version"]),
            created_at=str(payload["created_at"]),
            tackle_version=str(payload.get("tackle_version", "unknown")),
            analysis_dir=str(payload["analysis_dir"]),
            result_dir=str(payload["result_dir"]),
            data_dir=str(payload["data_dir"]),
            analysis_name=str(payload["analysis_name"]),
            run_name=str(payload["run_name"]),
            conf_original=str(payload["conf_original"]),
            conf_copied=str(payload["conf_copied"]),
            conf_sha256=str(payload["conf_sha256"]),
            build_spec=build_spec,
        )


class AnalysisDir:
    def __init__(self, path: str | Path):
        self.path = Path(path)

    @property
    def analysis_json_path(self) -> Path:
        return self.path / "analysis.json"

    @property
    def manifest_sqlite_path(self) -> Path:
        return self.path / "manifest.sqlite"

    @property
    def manifest_json_path(self) -> Path:
        return self.path / "manifest.json"

    @property
    def inputs_dir(self) -> Path:
        return self.path / "inputs"

    @property
    def matrices_dir(self) -> Path:
        return self.path / "matrices"

    def ensure_layout(self) -> None:
        _safe_mkdir(self.path)
        _safe_mkdir(self.inputs_dir)
        _safe_mkdir(self.matrices_dir)

    def write_analysis_config(
        self,
        *,
        conf_path: str | Path,
        name: str,
        result_dir: str | Path,
        data_dir: str | Path,
        build_spec: BuildSpec,
        force: bool = False,
    ) -> AnalysisConfig:
        self.ensure_layout()
        if self.analysis_json_path.exists() and not force:
            raise FileExistsError(f"analysis.json already exists at {self.analysis_json_path}")

        src_conf = Path(conf_path)
        if not src_conf.exists():
            raise FileNotFoundError(str(src_conf))

        copied_name = src_conf.name
        copied_conf = self.inputs_dir / copied_name
        shutil.copy2(src_conf, copied_conf)

        cfg = AnalysisConfig(
            schema_version=SCHEMA_VERSION,
            created_at=_utc_now_iso(),
            tackle_version=_tackle_version(),
            analysis_dir=str(self.path.resolve()),
            result_dir=str(Path(result_dir).resolve()),
            data_dir=str(Path(data_dir).resolve()),
            analysis_name=src_conf.stem,
            run_name=str(name),
            conf_original=str(src_conf.resolve()),
            conf_copied=str(copied_conf.relative_to(self.path)),
            conf_sha256=sha256_file(copied_conf),
            build_spec=build_spec,
        )

        with self.analysis_json_path.open("w", encoding="utf-8") as handle:
            json.dump(cfg.to_json(), handle, sort_keys=True, indent=2)

        return cfg

    def load_analysis_config(self) -> AnalysisConfig:
        if not self.analysis_json_path.exists():
            raise FileNotFoundError(f"Missing analysis.json at {self.analysis_json_path}")
        with self.analysis_json_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        return AnalysisConfig.from_json(payload)

    def resolve_copied_conf_path(self, cfg: AnalysisConfig) -> Path:
        return self.path / cfg.conf_copied

    def canonical_matrix_paths(self) -> dict[str, Path]:
        return {
            "area_tsv": self.matrices_dir / "matrix.area.tsv",
            "gct": self.matrices_dir / "matrix.gct",
            "mspc_tsv": self.matrices_dir / "matrix.mspc.tsv",
        }

    def relpath(self, path: str | Path) -> str:
        try:
            return str(Path(path).resolve().relative_to(self.path.resolve()))
        except Exception:
            return str(path)
