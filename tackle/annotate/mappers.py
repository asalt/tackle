from __future__ import annotations

import glob
import logging
import threading
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, Optional, Sequence

import pandas as pd


def _tackle_data_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "data"


def _truthy_membership_mask(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    norm = values.fillna("").astype(str).str.strip()
    norm_l = norm.str.lower()
    return (norm != "") & (~norm_l.isin(("0", "false", "na", "nan", "none")))


@dataclass(frozen=True, slots=True)
class LazyTableConfig:
    path: Path
    read_kwargs: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(
            self, "read_kwargs", MappingProxyType(dict(self.read_kwargs))
        )


@dataclass(slots=True)
class LazyTable:
    config: LazyTableConfig
    logger: Optional[logging.Logger] = None
    _df: Optional[pd.DataFrame] = field(default=None, init=False, repr=False)
    _lock: threading.Lock = field(default_factory=threading.Lock, init=False, repr=False)

    def read(self) -> pd.DataFrame:
        path = self.config.path
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        suffix = path.suffix.lower()
        read_kwargs = dict(self.config.read_kwargs)

        try:
            if suffix in (".tsv", ".tab", ".data"):
                return pd.read_table(path, **read_kwargs)
            if suffix in (".xlsx", ".xls"):
                return pd.read_excel(path, **read_kwargs)
        except Exception as e:
            raise RuntimeError(
                f"Failed to read {path}: {type(e).__name__}: {e}"
            ) from e

        raise ValueError(f"Unsupported file format for {path}")

    @property
    def df(self) -> pd.DataFrame:
        if self._df is None:
            with self._lock:
                if self._df is None:
                    if self.logger is not None:
                        try:
                            self.logger.info("Loading %s", self.config.path)
                        except Exception:
                            pass
                    self._df = self.read()
        return self._df


@dataclass(frozen=True, slots=True)
class GeneMapperConfig(LazyTableConfig):
    @classmethod
    def default(cls, *, data_dir: Optional[Path] = None) -> "GeneMapperConfig":
        base = data_dir or _tackle_data_dir()
        return cls(
            path=base / "genetable20201208_median_isoform_mass.tsv",
            read_kwargs={"dtype": {"GeneID": str, "TaxonID": str}, "index_col": "GeneID"},
        )


@dataclass(slots=True)
class GeneMapper(LazyTable):
    config: GeneMapperConfig

    _symbol: Optional[dict[str, str]] = field(default=None, init=False, repr=False)
    _funcat: Optional[dict[str, str]] = field(default=None, init=False, repr=False)
    _description: Optional[dict[str, str]] = field(default=None, init=False, repr=False)
    _taxon: Optional[dict[str, str]] = field(default=None, init=False, repr=False)

    @property
    def symbol(self) -> dict[str, str]:
        if self._symbol is None:
            series = self.df.get("GeneSymbol")
            if series is None or series.empty:
                self._symbol = {}
            else:
                self._symbol = series.astype(str).to_dict()
        return self._symbol

    @property
    def funcat(self) -> dict[str, str]:
        if self._funcat is None:
            series = self.df.get("FunCats")
            if series is None or series.empty:
                self._funcat = {}
            else:
                self._funcat = series.fillna("").astype(str).to_dict()
        return self._funcat

    @property
    def description(self) -> dict[str, str]:
        if self._description is None:
            series = self.df.get("GeneDescription")
            if series is None or series.empty:
                self._description = {}
            else:
                self._description = series.fillna("").astype(str).to_dict()
        return self._description

    @property
    def taxon(self) -> dict[str, str]:
        if self._taxon is None:
            series = self.df.get("TaxonID")
            if series is None or series.empty:
                self._taxon = {}
            else:
                self._taxon = series.astype(str).to_dict()
        return self._taxon


@dataclass(frozen=True, slots=True)
class HomologeneMapperConfig(LazyTableConfig):
    @classmethod
    def default(cls, *, data_dir: Optional[Path] = None) -> "HomologeneMapperConfig":
        base = data_dir or _tackle_data_dir()
        candidates = sorted(
            glob.glob(str(base / "homologene*.data")), reverse=True
        ) or sorted(glob.glob(str(base / "homologene*data")), reverse=True)
        if not candidates:
            raise FileNotFoundError(f"No homologene data found under {base}")
        return cls(
            path=Path(candidates[0]),
            read_kwargs={
                "header": None,
                "dtype": str,
                "names": (
                    "Homologene",
                    "TaxonID",
                    "GeneID",
                    "Symbol",
                    "ProteinGI",
                    "ProteinAccession",
                ),
            },
        )


@dataclass(slots=True)
class HomologeneMapper(LazyTable):
    config: HomologeneMapperConfig

    def map_to_human(
        self,
        gene_ids: Sequence[str],
        *,
        human_taxon: str = "9606",
    ) -> dict[str, Optional[str]]:
        gids = [str(x) for x in gene_ids]
        df = self.df

        query = df[df["GeneID"].isin(gids)]
        gid_to_homologene = (
            query[["GeneID", "Homologene"]].set_index("GeneID")["Homologene"].to_dict()
        )
        homologene_to_human = (
            df[df["TaxonID"] == str(human_taxon)][["GeneID", "Homologene"]]
            .set_index("Homologene")["GeneID"]
            .to_dict()
        )
        return {gid: homologene_to_human.get(gid_to_homologene.get(gid)) for gid in gids}


@dataclass(frozen=True, slots=True)
class AnnotationTableConfig(LazyTableConfig):
    @classmethod
    def default(cls, *, data_dir: Optional[Path] = None) -> "AnnotationTableConfig":
        base = data_dir or _tackle_data_dir()
        return cls(path=base / "combined_annotations_new.tsv", read_kwargs={"dtype": str})


@dataclass(slots=True)
class AnnotationTable(LazyTable):
    config: AnnotationTableConfig

    _categories: Optional[list[str]] = field(default=None, init=False, repr=False)

    @property
    def categories(self) -> list[str]:
        if self._categories is None:
            if self._df is not None and isinstance(self._df, pd.DataFrame):
                cols = list(self._df.columns)
            else:
                try:
                    with open(self.config.path, "r") as handle:
                        cols = handle.readline().rstrip("\n").split("\t")
                except OSError:
                    cols = []
            self._categories = [
                c for c in cols if c and c not in ("GeneID", "GeneSymbol")
            ]
        return self._categories

    def get_annotation(self, annotation: str) -> pd.DataFrame:
        df = self.df
        if annotation in (None, "", "_all"):
            return df
        if annotation not in df.columns:
            raise KeyError(f"Unknown annotation column: {annotation}")
        mask = _truthy_membership_mask(df[annotation])
        return df.loc[mask]

    def map_gene_ids(
        self,
        gene_ids: Sequence[str],
        *,
        field: str = "GeneID",
        taxon: str = "9606",
        fallback_to_human: bool = True,
        homologene: Optional[HomologeneMapper] = None,
        human_taxon: str = "9606",
    ) -> pd.DataFrame:
        gids = [str(x) for x in gene_ids]
        unique_gids = list(dict.fromkeys(gids))

        df = self.df
        if field not in df.columns:
            raise ValueError(f"{field} not found in annotation table")

        found = df[df[field].astype(str).isin(unique_gids)]
        drop_cols = [c for c in ("MitoCarta_Pathways",) if c in found.columns]
        found = found.drop(columns=drop_cols, errors="ignore")

        result = found.copy()
        missing = set(unique_gids) - set(found[field].astype(str).tolist())

        if (
            fallback_to_human
            and str(taxon) != str(human_taxon)
            and missing
        ):
            homologene = homologene or get_homologene_mapper()
            mapped = homologene.map_to_human(sorted(missing), human_taxon=human_taxon)
            mapped_ids = {k: v for k, v in mapped.items() if v is not None}
            if mapped_ids:
                map_df = pd.DataFrame(
                    {"original_gid": list(mapped_ids.keys()), "human_gid": list(mapped_ids.values())}
                )
                df_human = df[df["GeneID"].astype(str).isin(map_df["human_gid"].astype(str))]
                merged = map_df.merge(df_human, left_on="human_gid", right_on="GeneID", how="left")
                merged = merged.drop(columns=["GeneID", "human_gid"], errors="ignore").rename(
                    columns={"original_gid": "GeneID"}
                )
                merged = merged.drop(columns=drop_cols, errors="ignore")
                result = pd.concat([result, merged], ignore_index=True, sort=False)

        if "GeneID" not in result.columns and field == "GeneID":
            result = result.reset_index().rename(columns={"index": "GeneID"})

        if "GeneID" in result.columns:
            result["GeneID"] = result["GeneID"].astype(str)
            result = result.drop_duplicates(subset=["GeneID"], keep="first")

        front = [c for c in ("GeneID", "GeneSymbol") if c in result.columns]
        remaining = [c for c in result.columns if c not in front]
        return result[front + remaining]


@lru_cache(maxsize=1)
def get_gene_mapper() -> GeneMapper:
    return GeneMapper(GeneMapperConfig.default())


@lru_cache(maxsize=1)
def get_annotation_table() -> AnnotationTable:
    return AnnotationTable(AnnotationTableConfig.default())


@lru_cache(maxsize=1)
def get_homologene_mapper() -> HomologeneMapper:
    return HomologeneMapper(HomologeneMapperConfig.default())
