from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd


def _env_truthy(value: Any) -> bool:
    return str(value).strip().lower() in ("1", "true", "yes", "on")


def _stable_json_dumps(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _sha1_text(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class ClusterDbConfig:
    enabled: bool
    db_path: str
    schema_version: int = 1

    ENV_ENABLE: str = "TACKLE_CLUSTER_DB"
    ENV_PATH: str = "TACKLE_CLUSTER_DB_PATH"

    @staticmethod
    def default_db_path(analysis_outpath: str) -> str:
        return os.path.join(
            os.path.abspath(analysis_outpath), "cluster2", "cluster2.sqlite"
        )

    @classmethod
    def from_env(
        cls,
        analysis_outpath: str,
        *,
        enabled_override: Optional[bool] = None,
        path_override: Optional[str] = None,
    ) -> "ClusterDbConfig":
        enabled_env = _env_truthy(os.getenv(cls.ENV_ENABLE, ""))
        enabled = enabled_env if enabled_override is None else bool(enabled_override)

        env_path = os.getenv(cls.ENV_PATH) or None
        db_path = path_override or env_path or cls.default_db_path(analysis_outpath)

        return cls(enabled=enabled, db_path=db_path, schema_version=cls.schema_version)


class ClusterDb:
    def __init__(self, config: ClusterDbConfig, *, logger=None):
        self.config = config
        self._logger = logger
        self._conn: Optional[sqlite3.Connection] = None

        if self.config.enabled:
            self._open()

    def close(self) -> None:
        if self._conn is None:
            return
        try:
            self._conn.close()
        finally:
            self._conn = None

    def _open(self) -> None:
        db_path = self.config.db_path
        Path(db_path).parent.mkdir(parents=True, exist_ok=True)

        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA journal_mode=WAL;")
        conn.execute("PRAGMA synchronous=NORMAL;")
        conn.execute("PRAGMA foreign_keys=ON;")

        version = int(conn.execute("PRAGMA user_version;").fetchone()[0])
        if version not in (0, self.config.schema_version):
            conn.close()
            raise RuntimeError(
                f"Unsupported cluster DB schema version {version} (expected {self.config.schema_version})"
            )
        if version == 0:
            self._init_schema(conn)
            conn.execute(f"PRAGMA user_version={self.config.schema_version};")

        self._conn = conn

    @staticmethod
    def _init_schema(conn: sqlite3.Connection) -> None:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS runs (
              run_id TEXT PRIMARY KEY,
              grid_id TEXT NOT NULL,
              tsv_path TEXT NOT NULL,
              created_at REAL NOT NULL,
              cluster_func TEXT,
              nclusters INTEGER,
              seed INTEGER,
              linkage TEXT,
              cluster_fillna TEXT,
              z_score TEXT,
              z_score_by TEXT,
              standard_scale TEXT,
              n_genes INTEGER,
              n_clusters INTEGER,
              sil_mean REAL,
              sil_q10 REAL,
              sil_neg_frac REAL,
              params_json TEXT
            );

            CREATE UNIQUE INDEX IF NOT EXISTS idx_runs_tsv_path ON runs(tsv_path);
            CREATE INDEX IF NOT EXISTS idx_runs_grid_id ON runs(grid_id);

            CREATE TABLE IF NOT EXISTS memberships (
              run_id TEXT NOT NULL,
              gene_id TEXT NOT NULL,
              gene_symbol TEXT,
              cluster INTEGER,
              neighbor INTEGER,
              sil_width REAL,
              row_order INTEGER,
              PRIMARY KEY (run_id, gene_id),
              FOREIGN KEY (run_id) REFERENCES runs(run_id) ON DELETE CASCADE
            );

            CREATE INDEX IF NOT EXISTS idx_memberships_run_cluster ON memberships(run_id, cluster);
            CREATE INDEX IF NOT EXISTS idx_memberships_gene_id ON memberships(gene_id);

            CREATE TABLE IF NOT EXISTS clusters (
              run_id TEXT NOT NULL,
              cluster INTEGER NOT NULL,
              n_genes INTEGER,
              sil_mean REAL,
              sil_q10 REAL,
              sil_neg_frac REAL,
              sil_median REAL,
              PRIMARY KEY (run_id, cluster),
              FOREIGN KEY (run_id) REFERENCES runs(run_id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS ari (
              grid_id TEXT NOT NULL,
              run_id_a TEXT NOT NULL,
              run_id_b TEXT NOT NULL,
              n_genes INTEGER,
              ari REAL,
              PRIMARY KEY (grid_id, run_id_a, run_id_b),
              FOREIGN KEY (run_id_a) REFERENCES runs(run_id) ON DELETE CASCADE,
              FOREIGN KEY (run_id_b) REFERENCES runs(run_id) ON DELETE CASCADE
            );
            """
        )

    def _log(self, level: str, msg: str, *args) -> None:
        if self._logger is None:
            return
        fn = getattr(self._logger, level, None)
        if fn is None:
            return
        try:
            fn(msg, *args)
        except Exception:
            pass

    @staticmethod
    def _normalize_optional_str(value: Any) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, float) and not np.isfinite(value):
            return None
        text = str(value)
        return text if text != "" else None

    @staticmethod
    def _compute_silhouette_metrics(series: pd.Series) -> dict[str, Any]:
        values = pd.to_numeric(series, errors="coerce")
        values = values[np.isfinite(values)]
        if len(values) == 0:
            return {
                "sil_mean": None,
                "sil_q10": None,
                "sil_neg_frac": None,
                "sil_median": None,
            }
        sil_mean = float(values.mean())
        sil_median = float(values.median())
        sil_q10 = float(values.quantile(0.10))
        sil_neg_frac = float((values < 0).mean())
        return {
            "sil_mean": sil_mean,
            "sil_q10": sil_q10,
            "sil_neg_frac": sil_neg_frac,
            "sil_median": sil_median,
        }

    @staticmethod
    def _compute_run_id(rel_tsv_path: str) -> str:
        return _sha1_text(rel_tsv_path)

    @staticmethod
    def _compute_grid_id(rel_tsv_path: str, *, cluster_func: str, nclusters: Any) -> str:
        p = Path(rel_tsv_path)
        stem = p.stem
        suffix = f"_{cluster_func}_{nclusters}"
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)] + f"_{cluster_func}_*"
        grid_key = str(p.with_name(stem + p.suffix))
        return _sha1_text(grid_key)

    def ingest_cluster2_metrics(
        self,
        metrics_df: pd.DataFrame,
        *,
        analysis_outpath: str,
        tsv_path: str,
        cluster_func: str,
        nclusters: Any,
        seed: Any = None,
        linkage: Any = None,
        cluster_fillna: Any = None,
        z_score: Any = None,
        z_score_by: Any = None,
        standard_scale: Any = None,
        extra_params: Optional[dict[str, Any]] = None,
    ) -> Optional[str]:
        if not self.config.enabled:
            return None
        if self._conn is None:
            return None

        abs_outpath = os.path.abspath(analysis_outpath)
        abs_tsv_path = os.path.abspath(tsv_path)
        rel_tsv_path = (
            os.path.relpath(abs_tsv_path, abs_outpath)
            if abs_tsv_path.startswith(abs_outpath + os.sep)
            else abs_tsv_path
        )

        run_id = self._compute_run_id(rel_tsv_path)
        grid_id = self._compute_grid_id(
            rel_tsv_path, cluster_func=cluster_func, nclusters=nclusters
        )

        df = metrics_df.copy()
        df.columns = [str(c).strip() for c in df.columns]
        required = {"GeneID", "GeneSymbol", "cluster", "sil_width"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"cluster2 metrics missing required columns: {sorted(missing)}"
            )

        df["GeneID"] = df["GeneID"].astype(str)
        df["GeneSymbol"] = df["GeneSymbol"].astype(str)
        df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce").astype("Int64")
        if "neighbor" in df.columns:
            df["neighbor"] = pd.to_numeric(df["neighbor"], errors="coerce").astype(
                "Int64"
            )
        else:
            df["neighbor"] = pd.Series([pd.NA] * len(df), dtype="Int64")
        df["sil_width"] = pd.to_numeric(df["sil_width"], errors="coerce")
        df = df.reset_index(drop=True)
        df["row_order"] = pd.RangeIndex(start=1, stop=len(df) + 1, step=1)

        n_genes = int(len(df))
        n_clusters = int(df["cluster"].dropna().nunique())
        sil = self._compute_silhouette_metrics(df["sil_width"])

        def _finite_sil(series: pd.Series) -> pd.Series:
            vals = pd.to_numeric(series, errors="coerce")
            return vals[np.isfinite(vals)]

        def _sil_q10(series: pd.Series) -> float:
            vals = _finite_sil(series)
            return float(vals.quantile(0.10)) if len(vals) else float("nan")

        def _sil_neg_frac(series: pd.Series) -> float:
            vals = _finite_sil(series)
            return float((vals < 0).mean()) if len(vals) else float("nan")

        cluster_stats = (
            df.groupby("cluster", dropna=False)
            .agg(
                n_genes=("GeneID", "size"),
                sil_width_mean=("sil_width", "mean"),
                sil_width_median=("sil_width", "median"),
                sil_width_q10=("sil_width", _sil_q10),
                sil_width_neg_frac=("sil_width", _sil_neg_frac),
            )
            .reset_index()
        )

        params = {
            "cluster_func": self._normalize_optional_str(cluster_func),
            "nclusters": nclusters,
            "seed": seed,
            "linkage": self._normalize_optional_str(linkage),
            "cluster_fillna": self._normalize_optional_str(cluster_fillna),
            "z_score": self._normalize_optional_str(z_score),
            "z_score_by": self._normalize_optional_str(z_score_by),
            "standard_scale": self._normalize_optional_str(standard_scale),
        }
        if extra_params:
            params.update(extra_params)
        params_json = _stable_json_dumps(params)

        now = time.time()

        conn = self._conn
        assert conn is not None

        with conn:
            conn.execute(
                """
                INSERT INTO runs (
                  run_id, grid_id, tsv_path, created_at,
                  cluster_func, nclusters, seed, linkage, cluster_fillna, z_score, z_score_by, standard_scale,
                  n_genes, n_clusters, sil_mean, sil_q10, sil_neg_frac, params_json
                ) VALUES (
                  ?, ?, ?, ?,
                  ?, ?, ?, ?, ?, ?, ?, ?,
                  ?, ?, ?, ?, ?, ?
                )
                ON CONFLICT(run_id) DO UPDATE SET
                  grid_id=excluded.grid_id,
                  tsv_path=excluded.tsv_path,
                  created_at=excluded.created_at,
                  cluster_func=excluded.cluster_func,
                  nclusters=excluded.nclusters,
                  seed=excluded.seed,
                  linkage=excluded.linkage,
                  cluster_fillna=excluded.cluster_fillna,
                  z_score=excluded.z_score,
                  z_score_by=excluded.z_score_by,
                  standard_scale=excluded.standard_scale,
                  n_genes=excluded.n_genes,
                  n_clusters=excluded.n_clusters,
                  sil_mean=excluded.sil_mean,
                  sil_q10=excluded.sil_q10,
                  sil_neg_frac=excluded.sil_neg_frac,
                  params_json=excluded.params_json
                """,
                (
                    run_id,
                    grid_id,
                    rel_tsv_path,
                    now,
                    self._normalize_optional_str(cluster_func),
                    int(nclusters) if nclusters is not None else None,
                    int(seed) if seed is not None else None,
                    self._normalize_optional_str(linkage),
                    self._normalize_optional_str(cluster_fillna),
                    self._normalize_optional_str(z_score),
                    self._normalize_optional_str(z_score_by),
                    self._normalize_optional_str(standard_scale),
                    n_genes,
                    n_clusters,
                    sil["sil_mean"],
                    sil["sil_q10"],
                    sil["sil_neg_frac"],
                    params_json,
                ),
            )

            conn.execute("DELETE FROM memberships WHERE run_id = ?", (run_id,))
            conn.execute("DELETE FROM clusters WHERE run_id = ?", (run_id,))

            membership_rows = [
                (
                    run_id,
                    row["GeneID"],
                    row["GeneSymbol"],
                    None if pd.isna(row["cluster"]) else int(row["cluster"]),
                    None if pd.isna(row["neighbor"]) else int(row["neighbor"]),
                    None if not np.isfinite(row["sil_width"]) else float(row["sil_width"]),
                    int(row["row_order"]),
                )
                for row in df.to_dict(orient="records")
            ]
            conn.executemany(
                """
                INSERT INTO memberships (
                  run_id, gene_id, gene_symbol, cluster, neighbor, sil_width, row_order
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                membership_rows,
            )

            cluster_rows = [
                (
                    run_id,
                    None if pd.isna(row["cluster"]) else int(row["cluster"]),
                    int(row["n_genes"]),
                    None
                    if not np.isfinite(row["sil_width_mean"])
                    else float(row["sil_width_mean"]),
                    None
                    if not np.isfinite(row["sil_width_q10"])
                    else float(row["sil_width_q10"]),
                    None
                    if not np.isfinite(row["sil_width_neg_frac"])
                    else float(row["sil_width_neg_frac"]),
                    None
                    if not np.isfinite(row["sil_width_median"])
                    else float(row["sil_width_median"]),
                )
                for row in cluster_stats.to_dict(orient="records")
            ]
            conn.executemany(
                """
                INSERT INTO clusters (
                  run_id, cluster, n_genes, sil_mean, sil_q10, sil_neg_frac, sil_median
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                cluster_rows,
            )

            self._update_ari_for_run(conn, grid_id, run_id, df[["GeneID", "cluster"]])

        self._log(
            "info",
            "Cluster DB ingest: %s (n_genes=%s, sil_mean=%s, sil_q10=%s, sil_neg_frac=%s)",
            self.config.db_path,
            n_genes,
            sil["sil_mean"],
            sil["sil_q10"],
            sil["sil_neg_frac"],
        )
        return run_id

    @staticmethod
    def _update_ari_for_run(
        conn: sqlite3.Connection,
        grid_id: str,
        run_id: str,
        membership: pd.DataFrame,
    ) -> None:
        from sklearn.metrics import adjusted_rand_score

        current = membership.copy()
        current["GeneID"] = current["GeneID"].astype(str)
        current["cluster"] = pd.to_numeric(current["cluster"], errors="coerce").astype(
            "Int64"
        )
        current = current.dropna(subset=["cluster"])
        if current.empty:
            return
        current_map = dict(zip(current["GeneID"], current["cluster"].astype(int)))

        other_runs = [
            row[0]
            for row in conn.execute(
                "SELECT run_id FROM runs WHERE grid_id = ? AND run_id != ?",
                (grid_id, run_id),
            ).fetchall()
        ]
        if not other_runs:
            return

        for other_id in other_runs:
            other_rows = conn.execute(
                "SELECT gene_id, cluster FROM memberships WHERE run_id = ? AND cluster IS NOT NULL",
                (other_id,),
            ).fetchall()
            if not other_rows:
                continue
            other_map = {str(g): int(c) for g, c in other_rows}
            common = sorted(set(current_map) & set(other_map))
            if len(common) < 2:
                continue
            labels_a = [current_map[g] for g in common]
            labels_b = [other_map[g] for g in common]
            ari = float(adjusted_rand_score(labels_a, labels_b))

            run_a, run_b = (run_id, other_id) if run_id < other_id else (other_id, run_id)
            conn.execute(
                """
                INSERT INTO ari (grid_id, run_id_a, run_id_b, n_genes, ari)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(grid_id, run_id_a, run_id_b) DO UPDATE SET
                  n_genes=excluded.n_genes,
                  ari=excluded.ari
                """,
                (grid_id, run_a, run_b, int(len(common)), ari),
            )

