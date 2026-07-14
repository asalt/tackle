from __future__ import annotations

import hashlib
import json
import os
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from .gct_io import (
    dataframe_content_hash,
    gct_io_contracts,
    read_gctx_content_hash,
    selected_content_hash_algorithm,
    write_gctx,
)
from .pca_replay_rmd import render_pca_replay_rmd, render_pca_replay_sh


PCA_REPLAY_CONTRACT_VERSION = 4
PCA_REPLAY_MATRIX_DTYPE = "float64"


@dataclass(frozen=True)
class PcaReplayFiles:
    replay_dir: Path
    run_id: str
    gct_path: Path
    context_path: Path
    rmd_path: Path
    stats_r_path: Path
    pointer_path: Path


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    if hasattr(value, "tolist"):
        try:
            return value.tolist()
        except Exception:
            pass
    return str(value)


def _hash_frame(frame: pd.DataFrame) -> str:
    return dataframe_content_hash(
        frame,
        serialization_contract="tackle-pca-replay-frame-v2",
        float_format="%.17g",
    )


def _clean_gct_metadata(frame: pd.DataFrame, ids: list[str]) -> pd.DataFrame:
    result = frame.copy()
    result.index = result.index.astype(str)
    result = result.reindex(ids)
    if "id" in result.columns:
        result = result.drop(columns=["id"])
    for column in result.columns:
        series = result[column]
        if isinstance(series.dtype, pd.CategoricalDtype):
            result[column] = series.astype(object).where(series.notna(), "").astype(str)
        elif pd.api.types.is_datetime64_any_dtype(series):
            result[column] = series.map(lambda x: "" if pd.isna(x) else x.isoformat())
        elif not pd.api.types.is_numeric_dtype(series):
            result[column] = series.map(lambda x: "" if pd.isna(x) else str(x))
    result.insert(0, "id", ids)
    result.index = ids
    return result


def _normalize_mapping(mapping: Mapping[str, Any] | None) -> dict[str, Any]:
    return json.loads(
        json.dumps(dict(mapping or {}), sort_keys=True, default=_json_default)
    )


def write_pca2_replay(
    *,
    pca_matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    feature_symbols: Mapping[Any, Any] | None,
    pca2_outname: str,
    analysis_outpath: str,
    preprocessing: Mapping[str, Any],
    plot_parameters: Mapping[str, Any],
    data_parameters: Mapping[str, Any],
    metadata_colors: Mapping[str, Mapping[str, str]] | None = None,
    separation_testing: Mapping[str, Any] | None = None,
) -> PcaReplayFiles:
    """Write an exact post-preprocessing, pre-prcomp pca2 replay bundle."""

    matrix = pca_matrix.copy().astype(float)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix.index.name = "sample"
    matrix.columns.name = "feature"
    if matrix.index.has_duplicates:
        raise ValueError("The pre-SVD PCA matrix contains duplicate sample identifiers.")
    if matrix.columns.has_duplicates:
        raise ValueError("The pre-SVD PCA matrix contains duplicate feature identifiers.")
    if matrix.isna().any().any():
        raise ValueError("The final pre-SVD PCA matrix unexpectedly contains NA values.")

    samples = list(matrix.index)
    features = list(matrix.columns)
    matrix_hash = _hash_frame(matrix)
    identity = {
        "contract_version": PCA_REPLAY_CONTRACT_VERSION,
        "matrix_hash_algorithm": selected_content_hash_algorithm(),
        "matrix_hash": matrix_hash,
        "preprocessing": _normalize_mapping(preprocessing),
        "plot_parameters": _normalize_mapping(plot_parameters),
        "data_parameters": _normalize_mapping(data_parameters),
        "separation_testing": _normalize_mapping(separation_testing),
    }
    run_id = hashlib.sha1(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:10]

    replay_root = Path(pca2_outname).expanduser().resolve().parent / "replay"
    replay_dir = replay_root / run_id
    replay_dir.mkdir(parents=True, exist_ok=True)

    sample_metadata = sample_metadata.copy()
    sample_metadata.index = sample_metadata.index.astype(str)
    sample_rdesc = _clean_gct_metadata(sample_metadata, samples)
    if "sample" not in sample_rdesc.columns:
        sample_rdesc.insert(1, "sample", samples)

    symbol_lookup = {
        str(key): "" if value is None or pd.isna(value) else str(value)
        for key, value in (feature_symbols or {}).items()
    }
    feature_cdesc = pd.DataFrame(
        {
            "GeneID": features,
            "GeneSymbol": [symbol_lookup.get(feature, "") for feature in features],
        },
        index=features,
    )
    feature_cdesc = _clean_gct_metadata(feature_cdesc, features)

    gct_path = write_gctx(
        matrix,
        replay_dir / "pca_input_pre_svd.gctx",
        row_metadata=sample_rdesc,
        col_metadata=feature_cdesc,
        matrix_dtype=PCA_REPLAY_MATRIX_DTYPE,
        content_addressed=True,
    )
    context_path = replay_dir / "pca_replay_context.json"
    rmd_path = replay_dir / "replot.Rmd"
    stats_r_path = replay_dir / "pca_stats.R"
    render_path = replay_dir / "render.sh"
    readme_path = replay_dir / "README.txt"

    context = {
        "replay_contract_version": PCA_REPLAY_CONTRACT_VERSION,
        "replay_scope": "pca2",
        "run_id": run_id,
        "created_at": datetime.now().astimezone().isoformat(),
        "authoritative_input": {
            "path": os.path.relpath(gct_path, replay_dir).replace("\\", "/"),
            "role": "final post-fill, post-center/scale matrix immediately before prcomp",
            "orientation": "rows are samples; columns are features",
            "storage_format": "gctx",
            "matrix_dtype": PCA_REPLAY_MATRIX_DTYPE,
            "gctx_writer_contract": gct_io_contracts()["gctx"]["writer_contract"],
            "gctx_content_hash_algorithm": selected_content_hash_algorithm(),
            "gctx_content_hash": read_gctx_content_hash(gct_path),
            "matrix_hash_algorithm": selected_content_hash_algorithm(),
            "matrix_hash": matrix_hash,
        },
        "matrix": {
            "n_samples": int(matrix.shape[0]),
            "n_features": int(matrix.shape[1]),
            "contains_na": False,
        },
        "sample_ids": samples,
        "feature_ids": features,
        "preprocessing": _normalize_mapping(preprocessing),
        "prcomp_arguments": {"center": False, "scale.": False},
        "plot_parameters": _normalize_mapping(plot_parameters),
        "data_parameters": _normalize_mapping(data_parameters),
        "metadata_colors": _normalize_mapping(metadata_colors),
        "separation_testing": _normalize_mapping(separation_testing),
        "replot_default_pc_pairs": [[1, 2], [1, 3], [2, 3]],
    }
    context_path.write_text(
        json.dumps(context, sort_keys=True, indent=2, default=_json_default),
        encoding="utf-8",
    )
    title = f"Tackle pca2 replay: {Path(pca2_outname).name}"
    rmd_path.write_text(
        render_pca_replay_rmd(
            title=title,
            gct_relpath=os.path.relpath(gct_path, replay_dir).replace("\\", "/"),
        ),
        encoding="utf-8",
    )
    shutil.copyfile(Path(__file__).resolve().parent / "R" / "pca_stats.R", stats_r_path)
    render_path.write_text(render_pca_replay_sh(), encoding="utf-8")
    render_path.chmod(0o755)
    readme_path.write_text(
        "This bundle stores the exact sample-by-feature matrix passed to prcomp by tackle pca2.\n"
        "Run ./render.sh from this directory to render replot.Rmd.\n"
        "The replay must not center, scale, fill, normalize, or transpose the GCTX matrix.\n"
        "pca_stats.R contains the package-independent Euclidean R2 and Welch-James replay code.\n",
        encoding="utf-8",
    )

    pointer_path = Path(analysis_outpath).expanduser().resolve() / "context" / "last_pca2_replay.json"
    pointer_path.parent.mkdir(parents=True, exist_ok=True)
    pointer_payload = {
        "replay_contract_version": PCA_REPLAY_CONTRACT_VERSION,
        "run_id": run_id,
        "replay_dir": str(replay_dir),
        "gct_path": str(gct_path),
        "context_path": str(context_path),
        "rmd_path": str(rmd_path),
        "stats_r_path": str(stats_r_path),
    }
    pointer_path.write_text(
        json.dumps(pointer_payload, sort_keys=True, indent=2), encoding="utf-8"
    )

    return PcaReplayFiles(
        replay_dir=replay_dir,
        run_id=run_id,
        gct_path=gct_path,
        context_path=context_path,
        rmd_path=rmd_path,
        stats_r_path=stats_r_path,
        pointer_path=pointer_path,
    )
