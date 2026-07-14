from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from .gct_io import (
    dataframe_content_hash,
    gct_io_contracts,
    read_gctx_content_hash,
    selected_content_hash_algorithm,
    write_gctx,
)
from .limma_replay_rmd import render_limma_replay_rmd, render_limma_replay_sh
from .utils import _get_logger

logger = _get_logger(__name__)

REPLAY_CONTRACT_VERSION = 4
REPLAY_GCT_PRECISION = 17  # Retained for legacy text-GCT discovery/context readers.
REPLAY_MATRIX_DTYPE = "float64"


@dataclass(frozen=True)
class LimmaReplayFiles:
    volcano_dir: Path
    replay_dir: Path
    run_id: str
    gct_path: Path
    context_path: Path
    pointer_path: Path
    pre_impute_gct_path: Optional[Path] = None


def _sanitize_name(value: str) -> str:
    # Mirror tackle.statmodels.limma_runner._sanitize_name
    s = str(value)
    s = s.replace(":", "_")
    s = s.replace(" ", "_")
    s = s.replace("-", "_")
    s = s.replace("+", "_")
    s = s.replace("/", "_")
    s = s.replace("?", "qmk")
    return s


def _json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "item"):
        try:
            return obj.item()
        except Exception:
            pass
    if hasattr(obj, "tolist"):
        try:
            return obj.tolist()
        except Exception:
            pass
    return str(obj)


def _relative_path(path: Path, start: Path) -> str:
    rel = os.path.relpath(Path(path).resolve(), Path(start).resolve()).replace("\\", "/")
    return "." if rel == "." else rel


def _compute_run_id(params: Mapping[str, Any]) -> str:
    payload = json.dumps(params, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:10]


def _hash_frame(frame: pd.DataFrame) -> str:
    return dataframe_content_hash(
        frame,
        serialization_contract="tackle-limma-replay-frame-v2",
        float_format="%.17g",
    )


def _gct_output_base(out_path: Path) -> Path:
    base = str(out_path)
    base_lower = base.lower()
    for suffix in (".gct.gz", ".gctx", ".gct"):
        if base_lower.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return Path(base)


def _find_existing_gct(out_path: Path) -> Optional[Path]:
    out_base = _gct_output_base(Path(out_path))
    candidates = []
    for suffix in (".gctx", ".gct"):
        exact = out_base.with_suffix(suffix)
        if exact.exists() and exact.stat().st_size > 0:
            candidates.append(exact)
        for path in out_base.parent.glob(out_base.name + "*" + suffix):
            if not path.exists() or path.stat().st_size <= 0:
                continue
            # ``limma_input*`` also matches the separately written
            # ``limma_input_pre_impute*`` audit artifact. Never select that as
            # the authoritative modeled matrix.
            if (
                out_base.name == "limma_input"
                and path.name.startswith("limma_input_pre_impute")
            ):
                continue
            candidates.append(path)
    if not candidates:
        return None
    return max(set(candidates), key=lambda p: p.stat().st_mtime)


def _split_result_label(label: str) -> Tuple[str, str]:
    text = str(label)
    if "=" in text:
        left, right = text.split("=", 1)
        return left.strip(), right.strip()
    return text.strip(), text.strip()


def parse_limma_result_labels(labels: Iterable[str]) -> Tuple[List[str], List[str]]:
    """Return (direct_coef_list, contrast_expr_list) from run_limma_pipeline labels."""
    direct_coefs: List[str] = []
    contrast_exprs: List[str] = []
    for label in labels:
        left, right = _split_result_label(label)
        if left.startswith("coef_"):
            direct_coefs.append(right)
        else:
            contrast_exprs.append(right)

    # Stable dedupe preserving order
    def _uniq(seq: Sequence[str]) -> List[str]:
        out: List[str] = []
        seen = set()
        for x in seq:
            if x in seen:
                continue
            seen.add(x)
            out.append(x)
        return out

    return _uniq(direct_coefs), _uniq(contrast_exprs)


def _extract_limma_input_from_results(
    results: Mapping[str, pd.DataFrame],
    *,
    sample_ids: Sequence[str],
) -> pd.DataFrame:
    if not results:
        raise ValueError("No limma results provided (expected dict of contrast->DataFrame).")
    first_key = next(iter(results.keys()))
    df = results[first_key]
    sample_cols = [c for c in sample_ids if c in df.columns]
    if not sample_cols:
        raise ValueError(
            "Unable to locate expression columns in limma result table. "
            "Expected sample columns to be present."
        )
    return df[sample_cols].copy()


def _make_rdesc(
    gene_ids: Sequence[Any],
    *,
    gene_symbols: Optional[Mapping[Any, Any]] = None,
    gene_descriptions: Optional[Mapping[str, str]] = None,
    funcats: Optional[Mapping[Any, Any]] = None,
) -> pd.DataFrame:
    rows: Dict[str, List[Any]] = {"GeneID": [str(g) for g in gene_ids]}
    if gene_symbols is not None:
        rows["GeneSymbol"] = [str(gene_symbols.get(g, "")) if gene_symbols.get(g, "") is not None else "" for g in gene_ids]
    if funcats is not None:
        rows["FunCats"] = [str(funcats.get(g, "")) if funcats.get(g, "") is not None else "" for g in gene_ids]
    if gene_descriptions is not None:
        rows["GeneDescription"] = [str(gene_descriptions.get(str(g), "")) for g in gene_ids]
    rdesc = pd.DataFrame(rows)
    rdesc["id"] = rdesc["GeneID"]
    rdesc.index = rdesc["GeneID"]
    return rdesc


def _write_gct(
    *,
    out_path: Path,
    mat: pd.DataFrame,
    cdesc: pd.DataFrame,
    rdesc: pd.DataFrame,
    precision: int = REPLAY_GCT_PRECISION,
) -> Path:
    """Write a content-addressed float64 replay GCTX and return its exact path."""

    _ = precision  # Compatibility with older internal callers/tests.
    out_path = _gct_output_base(Path(out_path)).with_suffix(".gctx")
    return write_gctx(
        mat,
        out_path,
        row_metadata=rdesc,
        col_metadata=cdesc,
        matrix_dtype=REPLAY_MATRIX_DTYPE,
        content_addressed=True,
    )


def _write_gct_reuse_existing(
    *,
    out_path: Path,
    mat: pd.DataFrame,
    cdesc: pd.DataFrame,
    rdesc: pd.DataFrame,
    precision: int = REPLAY_GCT_PRECISION,
) -> Path:
    # The first-class writer checks the full root-attribute digest and its
    # recorded algorithm, then skips I/O only for the exact requested payload.
    # Broad filename/mtime reuse is reserved for legacy replay discovery.
    return _write_gct(
        out_path=out_path,
        mat=mat,
        cdesc=cdesc,
        rdesc=rdesc,
        precision=precision,
    )


def _write_json(path: Path, payload: Mapping[str, Any], *, force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        raise FileExistsError(f"{path} already exists (use --force to overwrite)")
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n", encoding="utf-8")


def _write_text(path: Path, content: str, *, force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        raise FileExistsError(f"{path} already exists (use --force to overwrite)")
    path.write_text(content, encoding="utf-8")




def write_limma_replay_files(
    *,
    analysis_dir: str,
    volcano_dir: str,
    results: Mapping[str, pd.DataFrame],
    sample_metadata: pd.DataFrame,
    expression_matrix: Optional[pd.DataFrame] = None,
    expression_matrix_source: Optional[str] = None,
    pre_impute_expression_matrix: Optional[pd.DataFrame] = None,
    gene_symbols: Optional[Mapping[Any, Any]] = None,
    gene_descriptions: Optional[Mapping[str, str]] = None,
    funcats: Optional[Mapping[Any, Any]] = None,
    group: Optional[str] = None,
    formula_in: Optional[str] = None,
    formula_effective: Optional[str] = None,
    formula_rewritten: Optional[str] = None,
    contrasts_spec: Optional[str] = None,
    block: Optional[str] = None,
    limma_robust: bool = True,
    limma_trend: bool = True,
    impute_missing_values: bool = False,
    fill_na_zero: bool = False,
    normtype: Optional[str] = None,
    non_zeros: Optional[Any] = None,
    taxon: Optional[str] = None,
    batch_applied: Optional[bool] = None,
    batch_method: Optional[str] = None,
    imputation_backend: Optional[str] = None,
    gaussian_method: Optional[str] = None,
    lupine_mode: Optional[str] = None,
    force: bool = False,
) -> LimmaReplayFiles:
    analysis_root = Path(analysis_dir).expanduser().resolve()
    volcano_root = Path(volcano_dir).expanduser().resolve()

    if expression_matrix is not None:
        edata = expression_matrix.copy()
        matrix_source = expression_matrix_source or "provided_expression_matrix"
    else:
        sample_ids = list(sample_metadata.index.astype(str))
        edata = _extract_limma_input_from_results(results, sample_ids=sample_ids)
        matrix_source = expression_matrix_source or "limma_result_table_sample_columns"

    # Make a replay-specific cdesc that matches the matrix columns (and include injected covariates if present).
    cdesc = sample_metadata.copy()
    if list(cdesc.index.astype(str)) != list(edata.columns.astype(str)):
        cdesc = cdesc.reindex(list(edata.columns.astype(str)))
    if "id" not in cdesc.columns:
        cdesc["id"] = list(cdesc.index.astype(str))

    # Build rdesc (optional but helpful for downstream tables).
    rdesc = _make_rdesc(
        list(edata.index),
        gene_symbols=gene_symbols,
        gene_descriptions=gene_descriptions,
        funcats=funcats,
    )

    # Parse resolved coefficient/contrast lists from result labels.
    direct_coefs, contrast_exprs = parse_limma_result_labels(results.keys())

    pre_edata = None
    pre_impute_shape: Optional[List[int]] = None
    if pre_impute_expression_matrix is not None:
        pre_edata = pre_impute_expression_matrix.copy()
        pre_edata = pre_edata.reindex(index=edata.index, columns=edata.columns)
        pre_impute_shape = [int(pre_edata.shape[0]), int(pre_edata.shape[1])]

    matrix_hash = _hash_frame(edata)
    pre_impute_matrix_hash = _hash_frame(pre_edata) if pre_edata is not None else None
    cdesc_hash = _hash_frame(cdesc)
    rdesc_hash = _hash_frame(rdesc)

    # Stable run id derived from the resolved replay inputs (not timestamps).
    run_params = {
        "replay_contract_version": REPLAY_CONTRACT_VERSION,
        "matrix_storage_format": "gctx",
        "matrix_storage_dtype": REPLAY_MATRIX_DTYPE,
        "gctx_writer_contract": gct_io_contracts()["gctx"]["writer_contract"],
        "matrix_hash": matrix_hash,
        "frame_hash_algorithm": selected_content_hash_algorithm(),
        "pre_impute_matrix_hash": pre_impute_matrix_hash,
        "cdesc_hash": cdesc_hash,
        "rdesc_hash": rdesc_hash,
        "group": group,
        "formula_effective": formula_effective,
        "formula_rewritten": formula_rewritten,
        "contrasts_spec": contrasts_spec,
        "block": block,
        "limma_robust": bool(limma_robust),
        "limma_trend": bool(limma_trend),
        "impute_missing_values": bool(impute_missing_values),
        "fill_na_zero": bool(fill_na_zero),
        "imputation_backend": imputation_backend,
        "gaussian_method": gaussian_method if imputation_backend == "gaussian" else None,
        "lupine_mode": lupine_mode if imputation_backend == "lupine" else None,
        "direct_coefs": direct_coefs,
        "contrast_expressions": contrast_exprs,
    }
    run_id = _compute_run_id(run_params)

    replay_dir = volcano_root / "replay" / run_id
    replay_dir.mkdir(parents=True, exist_ok=True)

    gct_path = _write_gct_reuse_existing(
        out_path=replay_dir / "limma_input.gctx",
        mat=edata,
        cdesc=cdesc,
        rdesc=rdesc,
        precision=REPLAY_GCT_PRECISION,
    )
    pre_impute_gct_path: Optional[Path] = None
    if pre_edata is not None:
        pre_impute_gct_path = _write_gct_reuse_existing(
            out_path=replay_dir / "limma_input_pre_impute.gctx",
            mat=pre_edata,
            cdesc=cdesc,
            rdesc=rdesc,
            precision=REPLAY_GCT_PRECISION,
        )
    replay_rmd_path = replay_dir / "replay_explore.Rmd"
    replay_render_path = replay_dir / "render_replay_explore.sh"
    replay_html_path = replay_dir / "replay_explore.html"

    now = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    context = {
        "generated_at": now,
        "replay_contract_version": REPLAY_CONTRACT_VERSION,
        "replay_scope": "limma_only",
        "path_base": "replay_dir",
        "analysis_dir": _relative_path(analysis_root, replay_dir),
        "volcano_dir": _relative_path(volcano_root, replay_dir),
        "replay_dir": ".",
        "run_id": run_id,
        "matrix_shape": [int(edata.shape[0]), int(edata.shape[1])],
        "matrix_hash": matrix_hash,
        "frame_hash_algorithm": selected_content_hash_algorithm(),
        "expression_matrix_source": matrix_source,
        "pre_impute_matrix_hash": pre_impute_matrix_hash,
        "cdesc_hash": cdesc_hash,
        "rdesc_hash": rdesc_hash,
        "stored_matrix_role": (
            "limma_input_imputed" if bool(impute_missing_values) else "limma_input"
        ),
        "stored_matrix_is_authoritative": True,
        "matrix_storage_format": "gctx",
        "matrix_storage_dtype": REPLAY_MATRIX_DTYPE,
        "gctx_writer_contract": gct_io_contracts()["gctx"]["writer_contract"],
        "gctx_content_hash_algorithm": selected_content_hash_algorithm(),
        "gctx_content_hash": read_gctx_content_hash(gct_path),
        "pre_impute_gctx_content_hash": (
            read_gctx_content_hash(pre_impute_gct_path)
            if pre_impute_gct_path is not None
            else None
        ),
        "has_pre_impute_matrix": pre_impute_gct_path is not None,
        "pre_impute_matrix_role": (
            "audit_only" if pre_impute_gct_path is not None else None
        ),
        "pre_impute_matrix_shape": pre_impute_shape,
        "group": group,
        "formula_in": formula_in,
        "formula_effective": formula_effective,
        "formula_rewritten": formula_rewritten,
        "contrasts_spec": contrasts_spec,
        "block": block,
        "limma_robust": bool(limma_robust),
        "limma_trend": bool(limma_trend),
        "impute_missing_values": bool(impute_missing_values),
        "fill_na_zero": bool(fill_na_zero),
        "imputation_backend": imputation_backend,
        "gaussian_method": gaussian_method if imputation_backend == "gaussian" else None,
        "lupine_mode": lupine_mode if imputation_backend == "lupine" else None,
        "normtype": normtype,
        "non_zeros": non_zeros,
        "taxon": taxon,
        "batch_applied": batch_applied,
        "batch_method": batch_method,
        "gct_path": _relative_path(gct_path, replay_dir),
        "pre_impute_gct_path": (
            _relative_path(pre_impute_gct_path, replay_dir)
            if pre_impute_gct_path
            else None
        ),
        "recompute_imputation_supported": False,
        "recompute_imputation_note": (
            "Disabled by design: limma replay always consumes the stored modeled matrix and never re-imputes."
        ),
        "pca_included": False,
        "pca_replay_note": (
            "PCA is outside this limma replay contract and requires its own saved input and pca2 options."
        ),
        "replay_explore_rmd_path": _relative_path(replay_rmd_path, replay_dir),
        "replay_explore_render_path": _relative_path(replay_render_path, replay_dir),
        "replay_explore_html_path": _relative_path(replay_html_path, replay_dir),
        "result_labels": list(results.keys()),
        "direct_coefs": direct_coefs,
        "contrast_expressions": contrast_exprs,
    }
    context_path = replay_dir / "limma_replay_context.json"
    _write_json(context_path, context, force=force)
    _write_text(
        replay_rmd_path,
        render_limma_replay_rmd(
            title="Tackle limma replay",
            generated_at=now,
            gct_relpath=gct_path.name,
            context_relpath=context_path.name,
            outputs_rel_dir="replay_outputs",
        ),
        force=force,
    )
    _write_text(
        replay_render_path,
        render_limma_replay_sh(
            rmd_name=replay_rmd_path.name,
            html_name=replay_html_path.name,
        ),
        force=force,
    )

    # Update "last run" pointer under analysis context.
    pointer_path = analysis_root / "context" / "last_volcano_replay.json"
    pointer_payload = {
        "updated_at": now,
        "path_base": "analysis_dir",
        "analysis_dir": ".",
        "volcano_dir": _relative_path(volcano_root, analysis_root),
        "replay_dir": _relative_path(replay_dir, analysis_root),
        "run_id": run_id,
        "gct_path": _relative_path(gct_path, analysis_root),
        "pre_impute_gct_path": (
            _relative_path(pre_impute_gct_path, analysis_root)
            if pre_impute_gct_path
            else None
        ),
        "context_path": _relative_path(context_path, analysis_root),
    }
    _write_json(pointer_path, pointer_payload, force=True)

    logger.info("Limma replay: wrote %s and %s", gct_path, context_path)
    return LimmaReplayFiles(
        volcano_dir=volcano_root,
        replay_dir=replay_dir,
        run_id=run_id,
        gct_path=gct_path,
        context_path=context_path,
        pointer_path=pointer_path,
        pre_impute_gct_path=pre_impute_gct_path,
    )


def discover_replay_contexts(base_dir: Path) -> List[Path]:
    root = Path(base_dir).expanduser().resolve()
    return sorted(root.rglob("limma_replay_context.json"))


def resolve_replay_dir(
    *,
    analysis_dir: Path,
    volcano_dir: Optional[Path] = None,
    run_id: Optional[str] = None,
) -> Path:
    analysis_dir = Path(analysis_dir).expanduser().resolve()

    if volcano_dir is not None:
        vdir = Path(volcano_dir).expanduser().resolve()
        if (vdir / "limma_replay_context.json").exists():
            return vdir
        replay_root = vdir / "replay"
        if replay_root.exists():
            candidates = [p for p in replay_root.iterdir() if p.is_dir()]
            candidates.sort(key=lambda p: p.name)
            if run_id:
                match = replay_root / str(run_id)
                if match.exists() and match.is_dir():
                    return match
                raise FileNotFoundError(f"Replay run_id not found under {replay_root}: {run_id}")
            if len(candidates) == 1:
                return candidates[0]
            if not candidates:
                raise FileNotFoundError(f"No replay runs found under: {replay_root}")
            raise RuntimeError(
                "Multiple replay runs found; pass --run-id to disambiguate: "
                + ", ".join(p.name for p in candidates[:12])
            )
        raise FileNotFoundError(f"volcano_dir does not contain replay data: {vdir}")

    pointer = analysis_dir / "context" / "last_volcano_replay.json"
    if pointer.exists():
        payload = json.loads(pointer.read_text(encoding="utf-8"))
        replay_dir = Path(payload.get("replay_dir", "")).expanduser()
        if not replay_dir.is_absolute():
            replay_dir = (analysis_dir / replay_dir).resolve()
        if replay_dir.exists():
            return replay_dir

    contexts = discover_replay_contexts(analysis_dir)
    if len(contexts) == 1:
        return contexts[0].parent
    if not contexts:
        raise FileNotFoundError(
            "No limma replay contexts found. Run a limma volcano first (it should write replay inputs under volcano/**/replay/)."
        )
    rels = [str(p.relative_to(analysis_dir)) for p in contexts[:20]]
    raise RuntimeError(
        "Multiple limma replay contexts found; pass --volcano-dir to select one.\n- "
        + "\n- ".join(rels)
    )


def validate_replay_dir(replay_dir: Path) -> Tuple[Path, Path]:
    replay_dir = Path(replay_dir).expanduser().resolve()
    context_path = replay_dir / "limma_replay_context.json"
    if not context_path.exists():
        raise FileNotFoundError(str(context_path))

    # The context records the exact content-addressed GCTX (or legacy GCT), so
    # it is the authoritative selector. This also prevents the pre-imputation
    # audit matrix from winning a broad fallback glob.
    context = json.loads(context_path.read_text(encoding="utf-8"))
    context_gct = context.get("gct_path")
    gct_path: Optional[Path] = None
    if context_gct:
        candidate = Path(str(context_gct)).expanduser()
        if not candidate.is_absolute():
            candidate = replay_dir / candidate
        candidate = candidate.resolve()
        if (
            candidate.exists()
            and candidate.is_file()
            and not candidate.name.startswith("limma_input_pre_impute")
        ):
            gct_path = candidate

    if gct_path is None:
        candidates = []
        for suffix in (".gctx", ".gct"):
            exact = replay_dir / f"limma_input{suffix}"
            if exact.exists():
                candidates.append(exact)
            candidates.extend(
                path
                for path in replay_dir.glob(f"limma_input*{suffix}")
                if not path.name.startswith("limma_input_pre_impute")
            )
        candidates = [path for path in set(candidates) if path.is_file()]
        if candidates:
            gct_path = max(candidates, key=lambda p: p.stat().st_mtime)
        else:
            raise FileNotFoundError(str(replay_dir / "limma_input.gctx"))
    return gct_path, context_path


def sanitize_for_filename(text: str) -> str:
    s = str(text)
    s = re.sub(r"\s+", "_", s.strip())
    s = re.sub(r"[^\w.\-]+", "_", s)
    s = s.strip("._-")
    return s or "limma"
