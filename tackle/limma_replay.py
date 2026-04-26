from __future__ import annotations

import hashlib
import json
import os
import re
import textwrap
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from .utils import _get_logger

logger = _get_logger(__name__)


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


def _compute_run_id(params: Mapping[str, Any]) -> str:
    payload = json.dumps(params, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:10]


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
    precision: int = 4,
) -> Path:
    """Write a .gct using cmapR via rpy2; returns the written .gct path."""
    from rpy2 import robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    cmapR = importr("cmapR")
    _ = cmapR  # quiet lint

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    base = str(out_path)
    base_lower = base.lower()
    for suffix in (".gct.gz", ".gct"):
        if base_lower.endswith(suffix):
            base = base[: -len(suffix)]
            break
    out_base = Path(base)

    # Ensure aligned ordering.
    cid = [str(c) for c in mat.columns]
    rid = [str(r) for r in mat.index]
    cdesc = cdesc.copy()
    if "id" not in cdesc.columns:
        cdesc["id"] = [str(x) for x in cdesc.index]
    if list(cdesc.index.astype(str)) != cid:
        cdesc = cdesc.reindex(cid)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        robjects.r.assign("m", mat.astype(float))
        robjects.r.assign("rid", rid)
        robjects.r.assign("cid", cid)
        robjects.r.assign("cdesc", cdesc)
        robjects.r.assign("rdesc", rdesc)
        robjects.r.assign("outname", os.path.normpath(os.path.abspath(str(out_base))))

    robjects.r(
        """
my_ds <- new("GCT",
             mat=as.matrix(m),
             rid=rid,
             cid=cid,
             cdesc=cdesc,
             rdesc=as.data.frame(rdesc))
write_gct(my_ds, outname, precision=as.integer(%d))
"""
        % int(precision)
    )

    written = out_base.with_suffix(".gct")
    if written.exists():
        return written
    # cmapR may append timestamps; best-effort find newest candidate.
    candidates = list(out_base.parent.glob(out_base.name + "*.gct"))
    if candidates:
        return max(candidates, key=lambda p: p.stat().st_mtime)
    raise FileNotFoundError(f"Failed to write GCT under: {out_base.parent}")


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


def _render_replay_explore_rmd(
    *,
    gct_relpath: str,
    context_relpath: str,
    pre_impute_gct_relpath: Optional[str] = None,
) -> str:
    pre_impute_gct_literal = json.dumps(pre_impute_gct_relpath) if pre_impute_gct_relpath else "NULL"
    return textwrap.dedent(
        f"""\
        ---
        title: "Volcano Replay Explore"
        output:
          html_document:
            toc: true
            toc_float: true
        ---

        ```{{r setup, include=FALSE}}
        knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
        ```

        ```{{r load}}
        if (!requireNamespace("jsonlite", quietly = TRUE)) {{
          stop("Missing dependency: jsonlite")
        }}
        if (!requireNamespace("cmapR", quietly = TRUE)) {{
          stop("Missing dependency: cmapR")
        }}
        ctx <- jsonlite::read_json("{context_relpath}", simplifyVector = TRUE)
        suppressPackageStartupMessages(library(cmapR))

        parse_gct_any <- function(path) {{
          exports <- getNamespaceExports("cmapR")
          if ("parse.gctx" %in% exports) return(cmapR::parse.gctx(path))
          if ("parse_gctx" %in% exports) return(cmapR::parse_gctx(path))
          stop("cmapR does not export parse.gctx / parse_gctx")
        }}

        stored_ds <- parse_gct_any("{gct_relpath}")
        mat_stored <- stored_ds@mat
        cdesc <- stored_ds@cdesc
        rdesc <- stored_ds@rdesc
        mat <- mat_stored

        pre_impute_gct_path <- {pre_impute_gct_literal}
        mat_pre_impute <- NULL
        if (!is.null(pre_impute_gct_path) && nzchar(pre_impute_gct_path) && file.exists(pre_impute_gct_path)) {{
          pre_ds <- parse_gct_any(pre_impute_gct_path)
          mat_pre_impute <- pre_ds@mat
        }}
        ```

        ## Replay Context

        ```{{r context}}
        key_fields <- c(
          "run_id", "matrix_shape", "stored_matrix_role",
          "group", "formula_in", "formula_effective", "formula_rewritten",
          "contrasts_spec", "block",
          "has_pre_impute_matrix", "recompute_imputation_supported",
          "impute_missing_values", "imputation_backend", "gaussian_method", "lupine_mode",
          "limma_robust", "limma_trend",
          "normtype", "non_zeros", "taxon", "batch_applied", "batch_method"
        )
        present <- intersect(key_fields, names(ctx))
        as.data.frame(t(unlist(ctx[present])), stringsAsFactors = FALSE)
        ```

        ## Stored Matrix Replay

        ```{{r stored-matrix-replay}}
        cat("The stored limma_input.gct matrix is the authoritative replay matrix.\\n")
        cat("All summary plots below use mat_stored as written by tackle.\\n")
        if (!is.null(mat_pre_impute)) {{
          cat("A pre-imputation matrix is also available at: ", pre_impute_gct_path, "\\n", sep = "")
          cat("Pre-impute matrix dimensions: ", paste(dim(mat_pre_impute), collapse = " x "), "\\n", sep = "")
        }} else {{
          cat("No pre-imputation matrix was exported for this replay bundle.\\n")
        }}
        ```

        ## Optional Gaussian Recompute

        ```{{r optional-gaussian-recompute}}
        deterministic_gaussian_impute <- function(x, downshift = 1.8, width = 0.3, seed = 1234) {{
          out <- as.matrix(x)
          out[out == 0] <- NA_real_
          obs <- as.numeric(out[is.finite(out)])
          miss <- is.na(out)
          if (!any(miss) || length(obs) < 2) return(out)
          obs_sd <- stats::sd(obs)
          if (!is.finite(obs_sd) || obs_sd <= 0) obs_sd <- 1
          set.seed(seed)
          out[miss] <- stats::rnorm(
            sum(miss),
            mean = mean(obs) - downshift * obs_sd,
            sd = max(.Machine$double.eps, width * obs_sd)
          )
          out
        }}

        can_recompute <- isTRUE(ctx$impute_missing_values) &&
          identical(ctx$imputation_backend, "gaussian") &&
          !is.null(mat_pre_impute)

        if (can_recompute) {{
          method <- ctx$gaussian_method
          if (is.null(method) || !nzchar(method)) method <- "legacy"
          width <- if (identical(method, "mqish")) 0.3 else 0.8
          mat_recomputed <- deterministic_gaussian_impute(
            mat_pre_impute,
            downshift = 1.8,
            width = width,
            seed = 1234
          )
          common_rows <- intersect(rownames(mat_recomputed), rownames(mat_stored))
          common_cols <- intersect(colnames(mat_recomputed), colnames(mat_stored))
          delta <- mat_recomputed[common_rows, common_cols, drop = FALSE] -
            mat_stored[common_rows, common_cols, drop = FALSE]
          finite_delta <- as.numeric(delta[is.finite(delta)])
          cat("R-side deterministic Gaussian recompute completed for inspection.\\n")
          cat("Exact replay should use mat_stored; this chunk is a reproducibility diagnostic.\\n")
          if (length(finite_delta) > 0) {{
            knitr::kable(data.frame(
              metric = c("mean_abs_delta", "max_abs_delta", "correlation"),
              value = c(
                mean(abs(finite_delta)),
                max(abs(finite_delta)),
                suppressWarnings(stats::cor(
                  as.numeric(mat_recomputed[common_rows, common_cols]),
                  as.numeric(mat_stored[common_rows, common_cols]),
                  use = "complete.obs"
                ))
              )
            ))
          }} else {{
            cat("No finite overlapping values were available for comparison.\\n")
          }}
        }} else {{
          cat("Gaussian recompute skipped. This requires imputation_backend == 'gaussian' and limma_input_pre_impute.gct.\\n")
          if (identical(ctx$imputation_backend, "lupine")) {{
            cat("Lupine replay uses the stored matrix; model-side recompute is intentionally not attempted here.\\n")
          }}
        }}
        ```

        ## Matrix Summary

        ```{{r matrix-summary}}
        summary(as.vector(mat))
        sample_summary <- data.frame(
          sample = colnames(mat),
          mean = colMeans(mat, na.rm = TRUE),
          median = apply(mat, 2, median, na.rm = TRUE),
          sd = apply(mat, 2, sd, na.rm = TRUE),
          na_count = colSums(is.na(mat)),
          stringsAsFactors = FALSE
        )
        knitr::kable(sample_summary)
        ```

        ## Global Distribution

        ```{{r global-distribution}}
        hist(as.vector(mat), breaks = 80, main = "Global expression distribution", xlab = "Expression")
        lines(density(as.vector(mat), na.rm = TRUE), col = "steelblue", lwd = 2)
        ```

        ## Per-sample Distributions

        ```{{r sample-distributions, fig.width=10, fig.height=5}}
        boxplot(as.data.frame(mat), las = 2, outline = FALSE, main = "Per-sample distributions", ylab = "Expression")
        ```

        ## Sample Correlation

        ```{{r sample-correlation, fig.width=7, fig.height=7}}
        corr <- cor(mat, use = "pairwise.complete.obs")
        heatmap(corr, symm = TRUE, margins = c(8, 8), main = "Sample correlation")
        ```

        ## PCA

        ```{{r pca}}
        complete_rows <- complete.cases(mat)
        if (sum(complete_rows) >= 2 && ncol(mat) >= 2) {{
          mat_complete <- mat[complete_rows, , drop = FALSE]
          pcs <- prcomp(t(mat_complete), center = TRUE, scale. = TRUE)
          pct <- round(100 * summary(pcs)$importance[2, 1:2], 1)
          plot(
            pcs$x[, 1], pcs$x[, 2],
            pch = 19,
            xlab = paste0("PC1 (", pct[[1]], "%)"),
            ylab = paste0("PC2 (", pct[[2]], "%)"),
            main = "Sample PCA"
          )
          text(pcs$x[, 1], pcs$x[, 2], labels = rownames(pcs$x), pos = 3, cex = 0.7)
        }} else {{
          cat("Not enough complete rows to compute PCA.\\n")
        }}
        ```

        ## Top Variable Genes

        ```{{r top-variable-genes}}
        rv <- apply(mat, 1, var, na.rm = TRUE)
        rv <- rv[is.finite(rv)]
        if (length(rv) > 0) {{
          top_ids <- names(sort(rv, decreasing = TRUE))[seq_len(min(20, length(rv)))]
          gene_symbol <- if ("GeneSymbol" %in% colnames(rdesc)) rdesc[top_ids, "GeneSymbol"] else rep("", length(top_ids))
          out <- data.frame(
            GeneID = top_ids,
            GeneSymbol = gene_symbol,
            Variance = unname(rv[top_ids]),
            stringsAsFactors = FALSE
          )
          knitr::kable(out)
        }} else {{
          cat("No finite row variances were available.\\n")
        }}
        ```
        """
    )


def _render_replay_explore_sh(*, rmd_name: str, html_name: str) -> str:
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env bash
        set -euo pipefail
        cd "$(dirname "$0")"
        Rscript -e "rmarkdown::render('{rmd_name}', output_file = '{html_name}')"
        """
    )


def write_limma_replay_files(
    *,
    analysis_dir: str,
    volcano_dir: str,
    results: Mapping[str, pd.DataFrame],
    sample_metadata: pd.DataFrame,
    expression_matrix: Optional[pd.DataFrame] = None,
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
    else:
        sample_ids = list(sample_metadata.index.astype(str))
        edata = _extract_limma_input_from_results(results, sample_ids=sample_ids)

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

    # Stable run id derived from the resolved replay inputs (not timestamps).
    run_params = {
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

    gct_path = _write_gct(out_path=replay_dir / "limma_input.gct", mat=edata, cdesc=cdesc, rdesc=rdesc)
    pre_impute_gct_path: Optional[Path] = None
    pre_impute_shape: Optional[List[int]] = None
    if pre_impute_expression_matrix is not None:
        pre_edata = pre_impute_expression_matrix.copy()
        pre_edata = pre_edata.reindex(index=edata.index, columns=edata.columns)
        pre_impute_shape = [int(pre_edata.shape[0]), int(pre_edata.shape[1])]
        pre_impute_gct_path = _write_gct(
            out_path=replay_dir / "limma_input_pre_impute.gct",
            mat=pre_edata,
            cdesc=cdesc,
            rdesc=rdesc,
        )
    replay_rmd_path = replay_dir / "replay_explore.Rmd"
    replay_render_path = replay_dir / "render_replay_explore.sh"
    replay_html_path = replay_dir / "replay_explore.html"

    now = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    context = {
        "generated_at": now,
        "analysis_dir": str(analysis_root),
        "volcano_dir": str(volcano_root),
        "replay_dir": str(replay_dir),
        "run_id": run_id,
        "matrix_shape": [int(edata.shape[0]), int(edata.shape[1])],
        "stored_matrix_role": (
            "limma_input_imputed" if bool(impute_missing_values) else "limma_input"
        ),
        "stored_matrix_is_authoritative": True,
        "has_pre_impute_matrix": pre_impute_gct_path is not None,
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
        "gct_path": str(gct_path),
        "pre_impute_gct_path": str(pre_impute_gct_path) if pre_impute_gct_path else None,
        "recompute_imputation_supported": bool(
            impute_missing_values
            and imputation_backend == "gaussian"
            and pre_impute_gct_path is not None
        ),
        "recompute_imputation_note": (
            "R-side Gaussian recompute is a diagnostic; exact replay uses the stored limma_input.gct matrix."
            if imputation_backend == "gaussian"
            else "Stored matrix replay is authoritative; non-Gaussian backends are not recomputed in the Rmd."
        ),
        "replay_explore_rmd_path": str(replay_rmd_path),
        "replay_explore_render_path": str(replay_render_path),
        "replay_explore_html_path": str(replay_html_path),
        "result_labels": list(results.keys()),
        "direct_coefs": direct_coefs,
        "contrast_expressions": contrast_exprs,
    }
    context_path = replay_dir / "limma_replay_context.json"
    _write_json(context_path, context, force=force)
    _write_text(
        replay_rmd_path,
        _render_replay_explore_rmd(
            gct_relpath=gct_path.name,
            context_relpath=context_path.name,
            pre_impute_gct_relpath=pre_impute_gct_path.name if pre_impute_gct_path else None,
        ),
        force=force,
    )
    _write_text(
        replay_render_path,
        _render_replay_explore_sh(
            rmd_name=replay_rmd_path.name,
            html_name=replay_html_path.name,
        ),
        force=force,
    )

    # Update "last run" pointer under analysis context.
    pointer_path = analysis_root / "context" / "last_volcano_replay.json"
    pointer_payload = {
        "updated_at": now,
        "analysis_dir": str(analysis_root),
        "volcano_dir": str(volcano_root),
        "replay_dir": str(replay_dir),
        "run_id": run_id,
        "gct_path": str(gct_path),
        "pre_impute_gct_path": str(pre_impute_gct_path) if pre_impute_gct_path else None,
        "context_path": str(context_path),
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

    # GCT may be written as limma_input.gct (expected), or with a timestamp; prefer exact.
    gct_path = replay_dir / "limma_input.gct"
    if not gct_path.exists():
        candidates = sorted(replay_dir.glob("limma_input*.gct"))
        if candidates:
            gct_path = max(candidates, key=lambda p: p.stat().st_mtime)
        else:
            raise FileNotFoundError(str(gct_path))
    return gct_path, context_path


def sanitize_for_filename(text: str) -> str:
    s = str(text)
    s = re.sub(r"\s+", "_", s.strip())
    s = re.sub(r"[^\w.\-]+", "_", s)
    s = s.strip("._-")
    return s or "limma"
