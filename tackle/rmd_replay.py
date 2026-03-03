from __future__ import annotations

import shutil
import subprocess
import textwrap
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from .limma_replay import validate_replay_dir


@dataclass(frozen=True)
class RmdReplayBundle:
    report_dir: Path
    rmd_path: Path
    html_path: Optional[Path]
    gct_path: Path
    context_path: Path
    render_sh: Path


def _safe_write_text(path: Path, content: str, *, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"{path} already exists (use --force to overwrite)")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _copy_or_symlink(src: Path, dst: Path, *, copy: bool, force: bool) -> None:
    if (dst.exists() or dst.is_symlink()) and not force:
        raise FileExistsError(f"{dst} already exists (use --force to overwrite)")
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.parent.mkdir(parents=True, exist_ok=True)
    if copy:
        shutil.copy2(src, dst)
    else:
        dst.symlink_to(src.resolve())


def _render_limma_replay_rmd(
    *,
    title: str,
    generated_at: str,
    gct_relpath: str,
    context_relpath: str,
    outputs_rel_dir: str = "rmd_outputs",
) -> str:
    return textwrap.dedent(
        f"""\
        ---
        title: "{title}"
        date: "{generated_at}"
        output:
          html_document:
            toc: true
            toc_float: true
            df_print: paged
        ---

        ```{{r setup, include=FALSE}}
        knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
        ```

        ## Run Summary

        - Generated: `{generated_at}`
        - Inputs:
          - GCT: `{gct_relpath}`
          - Context: `{context_relpath}`

        ```{{r load-context}}
        if (!requireNamespace("jsonlite", quietly = TRUE)) {{
          stop("Missing dependency: jsonlite")
        }}
        ctx <- jsonlite::read_json("{context_relpath}", simplifyVector = TRUE)
        key_fields <- c(
          "analysis_dir", "volcano_dir", "run_id",
          "group", "formula_in", "formula_effective", "formula_rewritten",
          "contrasts_spec", "block",
          "limma_robust", "limma_trend",
          "impute_missing_values", "fill_na_zero",
          "normtype", "non_zeros", "taxon", "batch_applied", "batch_method"
        )
        present <- intersect(key_fields, names(ctx))
        as.data.frame(t(unlist(ctx[present])), stringsAsFactors = FALSE)
        ```

        ## Load Data

        ```{{r load-gct}}
        if (!requireNamespace("cmapR", quietly = TRUE)) {{
          stop("Missing dependency: cmapR")
        }}
        suppressPackageStartupMessages(library(cmapR))

        parse_gct_any <- function(path) {{
          exports <- getNamespaceExports("cmapR")
          if ("parse.gctx" %in% exports) return(cmapR::parse.gctx(path))
          if ("parse_gctx" %in% exports) return(cmapR::parse_gctx(path))
          stop("cmapR does not export parse.gctx / parse_gctx")
        }}

        ds <- parse_gct_any("{gct_relpath}")
        mat <- ds@mat
        cdesc <- ds@cdesc
        rdesc <- ds@rdesc

        if ("id" %in% colnames(cdesc)) {{
          rownames(cdesc) <- cdesc$id
        }}
        cdesc <- cdesc[colnames(mat), , drop = FALSE]
        stopifnot(all(rownames(cdesc) == colnames(mat)))

        if ("GeneID" %in% colnames(rdesc)) {{
          rownames(rdesc) <- rdesc$GeneID
        }} else if ("id" %in% colnames(rdesc)) {{
          rownames(rdesc) <- rdesc$id
        }}

        dim(mat)
        head(cdesc)
        head(rdesc)
        ```

        ## Input Validation

        ```{{r validate}}
        # Fail loudly if required design fields are missing.
        group_col <- if (!is.null(ctx$group) && nzchar(ctx$group)) ctx$group else NULL
        formula_str <- NULL
        if (!is.null(ctx$formula_rewritten) && nzchar(ctx$formula_rewritten)) {{
          formula_str <- ctx$formula_rewritten
        }} else if (!is.null(ctx$formula_effective) && nzchar(ctx$formula_effective)) {{
          formula_str <- ctx$formula_effective
        }}

        if (is.null(formula_str) && is.null(group_col)) {{
          stop("Context is missing both formula and group.")
        }}
        if (!is.null(group_col) && !(group_col %in% colnames(cdesc))) {{
          stop(paste0("Missing group column in cdesc: ", group_col))
        }}
        if (!is.null(ctx$block) && nzchar(ctx$block) && !(ctx$block %in% colnames(cdesc))) {{
          stop(paste0("Missing block column in cdesc: ", ctx$block))
        }}
        if (!is.null(group_col)) {{
          vc <- table(as.factor(cdesc[[group_col]]), useNA = "ifany")
          vc
        }}
        ```

        ## Preprocessing

        This report replays limma using the matrix stored in the input GCT (already normalized/filtered/transformed by the pipeline).

        ## Design + Contrasts

        ```{{r design-contrasts}}
        if (!requireNamespace("limma", quietly = TRUE)) {{
          stop("Missing dependency: limma")
        }}
        suppressPackageStartupMessages(library(limma))

        pheno <- as.data.frame(cdesc, stringsAsFactors = FALSE)

        sanitize_name <- function(x) {{
          x <- gsub(":", "_", x, fixed = TRUE)
          x <- gsub(" ", "_", x, fixed = TRUE)
          x <- gsub("-", "_", x, fixed = TRUE)
          x <- gsub("\\\\+", "_", x, fixed = FALSE)
          x <- gsub("/", "_", x, fixed = TRUE)
          x <- gsub("\\\\?", "qmk", x, fixed = FALSE)
          x
        }}

        if (is.null(formula_str)) {{
          mod <- model.matrix(as.formula(paste0("~0+", group_col)), pheno)
        }} else {{
          mod <- model.matrix(as.formula(formula_str), pheno)
        }}
        colnames(mod) <- make.names(sanitize_name(colnames(mod)))

        mod

        contrast_exprs <- ctx$contrast_expressions
        if (is.null(contrast_exprs)) contrast_exprs <- character()
        if (!is.character(contrast_exprs)) contrast_exprs <- as.character(contrast_exprs)

        contrasts_matrix <- NULL
        if (length(contrast_exprs) > 0) {{
          contrasts_matrix <- makeContrasts(contrasts = contrast_exprs, levels = mod)
          contrasts_matrix
        }}
        ```

        ## Limma Fit

        ```{{r fit}}
        robust_flag <- isTRUE(ctx$limma_robust)
        trend_flag <- isTRUE(ctx$limma_trend)

        block <- NULL
        cor <- NULL
        if (!is.null(ctx$block) && nzchar(ctx$block)) {{
          block <- as.factor(pheno[[ctx$block]])
          corfit <- duplicateCorrelation(mat, design = mod, block = block)
          cor <- corfit$consensus
          cor
        }}

        fit <- lmFit(mat, mod, block = block, correlation = cor)
        fit2_base <- eBayes(fit, robust = robust_flag, trend = trend_flag)

        fit2_con <- NULL
        if (!is.null(contrasts_matrix)) {{
          fit2_con <- contrasts.fit(fit, contrasts_matrix)
          fit2_con <- eBayes(fit2_con, robust = robust_flag, trend = trend_flag)
        }}
        ```

        ## Results

        ```{{r results}}
        out_dir <- "{outputs_rel_dir}"
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        safe_file <- function(x) {{
          x <- gsub("\\\\s+", "_", x)
          x <- gsub("[^A-Za-z0-9._-]+", "_", x)
          x <- gsub("^[_\\\\.-]+", "", x)
          x <- gsub("[_\\\\.-]+$", "", x)
          if (!nzchar(x)) x <- "limma"
          x
        }}

        add_gene_meta <- function(tt) {{
          tt$GeneID <- rownames(tt)
          tt <- cbind(tt, rdesc[tt$GeneID, setdiff(colnames(rdesc), c("GeneID", "id")), drop = FALSE])
          tt
        }}

        to_log2 <- function(x) {{
          x / log10(2)
        }}

        normalize_tt <- function(tt) {{
          # Match tackle output naming where possible.
          if ("logFC" %in% colnames(tt)) tt$log2_FC <- to_log2(tt$logFC)
          if ("CI.L" %in% colnames(tt)) tt$CI.L <- to_log2(tt$CI.L)
          if ("CI.R" %in% colnames(tt)) tt$CI.R <- to_log2(tt$CI.R)
          if ("adj.P.Val" %in% colnames(tt)) tt$pAdj <- tt$adj.P.Val
          if ("P.Value" %in% colnames(tt)) tt$pValue <- tt$P.Value
          tt
        }}

        labels <- ctx$result_labels
        if (is.null(labels)) {{
          labels <- c()
        }}
        if (!is.character(labels)) labels <- as.character(labels)

        tables <- list()
        for (lab in labels) {{
          parts <- strsplit(lab, "=", fixed = TRUE)[[1]]
          left <- parts[1]
          expr <- if (length(parts) >= 2) parts[2] else left

          if (startsWith(left, "coef_")) {{
            tt <- topTable(fit2_base, n = Inf, sort.by = "none", coef = expr, confint = TRUE)
          }} else {{
            if (is.null(fit2_con)) stop("Contrast results requested but contrast fit is NULL.")
            tt <- topTable(fit2_con, n = Inf, sort.by = "none", coef = expr, confint = TRUE)
          }}
          tt <- add_gene_meta(tt)
          tt <- normalize_tt(tt)
          out_path <- file.path(out_dir, paste0("limma_", safe_file(left), ".tsv"))
          write.table(tt, out_path, sep = "\\\\t", quote = FALSE, row.names = FALSE)
          tables[[lab]] <- tt
        }}

        # Inline preview of the first result (if present).
        if (length(tables) > 0) {{
          head(tables[[1]], 20)
        }} else {{
          "No result labels found in context."
        }}
        ```

        ## Diagnostic Plots

        ```{{r diagnostic-plots}}
        # Write one MD plot for the first available result label.
        if (length(labels) > 0) {{
          lab <- labels[[1]]
          parts <- strsplit(lab, "=", fixed = TRUE)[[1]]
          left <- parts[1]
          expr <- if (length(parts) >= 2) parts[2] else left
          png(file.path(out_dir, paste0("md_plot_", safe_file(left), ".png")), width = 1200, height = 900, res = 150)
          if (startsWith(left, "coef_")) {{
            plotMD(fit2_base, coef = expr, main = paste0("MD: ", expr))
          }} else {{
            plotMD(fit2_con, coef = expr, main = paste0("MD: ", expr))
          }}
          dev.off()
        }}

        # Optional PCA.
        if (requireNamespace("ggplot2", quietly = TRUE)) {{
          suppressPackageStartupMessages(library(ggplot2))
          pc <- prcomp(t(mat), scale. = TRUE)
          df <- data.frame(
            sample = rownames(pc$x),
            PC1 = pc$x[, 1],
            PC2 = pc$x[, 2]
          )
          if (!is.null(group_col) && group_col %in% colnames(pheno)) {{
            df$group <- as.factor(pheno[[group_col]])
          }}
          p <- ggplot(df, aes(PC1, PC2)) +
            geom_point(size = 3, alpha = 0.9) +
            theme_minimal()
          if ("group" %in% colnames(df)) {{
            p <- p + aes(color = group) + labs(color = group_col)
          }}
          ggsave(filename = file.path(out_dir, "pca_plot.png"), plot = p, width = 6, height = 5, dpi = 150)
        }}
        ```

        ## Provenance

        ```{{r context-dump, results='asis'}}
        cat("<details><summary><strong>Resolved context (JSON)</strong></summary><pre>")
        cat(jsonlite::toJSON(ctx, pretty = TRUE, auto_unbox = TRUE))
        cat("</pre></details>")
        ```

        ```{{r session-info}}
        sessionInfo()
        ```
        """
    ).strip() + "\n"


def _render_render_sh(*, rmd_name: str = "report.Rmd") -> str:
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env bash
        set -euo pipefail

        Rscript -e 'rmarkdown::render("{rmd_name}")'
        """
    )


def write_limma_replay_bundle(
    *,
    report_dir: str,
    replay_dir: str,
    title: Optional[str] = None,
    copy_inputs: bool = True,
    force: bool = False,
    render: bool = False,
) -> RmdReplayBundle:
    report_path = Path(report_dir).expanduser().resolve()
    report_path.mkdir(parents=True, exist_ok=True)
    if not report_path.is_dir():
        raise NotADirectoryError(str(report_path))

    replay_path = Path(replay_dir).expanduser().resolve()
    gct_src, ctx_src = validate_replay_dir(replay_path)

    now = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    title = title or f"Tackle limma replay: {replay_path.name}"

    gct_dst_rel = Path("data") / "limma_input.gct"
    ctx_dst_rel = Path("context") / "limma_replay_context.json"
    gct_dst = report_path / gct_dst_rel
    ctx_dst = report_path / ctx_dst_rel
    _copy_or_symlink(gct_src, gct_dst, copy=copy_inputs, force=force)
    _copy_or_symlink(ctx_src, ctx_dst, copy=copy_inputs, force=force)

    rmd_text = _render_limma_replay_rmd(
        title=title,
        generated_at=now,
        gct_relpath=str(gct_dst_rel).replace("\\", "/"),
        context_relpath=str(ctx_dst_rel).replace("\\", "/"),
    )
    rmd_path = report_path / "report.Rmd"
    _safe_write_text(rmd_path, rmd_text, force=force)

    render_sh_path = report_path / "render.sh"
    _safe_write_text(render_sh_path, _render_render_sh(rmd_name=rmd_path.name), force=force)
    try:
        render_sh_path.chmod(render_sh_path.stat().st_mode | 0o111)
    except OSError:
        pass

    html_path: Optional[Path] = None
    if render:
        cmd = ["Rscript", "-e", f'rmarkdown::render("{rmd_path.name}")']
        try:
            subprocess.run(cmd, cwd=str(report_path), check=True)
        except FileNotFoundError as e:
            raise RuntimeError("Rscript not found on PATH") from e
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to render Rmd (exit={e.returncode})") from e
        candidate = report_path / "report.html"
        if candidate.exists():
            html_path = candidate

    return RmdReplayBundle(
        report_dir=report_path,
        rmd_path=rmd_path,
        html_path=html_path,
        gct_path=gct_dst,
        context_path=ctx_dst,
        render_sh=render_sh_path,
    )

