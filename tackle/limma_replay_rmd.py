from __future__ import annotations

import textwrap


def render_limma_replay_rmd(
    *,
    title: str,
    generated_at: str,
    gct_relpath: str,
    context_relpath: str,
    outputs_rel_dir: str,
) -> str:
    """Render the single authoritative R Markdown implementation of limma replay."""
    title_yaml = str(title).replace("\\", "\\\\").replace('"', '\\"')
    generated_yaml = str(generated_at).replace("\\", "\\\\").replace('"', '\\"')

    return textwrap.dedent(
        f"""\
        ---
        title: "{title_yaml}"
        date: "{generated_yaml}"
        output:
          html_document:
            toc: true
            toc_float: true
            df_print: paged
        ---

        ```{{r setup, include=FALSE}}
        knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
        ```

        ## Replay Contract

        - Analysis replayed here: limma model, contrasts, empirical Bayes moderation, and result tables.
        - Authoritative input: `{gct_relpath}`, the matrix supplied to limma by tackle.
        - Imputation is **not** rerun. If the original analysis imputed values, they are already in the stored matrix.
        - PCA is intentionally excluded. Tackle `pca2` has a separate input matrix and preprocessing contract.

        ```{{r load-context}}
        if (!requireNamespace("jsonlite", quietly = TRUE)) {{
          stop("Missing dependency: jsonlite")
        }}
        ctx <- jsonlite::read_json("{context_relpath}", simplifyVector = TRUE)

        key_fields <- c(
          "replay_contract_version", "replay_scope", "analysis_dir", "volcano_dir", "run_id",
          "matrix_shape", "stored_matrix_role", "stored_matrix_is_authoritative",
          "expression_matrix_source",
          "group", "formula_in", "formula_effective", "formula_rewritten",
          "contrasts_spec", "block", "limma_robust", "limma_trend",
          "impute_missing_values", "imputation_backend", "gaussian_method", "lupine_mode",
          "has_pre_impute_matrix", "pre_impute_matrix_role", "recompute_imputation_supported",
          "normtype", "non_zeros", "taxon", "batch_applied", "batch_method",
          "pca_included"
        )
        present <- intersect(key_fields, names(ctx))
        as.data.frame(t(unlist(ctx[present])), stringsAsFactors = FALSE)
        ```

        ## Load the Modeled Matrix

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

        stored_ds <- parse_gct_any("{gct_relpath}")
        mat_stored <- stored_ds@mat
        storage.mode(mat_stored) <- "double"
        cdesc <- as.data.frame(stored_ds@cdesc, stringsAsFactors = FALSE)
        rdesc <- as.data.frame(stored_ds@rdesc, stringsAsFactors = FALSE)

        if ("id" %in% colnames(cdesc)) rownames(cdesc) <- as.character(cdesc$id)
        missing_samples <- setdiff(colnames(mat_stored), rownames(cdesc))
        if (length(missing_samples) > 0) {{
          stop(paste("Stored metadata is missing samples:", paste(missing_samples, collapse = ", ")))
        }}
        cdesc <- cdesc[colnames(mat_stored), , drop = FALSE]
        stopifnot(identical(rownames(cdesc), colnames(mat_stored)))

        if ("GeneID" %in% colnames(rdesc)) {{
          rownames(rdesc) <- as.character(rdesc$GeneID)
        }} else if ("id" %in% colnames(rdesc)) {{
          rownames(rdesc) <- as.character(rdesc$id)
        }}

        if (!isTRUE(ctx$stored_matrix_is_authoritative)) {{
          stop("Replay context does not mark the stored matrix as authoritative.")
        }}
        if (!is.null(ctx$matrix_shape)) {{
          expected_shape <- as.integer(unlist(ctx$matrix_shape, use.names = FALSE))
          if (!identical(as.integer(dim(mat_stored)), expected_shape)) {{
            stop(
              "Stored matrix shape does not match replay context: ",
              paste(dim(mat_stored), collapse = " x "), " versus ",
              paste(expected_shape, collapse = " x ")
            )
          }}
        }}

        mat <- mat_stored
        cat("Modeled matrix dimensions:", paste(dim(mat), collapse = " x "), "\\n")
        cat("Missing modeled values:", sum(is.na(mat)), "\\n")
        cat("No imputation is performed by this document.\\n")
        ```

        ## Design and Contrasts

        ```{{r design-contrasts}}
        if (!requireNamespace("limma", quietly = TRUE)) {{
          stop("Missing dependency: limma")
        }}
        suppressPackageStartupMessages(library(limma))

        pheno <- cdesc
        group_col <- if (!is.null(ctx$group) && length(ctx$group) == 1 && nzchar(ctx$group)) {{
          as.character(ctx$group)
        }} else {{
          NULL
        }}

        formula_str <- NULL
        if (!is.null(ctx$formula_rewritten) && length(ctx$formula_rewritten) == 1 && nzchar(ctx$formula_rewritten)) {{
          formula_str <- as.character(ctx$formula_rewritten)
        }} else if (!is.null(ctx$formula_effective) && length(ctx$formula_effective) == 1 && nzchar(ctx$formula_effective)) {{
          formula_str <- as.character(ctx$formula_effective)
        }}

        if (is.null(formula_str) && is.null(group_col)) {{
          stop("Context is missing both formula and group.")
        }}
        if (!is.null(group_col) && !(group_col %in% colnames(pheno))) {{
          stop(paste0("Missing group column in cdesc: ", group_col))
        }}
        if (!is.null(ctx$block) && length(ctx$block) == 1 && nzchar(ctx$block) && !(ctx$block %in% colnames(pheno))) {{
          stop(paste0("Missing block column in cdesc: ", ctx$block))
        }}
        if (!is.null(formula_str) && grepl("ns\\\\(", formula_str)) {{
          if (!requireNamespace("splines", quietly = TRUE)) stop("Missing dependency: splines")
          suppressPackageStartupMessages(library(splines))
        }}

        sanitize_name <- function(x) {{
          x <- gsub(":", "_", x, fixed = TRUE)
          x <- gsub(" ", "_", x, fixed = TRUE)
          x <- gsub("-", "_", x, fixed = TRUE)
          x <- gsub("\\\\+", "_", x)
          x <- gsub("/", "_", x, fixed = TRUE)
          x <- gsub("\\\\?", "qmk", x)
          x
        }}

        if (is.null(formula_str)) {{
          mod <- model.matrix(as.formula(paste0("~0+", group_col)), pheno)
        }} else {{
          mod <- model.matrix(as.formula(formula_str), pheno)
        }}
        colnames(mod) <- make.names(sanitize_name(colnames(mod)))

        as_character_vector <- function(x) {{
          if (is.null(x) || length(x) == 0) return(character())
          as.character(unlist(x, use.names = FALSE))
        }}
        direct_coefs <- as_character_vector(ctx$direct_coefs)
        contrast_exprs <- as_character_vector(ctx$contrast_expressions)
        result_labels <- as_character_vector(ctx$result_labels)
        if (length(result_labels) == 0) {{
          result_labels <- c(paste0("coef_", direct_coefs, "=", direct_coefs),
                             paste0(contrast_exprs, "=", contrast_exprs))
        }}
        if (length(result_labels) == 0) stop("Replay context contains no limma results to reproduce.")

        contrasts_matrix <- NULL
        if (length(contrast_exprs) > 0) {{
          contrasts_matrix <- makeContrasts(contrasts = contrast_exprs, levels = mod)
        }}

        mod
        if (!is.null(contrasts_matrix)) contrasts_matrix
        ```

        ## Limma Fit

        ```{{r limma-fit}}
        robust_flag <- isTRUE(ctx$limma_robust)
        trend_flag <- isTRUE(ctx$limma_trend)

        block <- NULL
        cor_value <- NULL
        if (!is.null(ctx$block) && length(ctx$block) == 1 && nzchar(ctx$block)) {{
          block <- as.factor(pheno[[ctx$block]])
          corfit <- duplicateCorrelation(mat_stored, design = mod, block = block)
          cor_value <- corfit$consensus
          cat("duplicateCorrelation consensus:", cor_value, "\\n")
        }}

        fit <- lmFit(
          as.matrix(mat_stored),
          mod,
          block = block,
          correlation = cor_value
        )
        fit2_base <- eBayes(fit, robust = robust_flag, trend = trend_flag)

        fit2_con <- NULL
        if (!is.null(contrasts_matrix)) {{
          fit2_con <- contrasts.fit(fit, contrasts_matrix)
          fit2_con <- eBayes(fit2_con, robust = robust_flag, trend = trend_flag)
        }}

        cat("limma fit complete; robust =", robust_flag, "; trend =", trend_flag, "\\n")
        ```

        ## Replayed Result Tables

        ```{{r replay-results}}
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

        split_result_label <- function(label) {{
          hit <- regexpr("=", label, fixed = TRUE)[1]
          if (hit < 0) return(list(left = label, expr = label))
          list(
            left = trimws(substr(label, 1, hit - 1)),
            expr = trimws(substr(label, hit + 1, nchar(label)))
          )
        }}

        normalize_tackle_columns <- function(tt) {{
          rename_one <- function(old, new) {{
            hit <- match(old, colnames(tt), nomatch = 0)
            if (hit > 0) colnames(tt)[hit] <<- new
          }}
          rename_one("logFC", "log2_FC")
          rename_one("adj.P.Val", "pAdj")
          rename_one("P.Value", "pValue")
          for (column in intersect(c("log2_FC", "CI.L", "CI.R"), colnames(tt))) {{
            tt[[column]] <- tt[[column]] / log10(2)
          }}
          tt
        }}

        add_replay_data <- function(tt) {{
          gene_ids <- rownames(tt)
          tt$GeneID <- gene_ids

          meta_cols <- setdiff(colnames(rdesc), c("GeneID", "id"))
          if (length(meta_cols) > 0) {{
            gene_meta <- rdesc[gene_ids, meta_cols, drop = FALSE]
            tt <- cbind(tt, gene_meta)
          }}

          expression_values <- as.data.frame(
            mat_stored[gene_ids, , drop = FALSE],
            check.names = FALSE
          )
          cbind(tt, expression_values)
        }}

        tables <- list()
        manifest <- data.frame(
          label = character(), kind = character(), coefficient = character(), path = character(),
          stringsAsFactors = FALSE
        )

        for (label in result_labels) {{
          spec <- split_result_label(label)
          if (startsWith(spec$left, "coef_")) {{
            tt <- topTable(fit2_base, n = Inf, sort.by = "none", coef = spec$expr, confint = TRUE)
            kind <- "coefficient"
          }} else {{
            if (is.null(fit2_con)) stop("Contrast result requested but the contrast fit is NULL: ", label)
            tt <- topTable(fit2_con, n = Inf, sort.by = "none", coef = spec$expr, confint = TRUE)
            kind <- "contrast"
          }}

          tt <- normalize_tackle_columns(tt)
          tt <- add_replay_data(tt)
          out_path <- file.path(out_dir, paste0("limma_", safe_file(spec$left), ".tsv"))
          write.table(tt, out_path, sep = "\\t", quote = FALSE, row.names = FALSE)
          tables[[label]] <- tt
          manifest <- rbind(
            manifest,
            data.frame(
              label = label,
              kind = kind,
              coefficient = spec$expr,
              path = out_path,
              stringsAsFactors = FALSE
            )
          )
        }}

        write.table(
          manifest,
          file.path(out_dir, "limma_replay_manifest.tsv"),
          sep = "\\t", quote = FALSE, row.names = FALSE
        )
        manifest
        if (length(tables) > 0) head(tables[[1]], 20)
        ```

        ## Limma Diagnostic

        ```{{r limma-diagnostic}}
        if (length(result_labels) > 0) {{
          spec <- split_result_label(result_labels[[1]])
          md_path <- file.path(out_dir, paste0("md_plot_", safe_file(spec$left), ".png"))
          png(md_path, width = 1200, height = 900, res = 150)
          if (startsWith(spec$left, "coef_")) {{
            plotMD(fit2_base, coef = spec$expr, main = paste0("MD: ", spec$expr))
          }} else {{
            plotMD(fit2_con, coef = spec$expr, main = paste0("MD: ", spec$expr))
          }}
          dev.off()
          md_path
        }}
        ```

        ## Modeled Matrix QC

        These summaries inspect the exact modeled matrix without changing it.

        ```{{r matrix-summary}}
        sample_summary <- data.frame(
          sample = colnames(mat_stored),
          mean = colMeans(mat_stored, na.rm = TRUE),
          median = apply(mat_stored, 2, median, na.rm = TRUE),
          sd = apply(mat_stored, 2, sd, na.rm = TRUE),
          na_count = colSums(is.na(mat_stored)),
          stringsAsFactors = FALSE
        )
        knitr::kable(sample_summary)
        ```

        ## PCA Scope

        PCA is not computed here. Replaying tackle `pca2` requires its own saved matrix plus the
        original `fillna`, `center`, `scale`, feature-selection, and annotation options; the limma
        modeled matrix and context do not establish that contract.

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


def render_limma_replay_sh(*, rmd_name: str, html_name: str | None = None) -> str:
    output_arg = "" if html_name is None else f", output_file = '{html_name}'"
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env bash
        set -euo pipefail
        cd "$(dirname "$0")"
        Rscript -e "rmarkdown::render('{rmd_name}'{output_arg})"
        """
    )
