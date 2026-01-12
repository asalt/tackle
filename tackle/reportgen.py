from __future__ import annotations

import shutil
import textwrap
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from .scriptgen import ConfigSummary, summarize_config


@dataclass(frozen=True)
class ReportBundle:
    report_dir: str
    rmd_path: str
    gct_path: str
    render_sh: str


def _safe_write_text(path: Path, content: str, *, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"{path} already exists (use --force to overwrite)")
    path.write_text(content)


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


def _render_report_rmd(
    *,
    title: str,
    gct_relpath: str,
    generated_at: str,
    tackle_command: Optional[str],
    config_summary: Optional[ConfigSummary],
) -> str:
    meta_lines = []
    if config_summary is not None:
        meta_lines.append(f"- Config: `{config_summary.analysis_name}`")
        meta_lines.append(f"- Samples: `{config_summary.total_samples}`")
        if config_summary.metadata_columns:
            meta_lines.append(
                "- Metadata columns: " + ", ".join(f"`{c}`" for c in config_summary.metadata_columns)
            )
        if config_summary.recommended_design_columns:
            meta_lines.append(
                "- Suggested design cols: "
                + ", ".join(f"`{c}`" for c in config_summary.recommended_design_columns)
            )
    if tackle_command:
        meta_lines.append(f"- Tackle command: `{tackle_command}`")
    meta_block = "\n".join(meta_lines) if meta_lines else ""

    # Keep dependencies minimal; cmapR is the only hard requirement to parse GCT.
    # Limma/ggplot2 are used in optional example sections.
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

        ## Provenance
        {meta_block}

        ## Load the GCT

        ```{{r setup, include=FALSE}}
        knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
        ```

        ```{{r load-gct}}
        library(cmapR)

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

        dim(mat)
        head(cdesc)
        head(rdesc)
        ```

        ## Quick QC

        ```{{r qc}}
        summary(as.numeric(mat))
        frac_na_per_sample <- colMeans(is.na(mat))
        summary(frac_na_per_sample)
        ```

        ## PCA (example)

        ```{{r pca}}
        if (!requireNamespace("ggplot2", quietly = TRUE)) {{
          stop("Install ggplot2 to run PCA plotting.")
        }}
        library(ggplot2)

        # TODO: choose a grouping column from cdesc
        group_col <- "{(config_summary.recommended_design_columns[0] if (config_summary and config_summary.recommended_design_columns) else 'group')}"
        if (!group_col %in% colnames(cdesc)) {{
          stop(paste("Set group_col to one of:", paste(colnames(cdesc), collapse = ", ")))
        }}
        grp <- as.factor(cdesc[[group_col]])

        # Filter rows with too many missing values
        keep <- rowSums(is.na(mat)) < ncol(mat) * 0.5
        m <- mat[keep, , drop = FALSE]

        # Simple imputation for PCA only: fill NA with per-row median
        row_med <- apply(m, 1, median, na.rm = TRUE)
        for (i in seq_len(nrow(m))) {{
          na_idx <- is.na(m[i, ])
          if (any(na_idx)) m[i, na_idx] <- row_med[i]
        }}

        pc <- prcomp(t(m), scale. = TRUE)
        pca_df <- data.frame(
          sample = rownames(pc$x),
          grp = grp,
          PC1 = pc$x[, 1],
          PC2 = pc$x[, 2]
        )

        ggplot(pca_df, aes(PC1, PC2, color = grp)) +
          geom_point(size = 3) +
          theme_minimal() +
          labs(color = group_col)
        ```

        ## Limma (example)

        ```{{r limma}}
        if (!requireNamespace("limma", quietly = TRUE)) {{
          stop("Install limma to run differential expression.")
        }}
        library(limma)

        # Reuse `m` from PCA section (filtered + imputed). For a stricter analysis,
        # consider handling missingness differently.
        design <- model.matrix(~ 0 + grp)
        colnames(design) <- levels(grp)

        fit <- lmFit(m, design)
        if (nlevels(grp) < 2) {{
          stop("Need at least 2 groups in grp for contrasts.")
        }}
        contrast_str <- paste(levels(grp)[1], "-", levels(grp)[2])
        contrast <- makeContrasts(contrasts = contrast_str, levels = design)
        fit2 <- contrasts.fit(fit, contrast)
        fit2 <- eBayes(fit2)

        top <- topTable(fit2, number = 20)
        head(top)
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


def write_report_bundle(
    *,
    report_dir: str,
    gct_path: str,
    title: Optional[str] = None,
    conf_path: Optional[str] = None,
    tackle_command: Optional[str] = None,
    copy_gct: bool = True,
    force: bool = False,
) -> ReportBundle:
    report_path = Path(report_dir)
    report_path.mkdir(parents=True, exist_ok=True)
    if not report_path.is_dir():
        raise NotADirectoryError(str(report_path))

    src_gct = Path(gct_path)
    if not src_gct.exists():
        raise FileNotFoundError(str(src_gct))

    config_summary = summarize_config(conf_path) if conf_path else None

    normalized_name = "normalized" + src_gct.suffix
    gct_dst_rel = Path("data") / normalized_name
    gct_dst = report_path / gct_dst_rel
    _copy_or_symlink(src_gct, gct_dst, copy=copy_gct, force=force)

    if conf_path:
        _copy_or_symlink(Path(conf_path), report_path / "config.conf", copy=True, force=force)

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    title = title or f"Tackle report: {src_gct.stem}"
    rmd = _render_report_rmd(
        title=title,
        gct_relpath=str(gct_dst_rel).replace("\\", "/"),
        generated_at=now,
        tackle_command=tackle_command,
        config_summary=config_summary,
    )

    rmd_path = report_path / "report.Rmd"
    _safe_write_text(rmd_path, rmd, force=force)

    render_sh_path = report_path / "render.sh"
    _safe_write_text(render_sh_path, _render_render_sh(rmd_name=rmd_path.name), force=force)
    try:
        mode = render_sh_path.stat().st_mode
        render_sh_path.chmod(mode | 0o111)
    except OSError:
        pass

    if tackle_command:
        cmd_path = report_path / "tackle_command.sh"
        cmd_content = textwrap.dedent(
            f"""\
            #!/usr/bin/env bash
            set -euo pipefail

            # Generated by: tackle make-report
            # This captures the command used when the report bundle was created.

            {tackle_command}
            """
        )
        _safe_write_text(cmd_path, cmd_content, force=force)
        try:
            mode = cmd_path.stat().st_mode
            cmd_path.chmod(mode | 0o111)
        except OSError:
            pass

    return ReportBundle(
        report_dir=str(report_path),
        rmd_path=str(rmd_path),
        gct_path=str(gct_dst),
        render_sh=str(render_sh_path),
    )
