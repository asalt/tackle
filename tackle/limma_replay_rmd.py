from __future__ import annotations

import json
import textwrap
from collections.abc import Mapping, Sequence
from typing import Any

from .template_rendering import render_text_template


def _r_string(value: Any) -> str:
    if value is None or str(value) == "":
        return "NULL"
    return json.dumps(str(value), ensure_ascii=False)


def _r_string_vector(values: Any) -> str:
    if values is None:
        return "character()"
    if isinstance(values, (str, bytes)):
        items = [values]
    elif isinstance(values, Sequence):
        items = list(values)
    else:
        items = [values]
    if not items:
        return "character()"
    rendered = ",\n  ".join(_r_string(value) for value in items)
    return f"c(\n  {rendered}\n)"


def _r_matrix_shape(value: Any) -> str:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return "NULL"
    dimensions = [int(dimension) for dimension in value]
    if not dimensions:
        return "NULL"
    return "c(" + ", ".join(f"{dimension}L" for dimension in dimensions) + ")"


def render_limma_replay_rmd(
    *,
    title: str,
    generated_at: str,
    gct_relpath: str,
    outputs_rel_dir: str,
    replay_context: Mapping[str, Any],
) -> str:
    """Render the single authoritative R Markdown implementation of limma replay."""
    title_yaml = str(title).replace("\\", "\\\\").replace('"', '\\"')
    generated_yaml = str(generated_at).replace("\\", "\\\\").replace('"', '\\"')
    formula_str = (
        replay_context.get("formula_rewritten")
        or replay_context.get("formula_effective")
        or replay_context.get("formula_in")
    )
    return (
        render_text_template(
            "limma_replay.Rmd.j2",
            title_yaml=title_yaml,
            generated_yaml=generated_yaml,
            gct_relpath=gct_relpath,
            outputs_rel_dir=outputs_rel_dir,
            expected_matrix_shape=_r_matrix_shape(replay_context.get("matrix_shape")),
            group_col=_r_string(replay_context.get("group")),
            formula_str=_r_string(formula_str),
            block_col=_r_string(replay_context.get("block")),
            robust_flag="TRUE" if bool(replay_context.get("limma_robust")) else "FALSE",
            trend_flag="TRUE" if bool(replay_context.get("limma_trend")) else "FALSE",
            direct_coefs=_r_string_vector(replay_context.get("direct_coefs")),
            contrast_exprs=_r_string_vector(
                replay_context.get("contrast_expressions")
            ),
            result_labels=_r_string_vector(replay_context.get("result_labels")),
        ).strip()
        + "\n"
    )


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
