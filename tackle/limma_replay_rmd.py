from __future__ import annotations

import textwrap

from .template_rendering import render_text_template


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
    return (
        render_text_template(
            "limma_replay.Rmd.j2",
            title_yaml=title_yaml,
            generated_yaml=generated_yaml,
            gct_relpath=gct_relpath,
            context_relpath=context_relpath,
            outputs_rel_dir=outputs_rel_dir,
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
