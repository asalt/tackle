from __future__ import annotations

import json
import textwrap

from .template_rendering import render_text_template


def render_pca_replay_rmd(
    *,
    title: str,
    gct_relpath: str = "pca_input_pre_svd.gctx",
    context_relpath: str = "pca_replay_context.json",
    include_separation: bool = False,
) -> str:
    """Render a standalone pca2 replay with tackle's score/biplot styling."""
    return render_text_template(
        "pca_replay.Rmd.j2",
        title_json=json.dumps(str(title)),
        gct_relpath=str(gct_relpath).replace("\\", "/"),
        context_relpath=str(context_relpath).replace("\\", "/"),
        include_separation=bool(include_separation),
    )


def render_pca_replay_sh(*, rmd_name: str = "replot.Rmd") -> str:
    return textwrap.dedent(
        f'''\
        #!/usr/bin/env bash
        set -euo pipefail
        cd "$(dirname "$0")"
        Rscript -e "rmarkdown::render('{rmd_name}')"
        '''
    )
