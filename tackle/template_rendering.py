from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, StrictUndefined


@lru_cache(maxsize=1)
def _template_environment() -> Environment:
    """Return the shared environment for non-HTML text templates."""
    return Environment(
        loader=FileSystemLoader(str(Path(__file__).resolve().parent / "templates")),
        autoescape=False,
        undefined=StrictUndefined,
        keep_trailing_newline=True,
        variable_start_string="[[[",
        variable_end_string="]]]",
        block_start_string="[%",
        block_end_string="%]",
        comment_start_string="[#",
        comment_end_string="#]",
    )


def render_text_template(template_name: str, **context: object) -> str:
    """Render a package-owned text template with strict context checking."""
    return _template_environment().get_template(template_name).render(**context)
