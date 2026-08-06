# Tackle templates

This directory contains editable source templates for generated scripts and reports.

- `limma_replay.Rmd.j2` is the authoritative limma replay document.
- `pca_replay.Rmd.j2` is the authoritative `pca2` replay document.
- The remaining `.j2` files generate shell scripts or HTML/Markdown reports.

Replay templates are rendered by `tackle/template_rendering.py`. They use
`[[[ value ]]]` for variables and `[% ... %]` for blocks so ordinary R braces
remain untouched.
