---
name: tackle-cli
description: Construct and debug the tackle CLI (proteomics analysis pipeline) in this repo, including the chained invocation pattern `tackle [GLOBAL_OPTS] EXPERIMENTS.conf subcmd1 ... subcmdN`, the standalone helper subcommands that run without a .conf (make-run/make-cluster/make-xls/etc.), and quick mapping from subcommands to the relevant Python/R implementation files.
---

# Tackle CLI

Keep this skill lightweight: focus on the non-obvious CLI structure (chaining + config-file position) and where to look in the codebase.

## Overview

Tackle is a chained Click CLI for proteomics analysis. Most runs load an experiment config once (`EXPERIMENTS.conf`) and then execute one or more subcommands that write artifacts into a structured results directory.

## CLI shape (the “config file terminator”)

- Treat `tackle` as a Click `@group(chain=True)` where `EXPERIMENTS.conf` is a required positional argument for analysis runs.
- Put **global options before** `EXPERIMENTS.conf`. After that, chain one or more subcommands.
- In a chained run, the `Data` object is loaded once (shared across subcommands).

Template:

```bash
tackle [GLOBAL_OPTS] EXPERIMENTS.conf <subcmd1> [subcmd1 opts] <subcmd2> [subcmd2 opts] ...
```

Example (typical workflow):

```bash
tackle --result-dir results --name run1 EXPERIMENTS.conf \
  export --level align \
  volcano --contrast A,B \
  pca2 --color plex --marker rep \
  umap --color plex
```

## Standalone subcommands and “help without a config”

The first positional argument is implemented as a “path-or-subcommand”:

- `tackle volcano` prints volcano help (it does not run; it exits after printing help).
- Some utilities run **standalone** by occupying the config position:
  - `make_config`, `make-run`, `make-cluster`, `make-xls`, `make-deck` / `make-pptx`, `make-html`, `make-rmd`,
    `make-search-report`, `cluster-summary`, `replot_gsea`.

Examples:

```bash
tackle make-run --conf EXPERIMENTS.conf --out tacklerun_EXPERIMENTS.sh
tackle make-cluster --conf EXPERIMENTS.conf --k-start 3 --k-end 20
```

## Output layout (high-level)

- Base output dir: `--result-dir` (default `./results`)
- Analysis root: `RESULT_DIR/CONF_STEM/NAME/` (NAME is optional `--name`)
- Common subfolders created by subcommands:
  - `export/`, `volcano/`, `pca/`, `umap/`
  - `report/deck/` (PPTX/table assets), `report/html/` (static overview), `report/rmd/` (limma replay bundle)

## Subcommands: quick mental map

- `export`: write matrices/tables.
- `volcano`: differential + volcano plots (limma). Notable toggles: `--limma-robust/--no-limma-robust`, `--limma-trend/--no-limma-trend`.
- `pca2`: PCA via R (ggplot/ggfortify). NA handling: `--fillna {min,avg}`; scaling: `--center/--scale`; optional vectors: `--show-loadings`.
- `umap`: UMAP via R (`uwot`). Writes embedding/params TSVs and plots under `umap/`.
- `make-html`: build a static HTML overview bundle (`index.html` + `assets/`) for quick browsing of key PNGs + lightweight volcano TSV previews.
- `make-deck` / `make-pptx`: build slide/table assets and optionally a `.pptx` deck (requires `python-pptx`).
- `make-rmd`: build a standalone limma replay Rmd bundle from an existing volcano run (copies `limma_input.gct` + `limma_replay_context.json` into `report/rmd/` and writes `report.Rmd` + `render.sh`).

Use `tackle --help` for the full list, and `tackle SUBCMD` to print a specific subcommand’s help quickly.

## Code map (where to edit)

- CLI/options/wiring: `tackle/main.py`
- Core data plumbing + ComBat + limma integration: `tackle/containers.py`
- Limma runner: `tackle/statmodels/limma_runner.py`
- Volcano: `tackle/volcanoplot.py`, `tackle/R/volcanoplot.R`
- PCA2: `tackle/R/pcaplot.R`
- UMAP: `tackle/R/umapplot.R`
- HTML overview report: `tackle/html_overview.py`, `tackle/templates/overview_report.html.j2`
- Slide deck generator: `tackle/deckgen.py` (used by `make-deck` / `make-pptx`)
- Limma replay bundle: `tackle/limma_replay.py`, `tackle/rmd_replay.py`
