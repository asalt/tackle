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

High-level sample similarity examples:

```bash
tackle [GLOBAL_OPTS] EXPERIMENTS.conf \
  correlation --metric l2 --cluster --linkage auto \
    --cut-by group --annotate \
  correlation --metric pearson --z-score --no-cluster \
    --cut-by fraction:treatment --legend-exclude replicate \
    --sample-exclude SAMPLE_TO_DROP --file-format .png
```

## Standalone subcommands and “help without a config”

The first positional argument is implemented as a “path-or-subcommand”:

- `tackle volcano` prints volcano help (it does not run; it exits after printing help).
- Some utilities run **standalone** by occupying the config position:
  - `make_config`, `make-run`, `make-cluster`, `make-xls`, `make-deck` / `make-pptx`, `make-html`, `make-rmd`,
    `make-search-report`, `cluster-summary`, `replot_gsea`, and `dev inspect-gctx`.

Examples:

```bash
tackle make-run --conf EXPERIMENTS.conf --out tacklerun_EXPERIMENTS.sh
tackle make-cluster --conf EXPERIMENTS.conf --k-start 3 --k-end 20
tackle dev inspect-gctx results/run/correlation/input.<hash>.gctx
tackle dev inspect-gctx --spec
```

Generated `make-run` scripts include editable `run_correlation` variants and a
`run_pca` default using the selected design column, PCs 1-3, sample labels,
encircles, and heteroscedastic separation tests.

## Output layout (high-level)

- Base output dir: `--result-dir` (default `./results`)
- Analysis root: `RESULT_DIR/CONF_STEM/NAME/` (NAME is optional `--name`)
- Common subfolders created by subcommands:
  - `export/`, `volcano/`, `pca/`, `correlation/`, `umap/`
  - `report/deck/` (PPTX/table assets), `report/html/` (static overview), `report/rmd/` (limma replay bundle)

## Subcommands: quick mental map

- `export`: write matrices/tables. `--level zscore` exports the shared detection-aware global feature z-score of the masked logged matrix.
- `volcano`: differential + volcano plots (limma). Notable toggles: `--limma-robust/--no-limma-robust`, `--limma-trend/--no-limma-trend`.
- `correlation`: sample-level L2/L1 distance or Pearson/Spearman correlation heatmaps from masked logged values. `--metric l2` is the default and shows pairwise-complete RMS distance; `l1` shows pairwise-complete mean absolute distance. The accepted aliases are `euclidean -> l2` and `manhattan -> l1`. Pearson/Spearman cells show raw r while clustering privately uses `sqrt(2 * (1 - r))` chord distance. `--linkage auto` resolves to `average` for L1 and `ward.D2` otherwise; explicit linkage values remain available. `--z-score` mirrors cluster2's detection-aware `myzscore`: temporarily fill nondetections at the observed feature minimum minus its SD, scale, then restore the NA mask before the pairwise-observed metric. `--sample-exclude`, `--legend-include`, and `--legend-exclude` mirror cluster2-style filtering. Also supports repeatable/colon-delimited `--cut-by` and opt-in `--annotate`; exports content-addressed metadata-rich GCTX inputs/summary matrices, one shared pairwise-count TSV, and settings.
- `pca2`: PCA via R (ggplot/ggfortify). NA handling: `--fillna {min,avg}`; scaling: `--center/--scale`; optional vectors: `--show-loadings`. `--test-by FIELD` opts into plane-wise Euclidean R2 plus heteroscedastic Johansen Welch-James tests, adjusted global/pairwise tables, and optional plot captions; captions name the actual adjustment (for example, `Holm-adjusted p`). The omnibus TSV is strictly factor-level, while the pairwise TSV contains one plotting-ready row for every group pair even when the factor has exactly two levels. Displayed planes also get paired geometry/ranking summaries of centroid distance, pooled RMS spread, standardized separation, R2, and adjusted significance. Its default leading-space test selects the largest PC1 through PCk block for which the omnibus and every required pairwise test are jointly estimable, merges an identical displayed plane instead of duplicating the hypothesis, and reports the exact tested components' explained-variance percentage. Without `--test-by`, inference stays off. The returned scree plot uses a preprocessing-only shared stem, so chained aesthetic/test variants reuse one scree artifact. Every run also writes `pca/pca2/replay/<run-id>/`, whose GCTX is the exact sample-by-feature matrix immediately before `prcomp`; the standalone `replot.Rmd` reuses it without filling, centering, scaling, normalizing, or transposing and recomputes the requested tests and pairwise summaries.
- `umap`: UMAP via R (`uwot`). Writes embedding/params TSVs and plots under `umap/`.
- `make-html`: build a static HTML overview bundle (`index.html` + `assets/`) for quick browsing of key PNGs and lightweight volcano TSV previews. Correlation plots discovered under `correlation/` appear in their own tab immediately after Metrics. In-page plots default to a 1200 px cap (`--plot-max-width-px 0` restores full width), while click-through retains the original resolution. `--pngquant` uses bounded concurrent subprocesses by default; control them with `--pngquant-workers` (`0` means auto, up to eight).
- `dev inspect-gctx`: inspect GCTX dimensions, metadata fields, HDF5 storage, and tackle provenance without loading matrix values. `--json` emits machine-readable output and `--spec` includes the writer/hash contract. Auto hashing prefers BLAKE3 and falls back to `hashlib` BLAKE2b-256; the actual algorithm is stored in each tackle-authored GCTX.
- `make-deck` / `make-pptx`: build slide/table assets and optionally a `.pptx` deck (requires `python-pptx`).
- `make-rmd`: build a standalone limma-only replay from the authoritative modeled matrix (no re-imputation or PCA); it copies `limma_input.gct` + `limma_replay_context.json` into `report/rmd/` and writes `report.Rmd` + `render.sh`.

Use `tackle --help` for the full list, and `tackle SUBCMD` to print a specific subcommand’s help quickly.

## Code map (where to edit)

- CLI/options/wiring: `tackle/main.py`
- Core data plumbing + ComBat + limma integration: `tackle/containers.py`
- Limma runner: `tackle/statmodels/limma_runner.py`
- Volcano: `tackle/volcanoplot.py`, `tackle/R/volcanoplot.R`
- Correlation heatmap: `tackle/correlationplot.py`, `tackle/R/correlationplot.R`
- Shared Python detection-aware z-score: `tackle/zscore.py` (also re-exported by `tackle.containers` for `export --level zscore`); cluster2's matching R implementation is in `tackle/R/clusterplot.R`
- GCT/GCTX writer, hashing contract, and inspection: `tackle/gct_io.py`
- PCA2 and exact pre-SVD replay: `tackle/R/pcaplot.R`, `tackle/pca_replay.py`, `tackle/pca_replay_rmd.py`
- UMAP: `tackle/R/umapplot.R`
- HTML overview report: `tackle/html_overview.py`, `tackle/templates/overview_report.html.j2`
- Slide deck generator: `tackle/deckgen.py` (used by `make-deck` / `make-pptx`)
- Limma replay bundle: `tackle/limma_replay.py`, `tackle/rmd_replay.py`
