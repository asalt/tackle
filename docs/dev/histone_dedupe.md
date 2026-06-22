# Histone exact-duplicate dedupe notes

The current histone dedupe step is intentionally not a biological sum or family collapse.
It handles the case where multiple histone GeneIDs have the exact same abundance vector,
which implies shared peptide evidence rather than distinguishable isoform measurements.

Implementation notes:

- The helper lives in `tackle/histone_dedupe.py`.
- The call surface is narrow: `self.data`, `HistoneInfo()`, and `self.outpath`.
- `self.data` is the long `GeneID`/`Metric` frame produced inside `containers.Data.load_data()`.
- Exact duplicate profiles are compared using `Metric == "AreaSum_dstrAdj"`.
- Dedupe is taxon-aware: identical vectors are only reduced within the same `TaxonID`.
- The selected row remains in downstream analyses and exports.
- Removed equivalent candidates are recorded in `<outpath>/context/histone_dedupe.tsv`.

Future cleanup targets:

- `containers.Data.load_data()` is still very long and mixes loading, filtering, normalization,
  metric accounting, row reduction, and output-prep concerns.
- Similar small helpers could be extracted for row-level cleanup steps, but histone dedupe should
  stay narrow because it is a data-shape-changing operation.
