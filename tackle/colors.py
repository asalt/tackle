from __future__ import annotations

# Stable discrete color mappings for common metadata/annotation values.
#
# These are passed through to R plotting, so values should be either R-recognized
# color names or hex strings.

RESERVED_VALUE_COLORS: dict[str, str] = {
    "True": "green",
    "False": "red",
    "NA": "grey",
    # Used by multiple annotation tracks (e.g. CYTO_NUC, ER_GOLGI).
    "BOTH": "#984ea3",
}

RESERVED_VALUE_COLORS_BY_FIELD: dict[str, dict[str, str]] = {
    "IDG": {
        "Tbio": "blue",
        "Tchem": "green",
        "Tclin": "pink",
        "Tdark": "black",
    },
    "CYTO_NUC": {
        "CYTOPLASM": "#1f78b4",
        "NUCLEUS": "#33a02c",
        "BOTH": RESERVED_VALUE_COLORS["BOTH"],
    },
    "ER_GOLGI": {
        "ENDORETICULUM": "#e31a1c",
        "GOLGI": "#ff7f00",
        "BOTH": RESERVED_VALUE_COLORS["BOTH"],
    },
    "MATRISOME": {
        "CORE": "#1b9e77",
        "ASSOCIATED": "#d95f02",
    },
    "SurfaceLabel": {
        "surface": "#e31a1c",
        "nonsurface": "#999999",
    },
    "SECRETED": {
        "SECRETED": "#4daf4a",
    },
    "IO": {
        "IO": "#a65628",
    },
    "glycomineN": {
        "N": "#66c2a5",
    },
    "glycomineO": {
        "O": "#fc8d62",
    },
    "CellMembrane": {
        "Cell Membrane": "#8da0cb",
    },
}

