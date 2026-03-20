import os
import re
import numpy as np
import pandas as pd

#############################################
# 1) DICTIONARY: SORT PRESETS
#############################################

# Each key in this dict is a valid `sort_by` option.
# For each, we can define:
#   - sort_func: a callable or a string column name
#   - ascending: bool or list[bool]
#   - up_query:  query string for "up" subset, or None if not applicable
#   - dn_query:  query string for "down" subset, or None if not applicable
#   - reverse_dn: whether to reverse the down subset after sorting (useful for certain workflows)
#
# You can add or remove columns easily. For "abs_signedlogp", we show how to compute the column on the fly.
sort_presets = {
    "pValue": {
        "sort_func": "pValue",   # just sort by the column "pValue"
        "ascending": True,       # smaller pValue => more significant
        "up_query": "log2_FC > 0",
        "dn_query": "log2_FC < 0",
        "reverse_dn": True,
    },
    "pAdj": {
        "sort_func": "pAdj",   # adjusted p-value
        "ascending": True,
        "up_query": "log2_FC > 0",
        "dn_query": "log2_FC < 0",
        "reverse_dn": True,
    },
    "log2_FC": {
        "sort_func": lambda df: (df["log2_FC"] * 1/(df["pValue"]+1e-9)),  # just sort by the column "log2_FC"
        "ascending": False,      # bigger FC => more interesting
        "up_query": "log2_FC > 0",
        "dn_query": "log2_FC < 0",
        "reverse_dn": False,
    },
    "signedlogP": {
        #"sort_func": "signedlogP",  # sort by "signedlogp" column
        "sort_func": lambda df: df["signedlogP"].abs(),
        "ascending": False,         # typically bigger is more significant
        "up_query": "signedlogP > 0",
        "dn_query": "signedlogP < 0",
        "reverse_dn": False,
    },
    "abs_signedlogp": {
        # We can define a callable function for the sort key:
        "sort_func": lambda df: df["signedlogP"].abs(),
        "ascending": False,
        # If up/down is meaningless here, just set them to None
        "up_query": None,
        "dn_query": None,
        "reverse_dn": False,
    },
}

#############################################
# 2) SORT + TOP-N SELECTION FUNCTION
#############################################

def sort_and_select_topn(
    df: pd.DataFrame,
    sort_by: str = "pValue",
    direction: str = "both",   # 'both', 'up', or 'down'
    topn: int = 50,
    fc=2,
    pval_cutoff=.05,
    pval_type="pValue"

) -> pd.DataFrame:
    """
    Sort the DataFrame according to the selected 'sort_by' preset,
    then pick the top N for up, down, or both.

    :param df:        Input DataFrame.
    :param sort_by:   Key in the `sort_presets` dict.
    :param direction: 'both', 'up', or 'down'
    :param topn:      Number of entries to keep (split in half if direction=='both').
    :return:          A new DataFrame with top (up, down, or both) rows.
    """
    
    # If not recognized, fallback or raise
    if sort_by not in sort_presets:
        print(f"[WARNING] sort_by='{sort_by}' not recognized; using 'pValue' as default.")
        sort_by = "pValue"
    
    preset = sort_presets[sort_by]

    # Apply optional fold-change / p-value cutoffs. The CLI defaults (fc=0, pval_cutoff=1)
    # are intended to mean "no cutoff".
    try:
        fc_val = float(fc)
    except Exception:
        fc_val = None
    # Interpret `fc` as a fold-change ratio (e.g. 2 => abs(log2_FC) > log2(2) == 1).
    if fc_val is not None and fc_val > 1 and "log2_FC" in df.columns:
        df = df[abs(df["log2_FC"]) > np.log2(fc_val)]

    try:
        pval_val = float(pval_cutoff)
    except Exception:
        pval_val = None
    if pval_val is not None and pval_val < 1:
        df = df[df[pval_type] < pval_val]
    
    # Retrieve the sort key: can be a string (column name) or a callable
    sort_key = preset["sort_func"]
    ascending = preset["ascending"]
    up_query = preset["up_query"]
    dn_query = preset["dn_query"]
    reverse_dn = preset["reverse_dn"]
    
    def sort_df(df_sub: pd.DataFrame) -> pd.DataFrame:
        if callable(sort_key):
            df_sub = df_sub.assign(__sortkey=sort_key(df_sub))
            sort_cols = ["__sortkey"]
        else:
            sort_cols = [sort_key]
        return df_sub.sort_values(by=sort_cols, ascending=ascending)

    # If direction is "both", we pick topn//2 from up and topn//2 from down and
    # then backfill (from the remaining best rows) to reach `topn` when possible.
    # Otherwise we pick topn from just up or down.
    if direction == "both":
        n_each = int(int(topn) // 2)
    else:
        n_each = topn
    
    def subset_and_sort(df_sub: pd.DataFrame, query_expr: str | None) -> pd.DataFrame:
        """
        Subset by a query_expr if provided, then sort, then pick top n_each rows.
        """
        if query_expr:
            df_sub = df_sub.query(query_expr)

        df_sub = sort_df(df_sub).head(n_each)
        
        # Some workflows like to reverse the 'down' group after sorting
        # (e.g. so the "most negative" is last). You can also skip this step.
        return df_sub
    
    # Re-compute subsets for up/down/both
    if direction == "both":
        # Treat log2_FC == 0 as neither up nor down for balanced selection/backfill.
        df_dir = df.query("log2_FC != 0") if "log2_FC" in df.columns else df

        up_part = subset_and_sort(df_dir, up_query) if up_query else df_dir.head(n_each)
        dn_part = subset_and_sort(df_dir, dn_query) if dn_query else df_dir.tail(n_each)
        if reverse_dn:
            dn_part = dn_part.iloc[::-1]  # Reverse
        combined = pd.concat([up_part, dn_part], axis=0)

        # Backfill to hit the requested total `topn` whenever possible.
        missing = int(topn) - int(len(combined))
        if missing > 0 and not df_dir.empty:
            selected = set()
            if "GeneID" in combined.columns:
                selected = set(combined["GeneID"].astype(str).tolist())

            candidates = sort_df(df_dir)
            if selected and "GeneID" in candidates.columns:
                candidates = candidates[~candidates["GeneID"].astype(str).isin(selected)]
            filler = candidates.head(missing)
            if not filler.empty:
                combined = pd.concat([combined, filler], axis=0)

        return combined
    
    elif direction == "up":
        return subset_and_sort(df, up_query)
    
    elif direction == "down":
        dn_part = subset_and_sort(df, dn_query)
        if reverse_dn:
            dn_part = dn_part.iloc[::-1]
        return dn_part
    
    else:
        raise ValueError(f"Unknown direction '{direction}'. Must be 'both', 'up', or 'down'.")


#############################################
# 3) SAMPLE "MAIN" LOGIC / WRAPPER
#############################################

def parse_volcano_file(volcano_file: str) -> pd.DataFrame:
    """
    Simple file loader with minimal column renaming if needed.
    Adjust to your real logic.
    """
    read_kws = {"dtype": {"GeneID": str, "TaxonID": str}}
    df = pd.read_table(volcano_file,
            **read_kws
            )

    # Example: fix column name
    if "p-value" in df.columns and "pValue" not in df.columns:
        df = df.rename(columns={"p-value": "pValue"})
    # if "p-value" in df.columns and "pValue" not in df.columns:
    #     df = df.rename(columns={"p-value": "pValue"})
    
    return df


def extract_name_group(volcano_file: str) -> str:
    """
    Extract some name from the filename, e.g. removing 'Batch...'.
    """
    basename = os.path.basename(volcano_file)
    if "Batch" in basename:
        basename = basename[basename.find("Batch") + 12 :]
    match = re.search(r"(?<=group)[_]?(.*)(?=\.tsv)", basename)
    return match.group(1) if match else basename


def process_file(
    volcano_file: str,
    sort_by: str,
    direction: str,
    topn: int,
    fc=0,
    pval_cutoff=1,
    pval_type="pValue",
    restrict_gene_ids: set[str] | None = None,
) -> pd.DataFrame:
    """
    Wrapper that:
      1) Loads the file
      2) Optionally restricts to a gene universe
      2) Sorts/picks topN
      3) Returns the final DataFrame
    """
    df = parse_volcano_file(volcano_file)
    if restrict_gene_ids and "GeneID" in df.columns:
        universe = {str(x) for x in restrict_gene_ids if x is not None}
        if universe:
            df = df[df["GeneID"].astype(str).isin(universe)]
    df_filtered = sort_and_select_topn(df, sort_by=sort_by, direction=direction, topn=topn,
            fc=fc,
            pval_cutoff=pval_cutoff,
            pval_type=pval_type
            )
    return df_filtered


def sort_files(
    volcano_files: list[str],
    X: pd.DataFrame,
    sort_by: str = "signedlogP",   # default
    direction: str = "both",       # default
    topn: int = 50,                # default
    fc=2,
    pval_cutoff=.05,
    pval_type="pValue"
) -> pd.DataFrame:
    """
    Example function:
      - Loads multiple volcano files
      - For each, picks top rows
      - Combines them
      - Returns X filtered by the relevant GeneIDs
    """
    print(f"Sorting by {sort_by}")
    print(f"topn: {topn}")
    
    # Normalize to a list
    if isinstance(volcano_files, str):
        volcano_files = [volcano_files]
    
    dfs_combined = []
    gene_universe = set(X.index.astype(str))
    
    for vf in volcano_files:
        df_filtered = process_file(
            vf,
            sort_by=sort_by,
            direction=direction,
            topn=topn,
            fc=fc,
            pval_cutoff=pval_cutoff,
            pval_type=pval_type,
            restrict_gene_ids=gene_universe,
        )
        if not df_filtered.empty:
            dfs_combined.append(df_filtered)
    
    if dfs_combined:
        final_df = pd.concat(dfs_combined, axis=0, ignore_index=True)
        gene_set = None
        if "GeneID" in final_df.columns:
            # Filter X by the set of gene IDs
            gene_set = final_df["GeneID"].astype(str).unique()
        if final_df.index.name == "GeneID":
            # raise NotImpmentedError("")
            pass
            # gene_set = final_df.index.astype(str).unique()
        if gene_set is not None and len(gene_set) > 0:
            #import ipdb; ipdb.set_trace()
            keep_genes = [g for g in gene_set if g in X.index]
            X = X.loc[keep_genes]
    
    return X
