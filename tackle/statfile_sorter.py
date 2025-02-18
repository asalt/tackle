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
    "log2_FC": {
        "sort_func": "log2_FC",  # just sort by the column "log2_FC"
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
    topn: int = 50
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
    
    # Retrieve the sort key: can be a string (column name) or a callable
    sort_key = preset["sort_func"]
    ascending = preset["ascending"]
    up_query = preset["up_query"]
    dn_query = preset["dn_query"]
    reverse_dn = preset["reverse_dn"]
    
    # If direction is "both", we pick topn//2 from up and topn//2 from down
    # Otherwise we pick topn from just up or down.
    if direction == "both":
        n_each = int(topn // 2)
    else:
        n_each = topn
    
    def subset_and_sort(df_sub: pd.DataFrame, query_expr: str | None) -> pd.DataFrame:
        """
        Subset by a query_expr if provided, then sort, then pick top n_each rows.
        """
        if query_expr:
            df_sub = df_sub.query(query_expr)
        
        if callable(sort_key):
            # Evaluate the function on df_sub to get a Series to sort by
            df_sub = df_sub.assign(__sortkey=sort_key(df_sub))
            sort_cols = ["__sortkey"]
        else:
            # Otherwise it's a column name
            sort_cols = [sort_key]
        
        df_sub = df_sub.sort_values(by=sort_cols, ascending=ascending).head(n_each)
        
        # Some workflows like to reverse the 'down' group after sorting
        # (e.g. so the "most negative" is last). You can also skip this step.
        return df_sub
    
    # Re-compute subsets for up/down/both
    if direction == "both":
        up_part = subset_and_sort(df, up_query) if up_query else df.head(n_each)
        dn_part = subset_and_sort(df, dn_query) if dn_query else df.tail(n_each)
        if reverse_dn:
            dn_part = dn_part.iloc[::-1]  # Reverse
        return pd.concat([up_part, dn_part], axis=0)
    
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
    topn: int
) -> pd.DataFrame:
    """
    Wrapper that:
      1) Loads the file
      2) Sorts/picks topN
      3) Returns the final DataFrame
    """
    df = parse_volcano_file(volcano_file)
    df_filtered = sort_and_select_topn(df, sort_by=sort_by, direction=direction, topn=topn)
    return df_filtered


def sort_files(
    volcano_files: list[str],
    X: pd.DataFrame,
    sort_by: str = "signedlogP",   # default
    direction: str = "both",       # default
    topn: int = 50                 # default
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
    
    for vf in volcano_files:
        df_filtered = process_file(vf, sort_by=sort_by, direction=direction, topn=topn)
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
