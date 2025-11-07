import os
from pathlib import Path

import pandas as pd

from tackle import statfile_sorter as sfs


def test_sort_files_with_multiple_volcano_files_filters_union(stub_data_obj, tmp_path):
    X = stub_data_obj.areas_log
    base = Path(tmp_path)
    vdir = base / "volcano"
    vdir.mkdir(parents=True, exist_ok=True)

    # Prepare two volcano files with distinct up/down sets
    genes = list(X.index)
    up1, dn1 = genes[:3], genes[3:6]
    up2, dn2 = genes[6:9], genes[9:12]

    df1 = pd.DataFrame(
        {
            "GeneID": up1 + dn1,
            "log2_FC": [3.0] * len(up1) + [-3.0] * len(dn1),
            "pAdj": [0.001] * 6,
            "pValue": [0.001] * 6,
        }
    )
    df2 = pd.DataFrame(
        {
            "GeneID": up2 + dn2,
            "log2_FC": [3.0] * len(up2) + [-3.0] * len(dn2),
            "pAdj": [0.001] * 6,
            "pValue": [0.001] * 6,
        }
    )
    f1 = vdir / "B_by_A.tsv"
    f2 = vdir / "C_by_A.tsv"
    df1.to_csv(f1, sep="\t", index=False)
    df2.to_csv(f2, sep="\t", index=False)

    Xf = sfs.sort_files(
        [str(f1), str(f2)],
        X.copy(),
        sort_by="pValue",
        direction="both",
        topn=6,
        fc=2.0,
        pval_cutoff=0.01,
        pval_type="pAdj",
    )

    expected = set(up1 + dn1 + up2 + dn2)
    assert set(Xf.index) == expected
    assert Xf.shape[1] == X.shape[1]

