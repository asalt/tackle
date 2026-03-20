import pandas as pd


def test_sort_files_restricts_volcano_to_X_before_topn(tmp_path):
    from tackle.statfile_sorter import sort_files

    # X contains only the "g*" genes.
    x_geneids = [f"g{i}" for i in range(200)]
    X = pd.DataFrame(index=x_geneids, data={"dummy": [0] * len(x_geneids)})
    X.index = X.index.astype(str)

    # Volcano TSV contains many "h*" genes that are *more significant* than any "g*" gene.
    # Without restricting to X before selecting topn, the selection would pick only h* genes
    # and then intersect with X -> 0 rows.
    extra_geneids = [f"h{i}" for i in range(200)]
    df = pd.DataFrame(
        {
            "GeneID": extra_geneids + x_geneids,
            "log2_FC": [1.0] * (len(extra_geneids) + len(x_geneids)),
            "pValue": [1e-6 + i * 1e-9 for i in range(len(extra_geneids))]
            + [0.1 + i * 1e-4 for i in range(len(x_geneids))],
            "pAdj": [1e-6 + i * 1e-9 for i in range(len(extra_geneids))]
            + [0.1 + i * 1e-4 for i in range(len(x_geneids))],
        }
    )
    volcano_path = tmp_path / "volcano.tsv"
    df.to_csv(volcano_path, sep="\t", index=False)

    out = sort_files(
        [str(volcano_path)],
        X,
        sort_by="pValue",
        direction="up",
        topn=100,
        fc=0,
        pval_cutoff=1,
        pval_type="pAdj",
    )

    assert out.shape[0] == 100
    assert set(out.index).issubset(set(X.index))
    assert not any(str(g).startswith("h") for g in out.index)


def test_sort_and_select_topn_both_backfills_to_total_topn_when_possible():
    from tackle.statfile_sorter import sort_and_select_topn

    df = pd.DataFrame(
        {
            "GeneID": ["u1", "u2", "u3"]
            + [f"d{i}" for i in range(1, 11)]
            + ["z1"],
            "log2_FC": [1.0, 2.0, 0.5] + [-1.0] * 10 + [0.0],
            "pAdj": [0.5] * 14,
            "pValue": [0.02, 0.03, 0.04]
            + [0.001 * i for i in range(1, 11)]
            + [0.0],
        }
    )

    out = sort_and_select_topn(
        df,
        sort_by="pValue",
        direction="both",
        topn=10,
        fc=0,
        pval_cutoff=1,
        pval_type="pAdj",
    )

    assert out.shape[0] == 10
    assert out["GeneID"].astype(str).nunique() == 10
    # log2_FC == 0 should not be used for balanced up/down selection/backfill.
    assert "z1" not in set(out["GeneID"].astype(str))
