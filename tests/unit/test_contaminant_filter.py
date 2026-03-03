import pandas as pd


def test_drop_geneid_prefix_drops_by_geneid_column():
    from tackle.containers import drop_geneid_prefix

    df = pd.DataFrame(
        {
            "GeneID": ["1", "Cont_foo", "contaminant123", "2"],
            "value": [10, 20, 30, 40],
        }
    )

    out = drop_geneid_prefix(df, "Cont", geneid_col="GeneID")
    assert out["GeneID"].tolist() == ["1", "2"]
    assert out["value"].tolist() == [10, 40]


def test_drop_geneid_prefix_drops_by_index_when_no_geneid_col():
    from tackle.containers import drop_geneid_prefix

    df = pd.DataFrame({"value": [1, 2, 3]}, index=["ContA", "x", "ContB"])
    out = drop_geneid_prefix(df, "Cont")

    assert out.index.tolist() == ["x"]
    assert out["value"].tolist() == [2]


def test_drop_geneid_prefix_inplace_mutates_dataframe():
    from tackle.containers import drop_geneid_prefix

    df = pd.DataFrame({"GeneID": ["ContA", "1"], "value": [10, 20]})
    drop_geneid_prefix(df, "cont", geneid_col="GeneID", inplace=True)

    assert df["GeneID"].tolist() == ["1"]
    assert df["value"].tolist() == [20]
