import pandas as pd

from tackle.containers import Data


def test_mspc_export_source_uses_every_loaded_geneid_and_returns_a_copy():
    data = Data.__new__(Data)
    data.data = pd.DataFrame(
        {
            "GeneID": ["g1", "g1", "g2", "g2"],
            "Metric": ["GeneSymbol", "area", "GeneSymbol", "area"],
            "S1": ["A", 10.0, "B", 0.0],
        }
    )
    data.df_filtered = (
        data.data.loc[data.data["GeneID"] == "g1"]
        .set_index(["GeneID", "Metric"])
    )

    source = data._mspc_export_source()

    assert list(source.index.get_level_values("GeneID").unique()) == ["g1", "g2"]
    source.loc[("g1", "area"), "S1"] = 99.0
    original = data.data.loc[
        (data.data["GeneID"] == "g1") & (data.data["Metric"] == "area"),
        "S1",
    ].iloc[0]
    assert original == 10.0
