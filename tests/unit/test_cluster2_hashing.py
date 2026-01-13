import pandas as pd

from tackle.cluster2.hashing import hash_cluster2_input_matrix


def test_hash_cluster2_input_matrix_ignores_gene_symbol():
    df = pd.DataFrame(
        {
            "GeneID": ["g1", "g2"],
            "GeneSymbol": ["A", "B"],
            "s1": [1.0, 2.0],
            "s2": [3.0, 4.0],
        }
    )
    h1 = hash_cluster2_input_matrix(df)
    df2 = df.copy()
    df2["GeneSymbol"] = ["X", "Y"]
    h2 = hash_cluster2_input_matrix(df2)
    assert h1 == h2


def test_hash_cluster2_input_matrix_changes_when_values_change():
    df = pd.DataFrame(
        {
            "GeneID": ["g1", "g2"],
            "GeneSymbol": ["A", "B"],
            "s1": [1.0, 2.0],
            "s2": [3.0, 4.0],
        }
    )
    h1 = hash_cluster2_input_matrix(df)
    df2 = df.copy()
    df2.loc[0, "s1"] = 9.0
    h2 = hash_cluster2_input_matrix(df2)
    assert h1 != h2

