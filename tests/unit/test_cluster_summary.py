from pathlib import Path

import pandas as pd

from tackle.cluster2.summary import summarize_cluster_files


def test_cluster_summary_tables(tmp_path: Path):
    cluster_file = tmp_path / "example_clusters.tsv"
    df = pd.DataFrame(
        {
            "GeneID": ["101", "102", "103", "104"],
            "GeneSymbol": ["A", "B", "C", "D"],
            "cluster": [1, 1, 2, 2],
            "sil_width": [0.1, -0.2, 0.3, 0.4],
        }
    )
    df.to_csv(cluster_file, sep="\t", index=False)

    tables = summarize_cluster_files([str(cluster_file)])
    runs = tables.runs
    clusters = tables.clusters

    assert len(runs) == 1
    assert runs.loc[0, "n_genes"] == 4
    assert runs.loc[0, "n_clusters"] == 2
    assert runs.loc[0, "sil_neg_frac"] == 0.25

    assert len(clusters) == 2
    c1 = clusters.query("cluster == 1").iloc[0]
    c2 = clusters.query("cluster == 2").iloc[0]
    assert c1["n_genes"] == 2
    assert c1["sil_neg_frac"] == 0.5
    assert c2["n_genes"] == 2
    assert c2["sil_neg_frac"] == 0.0
