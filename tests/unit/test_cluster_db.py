from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

from tackle.cluster2.db import ClusterDb, ClusterDbConfig


def test_cluster_db_ingest_and_ari(tmp_path: Path):
    analysis_outpath = tmp_path / "analysis"
    analysis_outpath.mkdir(parents=True, exist_ok=True)

    db_path = analysis_outpath / "cluster2" / "cluster2.sqlite"
    cfg = ClusterDbConfig(enabled=True, db_path=str(db_path))
    db = ClusterDb(cfg)

    tsv1 = analysis_outpath / "clustermap" / "tx" / "1more_kmeans_2.tsv"
    df1 = pd.DataFrame(
        {
            "GeneID": ["g1", "g2", "g3", "g4"],
            "GeneSymbol": ["A", "B", "C", "D"],
            "cluster": [1, 1, 2, 2],
            "neighbor": [2, 2, 1, 1],
            "sil_width": [0.2, 0.1, -0.1, 0.4],
        }
    )
    run1 = db.ingest_cluster2_metrics(
        df1,
        analysis_outpath=str(analysis_outpath),
        tsv_path=str(tsv1),
        cluster_func="kmeans",
        nclusters=2,
        seed=1234,
        linkage="ward.D2",
        cluster_fillna="min",
        z_score="0",
        z_score_by=None,
        standard_scale="None",
    )
    assert run1 is not None

    tsv2 = analysis_outpath / "clustermap" / "tx" / "1more_kmeans_3.tsv"
    df2 = pd.DataFrame(
        {
            "GeneID": ["g1", "g2", "g3", "g4"],
            "GeneSymbol": ["A", "B", "C", "D"],
            "cluster": [1, 2, 2, 3],
            "neighbor": [2, 1, 1, 2],
            "sil_width": [0.3, 0.0, 0.2, 0.1],
        }
    )
    run2 = db.ingest_cluster2_metrics(
        df2,
        analysis_outpath=str(analysis_outpath),
        tsv_path=str(tsv2),
        cluster_func="kmeans",
        nclusters=3,
        seed=1234,
        linkage="ward.D2",
        cluster_fillna="min",
        z_score="0",
        z_score_by=None,
        standard_scale="None",
    )
    assert run2 is not None

    db.close()

    # Validate DB contents
    import sqlite3

    conn = sqlite3.connect(str(db_path))
    try:
        runs = conn.execute("SELECT COUNT(*) FROM runs").fetchone()[0]
        memberships = conn.execute("SELECT COUNT(*) FROM memberships").fetchone()[0]
        clusters = conn.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]
        ari_rows = conn.execute("SELECT COUNT(*) FROM ari").fetchone()[0]
        assert runs == 2
        assert memberships == 8
        assert clusters >= 5  # 2 clusters + 3 clusters
        assert ari_rows == 1

        (ari_val,) = conn.execute("SELECT ari FROM ari").fetchone()
        expected = float(adjusted_rand_score(df1["cluster"].tolist(), df2["cluster"].tolist()))
        assert float(ari_val) == expected
    finally:
        conn.close()
