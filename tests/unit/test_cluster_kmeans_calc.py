import numpy as np
import pandas as pd

from tackle.clusterplot import calc_kmeans


def test_calc_kmeans_basic_shapes():
    rng = np.random.default_rng(0)
    samples = [f"s{i+1}" for i in range(6)]
    # Build three well-separated clusters of genes to avoid singletons
    gA = [f"gA{i+1}" for i in range(5)]
    gB = [f"gB{i+1}" for i in range(5)]
    gC = [f"gC{i+1}" for i in range(5)]
    A = rng.normal(loc=0.0, scale=0.2, size=(5, 6)) + 0.0
    B = rng.normal(loc=0.2, scale=0.2, size=(5, 6)) + 4.5
    C = rng.normal(loc=0.2, scale=0.2, size=(5, 6)) - 4.5
    data = pd.DataFrame(
        np.vstack([A, B, C]),
        index=(gA + gB + gC),
        columns=samples,
    )

    res = calc_kmeans(data, nclusters=3, seed=123)
    assert res["nclusters"] == 3
    clusters = res["clusters"]
    assert set(clusters.index) == set(data.index)
    assert clusters.nunique() == 3
    # silhouette scores should be available per gene
    sil = res["silhouette_scores"]
    assert hasattr(sil, "__len__")
    assert len(sil) == data.shape[0]
