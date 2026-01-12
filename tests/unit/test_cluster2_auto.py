import pandas as pd

from tackle.cluster2.auto import select_best_k, summarize_silhouette_df


def test_summarize_silhouette_df_handles_nans():
    df = pd.DataFrame(
        {
            "cluster": [1, 1, 2, 2],
            "sil_width": [0.5, 0.1, float("nan"), -0.2],
        }
    )
    summary = summarize_silhouette_df(df)
    assert summary["n_genes"] == 4
    assert summary["n_clusters"] == 2
    assert summary["sil_mean"] == (0.5 + 0.1 - 0.2) / 3
    expected_q10 = pd.Series([0.5, 0.1, -0.2]).quantile(0.10)
    assert summary["sil_q10"] == expected_q10
    assert summary["sil_neg_frac"] == 1 / 3


def test_select_best_k_prefers_mean_then_q10_then_neg_frac_then_smaller_k():
    candidates = [
        dict(nclusters=3, sil_mean=0.5, sil_q10=0.1, sil_neg_frac=0.2),
        dict(nclusters=4, sil_mean=0.6, sil_q10=0.05, sil_neg_frac=0.0),
        dict(nclusters=5, sil_mean=0.6, sil_q10=0.07, sil_neg_frac=0.4),
        dict(nclusters=6, sil_mean=0.6, sil_q10=0.07, sil_neg_frac=0.1),
        dict(nclusters=7, sil_mean=0.6, sil_q10=0.07, sil_neg_frac=0.1),
    ]
    assert select_best_k(candidates) == 6


def test_select_best_k_ignores_nan_scores():
    candidates = [
        dict(nclusters=2, sil_mean=float("nan"), sil_q10=float("nan"), sil_neg_frac=0.0),
        dict(nclusters=3, sil_mean=0.1, sil_q10=-0.5, sil_neg_frac=0.9),
    ]
    assert select_best_k(candidates) == 3
