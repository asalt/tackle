import pytest


@pytest.mark.skip(reason="Planned column-z-score missing-value integration coverage")
def test_column_zscore_preserves_plot_mask_and_fills_clustering_matrix(
    ctx, cluster2_r_prereqs
):
    """Exercise both missing-value matrices through the full cluster2 dispatcher.

    TODO: add a narrow capture hook for the R result so this test can inspect the
    displayed and clustering matrices without changing normal user output.
    """
    # Arrange: replace part of the stub matrix with zeros/NA values and set the
    # corresponding data_obj.mask entries, leaving enough detections per sample.

    # Act: run cluster2 with z_score="col", z_score_fillna=True,
    # show_missing_values=True, and cluster_fillna="min".

    # Assert: the displayed matrix restores NA at every originally masked cell,
    # while each observed sample column is centered in column-z-score space.

    # Assert: the separate matrix used for row/column clustering is finite at the
    # masked cells because myzscore filled them before calculating distances.

    # Assert: the written artifact uses the canonical z_column filename token.
    raise NotImplementedError
