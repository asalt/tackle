import numpy as np
import pandas as pd

from tackle.zscore import my_zscore


def test_my_zscore_uses_absence_to_anchor_sparse_constant_detections():
    values = pd.Series([np.nan, 10.0, 10.0, np.nan], index=list("ABCD"))

    result = my_zscore(values)

    assert result[["A", "D"]].isna().all()
    assert np.allclose(result[["B", "C"]], np.sqrt(3) / 2)


def test_my_zscore_handles_singletons_constants_and_all_missing():
    singleton = my_zscore(pd.Series([np.nan, 20.0, np.nan, np.nan]))
    constant = my_zscore(pd.Series([5.0, 5.0, 5.0]))
    all_missing = my_zscore(pd.Series([np.nan, np.nan]))

    assert np.isclose(singleton.iloc[1], 1.5)
    assert singleton.drop(index=1).isna().all()
    assert np.allclose(constant, 0.0)
    assert all_missing.isna().all()


def test_containers_reexports_shared_my_zscore_for_zscore_export():
    from tackle.containers import my_zscore as export_my_zscore

    assert export_my_zscore is my_zscore


def test_zscore_export_uses_logged_masked_matrix_and_shared_transform(
    tmp_path, monkeypatch
):
    from tackle import containers

    matrix = pd.DataFrame(
        {
            "S1": [1.0, np.nan],
            "S2": [3.0, 20.0],
            "S3": [np.nan, np.nan],
        },
        index=pd.Index(["g1", "g2"], name="GeneID"),
    )
    mask = matrix.isna()
    data = containers.Data.__new__(containers.Data)
    data._areas_log = matrix.copy()
    data._mask = mask.copy()
    data._export_ibaq_column_name = lambda: "area"
    data.ifot = False
    data.ifot_ki = False
    data.ifot_tf = False
    data.median = True
    data.trim_mean = False
    data.quantile75 = False
    data.quantile90 = False
    data.outpath_name = "test"
    data.taxon = "all"
    data.non_zeros = 1
    data.colors_only = False
    data.batch_applied = None
    data.batch_nonparametric = False
    data.outpath = str(tmp_path)
    monkeypatch.setattr(
        containers,
        "get_outname",
        lambda *_args, **_kwargs: str(tmp_path / "zscore_export"),
    )

    output = data._perform_data_export(level="zscore", force=True)
    exported = pd.read_table(output, dtype={"GeneID": str})
    expected = matrix.apply(my_zscore, axis=1)

    assert list(exported.columns) == ["GeneID", "Metric", "S1", "S2", "S3"]
    assert exported["Metric"].eq("zscore").all()
    assert np.allclose(
        exported[["S1", "S2", "S3"]],
        expected,
        equal_nan=True,
    )
