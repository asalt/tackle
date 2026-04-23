import pytest

from tackle.containers import Data


def _make_data(**overrides):
    data = Data.__new__(Data)
    data.ifot = False
    data.ifot_ki = False
    data.ifot_tf = False
    data.median = False
    data.trim_mean = False
    data.quantile75 = False
    data.quantile90 = False
    data.genefile_norm = None
    for key, value in overrides.items():
        setattr(data, key, value)
    return data


def test_export_ibaq_column_name_none():
    data = _make_data()
    assert data._export_ibaq_column_name() == "iBAQ_dstrAdj"


@pytest.mark.parametrize(
    "flag,expected",
    [
        ("median", "iBAQ_dstrAdj_MED"),
        ("ifot", "iBAQ_dstrAdj_FOT"),
        ("trim_mean", "iBAQ_dstrAdj_TMN"),
        ("ifot_ki", "FOT_KI"),
        ("ifot_tf", "FOT_TF"),
        ("quantile75", "iBAQ_dstrAdj_Q75"),
        ("quantile90", "iBAQ_dstrAdj_Q90"),
    ],
)
def test_export_ibaq_column_name_flags(flag, expected):
    data = _make_data(**{flag: True})
    assert data._export_ibaq_column_name() == expected


def test_export_ibaq_column_name_genefile_norm():
    data = _make_data(genefile_norm="genes.txt")
    assert data._export_ibaq_column_name() == "iBAQ_dstrAdj_GNORM"


def test_mspc_ibaq_export_cols_include_core_columns():
    data = _make_data()
    assert data._mspc_ibaq_export_cols() == [
        "iBAQ_dstrAdj",
        "iBAQ_dstrAdj_nonorm",
        "iBAQ_dstrAdj_MED",
        "iBAQ_dstrAdj_FOT",
    ]


@pytest.mark.parametrize(
    "overrides,extra_col",
    [
        ({"ifot_ki": True}, "FOT_KI"),
        ({"ifot_tf": True}, "FOT_TF"),
        ({"trim_mean": True}, "iBAQ_dstrAdj_TMN"),
        ({"quantile75": True}, "iBAQ_dstrAdj_Q75"),
        ({"quantile90": True}, "iBAQ_dstrAdj_Q90"),
        ({"genefile_norm": "genes.txt"}, "iBAQ_dstrAdj_GNORM"),
    ],
)
def test_mspc_ibaq_export_cols_include_active_extra_norm(overrides, extra_col):
    data = _make_data(**overrides)
    cols = data._mspc_ibaq_export_cols()
    assert cols[:4] == [
        "iBAQ_dstrAdj",
        "iBAQ_dstrAdj_nonorm",
        "iBAQ_dstrAdj_MED",
        "iBAQ_dstrAdj_FOT",
    ]
    assert cols[-1] == extra_col


@pytest.mark.parametrize(
    "label,expected",
    [
        ("none", "0"),
        ("None", "0"),
        ("NONE", "0"),
        ("", ""),
        ("TMT126", "TMT126"),
    ],
)
def test_export_label_token(label, expected):
    assert Data._export_label_token(label) == expected
