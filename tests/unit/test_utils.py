import matplotlib as mpl
import pandas as pd
import seaborn as sb

from tackle.utils import (
    get_default_color_mapping,
    get_color_mapping_with_defaults,
    normalize_metadata_str_values,
    parse_metadata,
    read_config,
)


def test_read_config_preserves_user_metadata_field_order(tmp_path):
    config_path = tmp_path / "ordered.conf"
    config_path.write_text(
        """\
[sample_a]
recno = 12345
runno = 1
searchno = 4
assay = DIA
tissue = brain
fraction = lEV
treatment = BR
group = lEV_BR
replicate = R1
""",
        encoding="utf-8",
    )

    config = read_config(str(config_path))
    metadata = parse_metadata(config)

    assert list(metadata.columns) == [
        "runno",
        "searchno",
        "label",
        "recno",
        "assay",
        "tissue",
        "fraction",
        "treatment",
        "group",
        "replicate",
    ]


def test_get_default_color_mapping_integer_like_series_uses_light_palette():
    series = pd.Series(["1", "2", "3"])

    mapping = get_default_color_mapping(series)

    numeric_series = pd.to_numeric(series, errors="coerce")
    unique_values = sorted(numeric_series.unique())
    expected_colors = [
        mpl.colors.rgb2hex(color)
        for color in sb.color_palette(
            palette="light:#4133", n_colors=len(unique_values)
        )
    ]
    expected_mapping = {
        str(value): color for value, color in zip(unique_values, expected_colors)
    }

    assert mapping == expected_mapping


def test_get_default_color_mapping_non_integer_series_uses_bright_palette():
    series = pd.Series(["apple", "banana", "cherry"])

    mapping = get_default_color_mapping(series)

    unique_values = sorted(series.unique())
    expected_colors = [
        mpl.colors.rgb2hex(color)
        for color in sb.color_palette(palette="bright", n_colors=len(unique_values))
    ]
    expected_mapping = {
        str(value): color for value, color in zip(unique_values, expected_colors)
    }

    assert mapping == expected_mapping


def test_get_default_color_mapping_applies_idg_reserved_colors():
    series = pd.Series(["Tbio", "Tchem", "Tclin", "Tdark"], name="IDG")

    mapping = get_default_color_mapping(series)

    assert mapping["Tbio"] == "blue"
    assert mapping["Tchem"] == "green"
    assert mapping["Tclin"] == "pink"
    assert mapping["Tdark"] == "black"


def test_get_default_color_mapping_applies_both_reserved_color():
    series = pd.Series(["CYTOPLASM", "BOTH"], name="CYTO_NUC")

    mapping = get_default_color_mapping(series)

    assert mapping["BOTH"] == "#984ea3"


def test_normalize_metadata_str_values_coerces_missing_and_blank_to_NA():
    series = pd.Series([None, "", "null", pd.NA, float("nan"), "true", "FALSE"])
    out = normalize_metadata_str_values(series)
    assert out.tolist() == ["NA", "NA", "NA", "NA", "NA", "True", "False"]


def test_normalize_metadata_str_values_preserves_none_category():
    series = pd.Series(["none", "None", None])
    out = normalize_metadata_str_values(series)
    assert out.tolist() == ["none", "None", "NA"]


def test_get_color_mapping_with_defaults_merges_user_mapping_keys_as_strings():
    series = pd.Series([1, 2, 3])
    mapping = get_color_mapping_with_defaults(series, overrides={"1": "#111111"})

    assert list(mapping.keys()) == ["1", "2", "3"]
    assert mapping["1"] == "#111111"
    assert mapping["2"].startswith("#")
    assert mapping["3"].startswith("#")
