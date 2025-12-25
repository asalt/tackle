import pytest

from tackle.volcanoplot import (
    _clean_group_label,
    _dedupe_comparison_label,
    _split_comparison_groups,
)


def test_dedupe_comparison_label_collapses_duplicate():
    comparison = (
        "PCL_Bi_MSC_minus_ctrl = (A+B)/2 - (C+D)/2="
        "PCL_Bi_MSC_minus_ctrl = (A+B)/2 - (C+D)/2"
    )

    assert (
        _dedupe_comparison_label(comparison)
        == "PCL_Bi_MSC_minus_ctrl = (A+B)/2 - (C+D)/2"
    )


@pytest.mark.parametrize(
    "comparison",
    ["A - B", "A-B", "A - B-C", "(A-B) - (C-D)", "A-1 - B"],
)
def test_split_comparison_groups_simple(comparison):
    if comparison == "A - B-C":
        expected = ("A", "B-C")
    elif comparison == "(A-B) - (C-D)":
        expected = ("(A-B)", "(C-D)")
    elif comparison == "A-1 - B":
        expected = ("A-1", "B")
    else:
        expected = ("A", "B")
    assert _split_comparison_groups(comparison) == expected


@pytest.mark.parametrize(
    "comparison",
    ["A - B - C", "A-B-C"],
)
def test_split_comparison_groups_ambiguous(comparison):
    assert _split_comparison_groups(comparison) is None


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("(A+B)", "A+B"),
        (" (A+B)/2 ", "(A+B)/2"),
    ],
)
def test_clean_group_label_parentheses(value, expected):
    assert _clean_group_label(value) == expected
