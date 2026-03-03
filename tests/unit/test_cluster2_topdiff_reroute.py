import pytest

from tackle.clusterplot_dispatcher import (
    _sanitize_topdiff_contrast,
    extract_contrast_from_volcano_filename,
)


@pytest.mark.parametrize(
    "basename, expected",
    [
        ("group_test.tsv", "test"),
        ("run1_volcano_nz1_groupA_minus_groupB.tsv", "A_minus_groupB"),
        (
            "volcano_nz1_dir_B_fna_T_group_treatB_minus_treatA_imv_T_lrob_F_ltrd_F.tsv",
            "treatB_minus_treatA",
        ),
        ("volcano_output.tsv", "volcano_output"),
    ],
)
def test_extract_contrast_from_volcano_filename(basename, expected):
    assert extract_contrast_from_volcano_filename(basename) == expected


def test_sanitize_topdiff_contrast_replaces_problematic_characters():
    raw = "A:B+C?D|E/F\\G"
    assert _sanitize_topdiff_contrast(raw) == "A_B_CqmkDorE_F_G"
