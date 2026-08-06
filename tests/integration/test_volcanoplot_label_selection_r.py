from pathlib import Path

import pytest

pytest.importorskip("rpy2")
from rpy2 import robjects


@pytest.fixture(scope="module")
def select_labels():
    r_file = Path(__file__).parents[2] / "tackle" / "R" / "volcanoplot.R"
    robjects.r["source"](str(r_file))
    return robjects.globalenv[".select_volcano_labels"]


@pytest.fixture(scope="module")
def density_label_sizes(select_labels):
    del select_labels
    return robjects.globalenv[".volcano_density_label_sizes"]


@pytest.fixture(scope="module")
def label_candidates():
    return robjects.r(
        """
        data.frame(
          log2_FC = c(-8, -6, -4, -2, 1, 3, 5, 7),
          pValue = c(.08, .07, .06, .01, .02, .03, .04, .05),
          row.names = c('L8', 'L6', 'L4', 'L2', 'R1', 'R3', 'R5', 'R7')
        )
        """
    )


def _selected(select_labels, label_candidates, **kwargs):
    result = select_labels(label_candidates, **kwargs)
    return list(result)


def test_no_side_override_preserves_balanced_log2_selection(
    select_labels, label_candidates
):
    assert _selected(
        select_labels,
        label_candidates,
        max_labels=4,
        number_by="log2_FC",
        direction="both",
    ) == ["L8", "L6", "R7", "R5"]


def test_left_override_replaces_only_left_selection(select_labels, label_candidates):
    assert _selected(
        select_labels,
        label_candidates,
        max_labels=4,
        max_labels_left=1,
        number_by="log2_FC",
        direction="both",
    ) == ["L8", "R7", "R5", "R3"]


def test_right_override_replaces_only_right_selection(select_labels, label_candidates):
    assert _selected(
        select_labels,
        label_candidates,
        max_labels=4,
        max_labels_right=1,
        number_by="log2_FC",
        direction="both",
    ) == ["L8", "L6", "L4", "R7"]


def test_zero_side_override_suppresses_that_side(select_labels, label_candidates):
    assert _selected(
        select_labels,
        label_candidates,
        max_labels=4,
        max_labels_left=0,
        number_by="log2_FC",
        direction="both",
    ) == ["R7", "R5", "R3", "R1"]


def test_side_override_uses_requested_ranking_metric(select_labels, label_candidates):
    assert _selected(
        select_labels,
        label_candidates,
        max_labels=4,
        max_labels_left=1,
        max_labels_right=1,
        number_by="pValue",
        direction="both",
    ) == ["L2", "R1"]


def test_two_side_overrides_ignore_overall_number(select_labels, label_candidates):
    assert _selected(
        select_labels,
        label_candidates,
        max_labels=1,
        max_labels_left=2,
        max_labels_right=3,
        number_by="log2_FC",
        direction="both",
    ) == ["L8", "L6", "R7", "R5", "R3"]


def test_hybrid_selection_covers_fc_and_pvalue_extremes_on_each_side(
    select_labels, label_candidates
):
    selected = _selected(
        select_labels,
        label_candidates,
        max_labels=4,
        number_by="hybrid",
        direction="both",
    )

    assert len(selected) == 4
    assert set(selected) == {"L8", "L2", "R7", "R1"}


def test_hybrid_selection_uses_one_budget_without_duplicate_labels(select_labels):
    candidates = robjects.r(
        """
        data.frame(
          log2_FC = c(10, 9, 2, 1),
          pValue = c(.20, .15, 1e-6, 1e-5),
          row.names = c('FC_TOP', 'FC_NEXT', 'P_TOP', 'P_NEXT')
        )
        """
    )

    selected = _selected(
        select_labels,
        candidates,
        max_labels=2,
        number_by="hybrid",
        direction="up",
    )

    assert len(selected) == 2
    assert set(selected) == {"P_TOP", "FC_TOP"}


def test_density_label_sizes_make_isolated_labels_larger(density_label_sizes):
    sizes = list(
        density_label_sizes(
            robjects.FloatVector([0.00, 0.02, 0.04, 1.00]),
            robjects.FloatVector([0.00, 0.02, 0.04, 1.00]),
            min_size=2.4,
            max_size=4.0,
        )
    )

    assert min(sizes) == pytest.approx(2.4)
    assert max(sizes) == pytest.approx(4.0)
    assert sizes[-1] > max(sizes[:-1])


def test_density_label_sizes_use_maximum_for_one_label(density_label_sizes):
    sizes = list(
        density_label_sizes(
            robjects.FloatVector([0.5]),
            robjects.FloatVector([0.5]),
            min_size=2.4,
            max_size=4.0,
        )
    )

    assert sizes == [4.0]
