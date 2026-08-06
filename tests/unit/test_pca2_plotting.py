import pandas as pd
import pytest

from tackle.pca_plotting import resolve_pca2_figsize, resolve_pca2_scree_figsize


def test_pca2_figsize_uses_compact_base_for_an_ordinary_legend():
    metadata = pd.DataFrame({"group": ["A", "B", "A", "B"]})

    assert resolve_pca2_figsize(
        None,
        sample_metadata=metadata,
        color="group",
        marker=None,
    ) == (6.0, 7.0)


def test_pca2_figsize_expands_for_long_and_numerous_legend_entries():
    metadata = pd.DataFrame(
        {
            "group": [f"long_experimental_group_label_{i}" for i in range(20)],
            "batch": [f"batch_{i % 4}" for i in range(20)],
        }
    )

    width, height = resolve_pca2_figsize(
        None,
        sample_metadata=metadata,
        color="group",
        marker="batch",
    )

    assert width > 6.0
    assert height > 7.0


def test_pca2_figsize_preserves_an_explicit_override():
    metadata = pd.DataFrame({"group": [str(i) for i in range(30)]})

    assert resolve_pca2_figsize(
        (8.5, 5.5),
        sample_metadata=metadata,
        color="group",
        marker=None,
    ) == (8.5, 5.5)


@pytest.mark.parametrize(
    ("component_count", "expected_width"),
    ((6, 6.0), (8, 7.0), (10, 8.0), (15, 10.5), (20, 13.0)),
)
def test_pca2_scree_figsize_expands_with_component_count(
    component_count,
    expected_width,
):
    assert resolve_pca2_scree_figsize(
        (6.0, 7.0),
        component_count=component_count,
        auto_size=True,
    ) == (expected_width, 7.0)


def test_pca2_scree_figsize_preserves_explicit_size():
    assert resolve_pca2_scree_figsize(
        (7.25, 5.5),
        component_count=20,
        auto_size=False,
    ) == (7.25, 5.5)
