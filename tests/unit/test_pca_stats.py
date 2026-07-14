import numpy as np
import pandas as pd

from tackle.pca_stats import (
    adjust_pvalues,
    analyze_pca_separation,
    euclidean_r_squared,
    format_pca_test_caption,
    pairwise_centroid_geometry,
    resolve_pca_test_scopes,
    welch_james_test,
)


def _reference_scores_and_groups():
    scores = np.array(
        [
            [0, 1],
            [1, 0],
            [2, 2],
            [1, 3],
            [3, 1],
            [2, 4],
            [4, 2],
            [3, 5],
            [5, 3],
            [6, 6],
            [4, 7],
            [6, 8],
            [7, 6],
            [9, 9],
            [8, 7],
            [10, 11],
            [7, 10],
            [11, 12],
        ],
        dtype=float,
    )
    groups = np.array(["A"] * 5 + ["B"] * 6 + ["C"] * 7, dtype=object)
    return scores, groups


def test_welch_james_matches_welchadf_reference_result():
    scores, groups = _reference_scores_and_groups()

    result = welch_james_test(scores, groups)

    assert result.status == "ok"
    assert np.isclose(result.statistic, 18.525774053380193)
    assert result.numerator_df == 4
    assert np.isclose(result.denominator_df, 9.711665712288147)
    assert np.isclose(result.p_value, 0.0001508437775382289)


def test_plane_r2_is_direct_euclidean_centroid_partition():
    scores = np.array([[0, 0], [2, 0], [8, 0], [10, 0]], dtype=float)
    groups = np.array(["A", "A", "B", "B"], dtype=object)

    # Grand-centered SS = 68; within-group SS = 4.
    assert np.isclose(euclidean_r_squared(scores, groups), 1 - 4 / 68)


def test_pairwise_centroid_geometry_uses_rms_within_group_radii():
    scores = np.array([[0, 0], [2, 0], [8, 0], [10, 0]], dtype=float)
    groups = np.array(["A", "A", "B", "B"], dtype=object)

    geometry = pairwise_centroid_geometry(scores, groups)

    assert geometry["geometry_group_a"] == "A"
    assert geometry["geometry_group_b"] == "B"
    assert np.isclose(geometry["centroid_distance"], 8.0)
    assert np.isclose(geometry["rms_radius_a"], 1.0)
    assert np.isclose(geometry["rms_radius_b"], 1.0)
    assert np.isclose(geometry["pooled_rms_radius"], 1.0)
    assert np.isclose(geometry["standardized_separation"], 8.0)


def test_default_scopes_include_displayed_planes_and_leading_candidates():
    scores = pd.DataFrame(
        {
            "PC1": [-2.0, -1.0, 1.0, 2.0],
            "PC2": [1.0, -1.0, -1.0, 1.0],
            "PC3": [0.5, -0.5, 0.25, -0.25],
            "PC4": [0.0, 0.0, 0.0, 0.0],
        }
    )

    scopes = resolve_pca_test_scopes((), scores=scores, max_pc=3)

    assert [scope.pcs for scope in scopes] == [
        ("PC1", "PC2"),
        ("PC1", "PC3"),
        ("PC2", "PC3"),
        ("PC1", "PC2", "PC3", "PC4"),
    ]
    assert scopes[-1].name == "leading_estimable_pcs"
    assert scopes[-1].selection == "leading_estimable"


def test_singular_full_scope_is_reported_without_pseudoinverse():
    scores = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 2.0],
            [2.0, 2.0, 4.0],
            [4.0, 0.0, 4.0],
            [5.0, 1.0, 6.0],
            [6.0, 2.0, 8.0],
        ]
    )
    groups = np.array(["A"] * 3 + ["B"] * 3, dtype=object)

    result = welch_james_test(scores, groups)

    assert result.status == "singular_covariance"
    assert np.isnan(result.statistic)


def test_pvalue_adjustments_match_standard_stepwise_definitions():
    pvalues = [0.01, 0.04, 0.03]

    assert np.allclose(adjust_pvalues(pvalues, "holm"), [0.03, 0.06, 0.06])
    assert np.allclose(adjust_pvalues(pvalues, "hochberg"), [0.03, 0.04, 0.04])
    assert np.allclose(adjust_pvalues(pvalues, "BH"), [0.03, 0.04, 0.04])


def test_analysis_writes_omnibus_and_adjusted_pairwise_families():
    values, groups = _reference_scores_and_groups()
    scores = pd.DataFrame(values, columns=["PC1", "PC2"])
    scores["PC3"] = np.linspace(-1.0, 1.0, len(scores))
    scores.index = [f"S{index:02d}" for index in range(len(scores))]
    metadata = pd.DataFrame({"group": groups}, index=scores.index)
    scopes = resolve_pca_test_scopes(
        ("1,2", "1,3"), scores=scores, max_pc=3
    )

    omnibus, pairwise = analyze_pca_separation(
        scores,
        metadata,
        group_fields=("group",),
        scopes=scopes,
        p_adjust_method="holm",
    )

    assert list(omnibus["scope"]) == ["PC1_PC2", "PC1_PC3"]
    assert omnibus["p_adj"].notna().all()
    assert omnibus["explained_variance_pct"].between(0, 100).all()
    assert len(pairwise) == 6
    assert pairwise["explained_variance_pct"].between(0, 100).all()
    assert pairwise["p_adj"].notna().all()
    assert pairwise["p_adj_all_scopes"].notna().all()
    assert (pairwise["p_adj_all_scopes"] >= pairwise["p_value"]).all()


def test_leading_scope_uses_largest_dimension_estimable_for_all_tests():
    rng = np.random.default_rng(20260713)
    scores = pd.DataFrame(
        rng.normal(size=(12, 12)),
        columns=[f"PC{index}" for index in range(1, 13)],
        index=[f"S{index:02d}" for index in range(12)],
    )
    metadata = pd.DataFrame(
        {"group": np.repeat(["A", "B", "C", "D"], 3)},
        index=scores.index,
    )
    scopes = resolve_pca_test_scopes(
        ("displayed", "leading"), scores=scores, max_pc=2
    )

    omnibus, pairwise = analyze_pca_separation(
        scores,
        metadata,
        group_fields=("group",),
        scopes=scopes,
        p_adjust_method="holm",
    )

    assert len(omnibus) == 1
    assert omnibus.iloc[0]["scope"] == "PC1_PC2"
    assert "is_leading_scope" not in omnibus.columns
    assert omnibus.iloc[0]["n_pcs"] == 2
    expected_variance = (
        scores[["PC1", "PC2"]].var().sum() / scores.var().sum() * 100
    )
    assert np.isclose(
        omnibus.iloc[0]["explained_variance_pct"], expected_variance
    )
    assert omnibus.iloc[0]["status"] == "ok"
    assert len(pairwise) == 6
    assert "is_leading_scope" not in pairwise.columns
    assert pairwise["n_pcs"].eq(2).all()
    assert pairwise["status"].eq("ok").all()


def test_caption_names_the_selected_pvalue_adjustment():
    row = {
        "group_field": "group",
        "r2": 0.625,
        "status": "ok",
        "numerator_df": 4.0,
        "denominator_df": 7.5,
        "welch_james_f": 3.25,
        "p_adj": 0.0123,
        "p_adjust_method": "holm",
    }

    caption = format_pca_test_caption(pd.DataFrame([row]))
    assert "Holm-adjusted p=0.0123" in caption

    row["p_adjust_method"] = "none"
    caption = format_pca_test_caption(pd.DataFrame([row]))
    assert "; p=0.0123" in caption
    assert "adjusted p" not in caption

    row.update(
        {
            "geometry_group_a": "A",
            "geometry_group_b": "B",
            "centroid_distance": 0.63,
            "rms_radius_a": 0.18,
            "rms_radius_b": 0.14,
            "pooled_rms_radius": np.sqrt((0.18**2 + 0.14**2) / 2),
            "standardized_separation": 3.85,
        }
    )
    caption = format_pca_test_caption(pd.DataFrame([row]))
    assert "Centroid distance = 0.63" in caption
    assert "RMS radii (A, B) = 0.18, 0.14" in caption
    assert "Standardized separation = 3.85" in caption
