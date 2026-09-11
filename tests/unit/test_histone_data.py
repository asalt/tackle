from types import SimpleNamespace

from click.testing import CliRunner
import pandas as pd

from tackle.containers import HistoneInfo
from tackle.histone_dedupe import apply_histone_collapse
from tackle.main import main


def test_requested_mouse_h3_genes_are_in_histone_group():
    expected = {
        "97908": "Hist1h3g",
        "260423": "Hist1h3f",
        "319148": "Hist1h3c",
        "319149": "Hist1h3d",
        "319150": "Hist1h3b",
        "319151": "Hist1h3e",
        "319152": "Hist1h3h",
        "319153": "Hist1h3i",
        "319154": "Hist2h3b",
    }
    observed = {}
    for info in HistoneInfo():
        observed.update(dict(zip(info.df["GeneID"], info.df["GeneSymbol"])))

    assert {geneid: observed.get(geneid) for geneid in expected} == expected


def _synthetic_histone_data():
    return pd.DataFrame(
        {
            "GeneID": ["10", "10", "10", "20", "20", "20", "30", "30", "30"],
            "Metric": ["TaxonID", "GeneSymbol", "AreaSum_dstrAdj"] * 3,
            "sample_a": ["10090", "H3a", 5.0, "10090", "H3b", 5.0, "10090", "X", 2.0],
            "sample_b": ["10090", "H3a", 7.0, "10090", "H3b", 7.0, "10090", "X", 3.0],
        }
    )


def test_histone_collapse_reports_only_the_removed_exact_duplicate():
    data = _synthetic_histone_data()
    histones = [
        SimpleNamespace(
            df=pd.DataFrame(
                {"GeneID": ["10", "20"], "GeneSymbol": ["H3a", "H3b"]}
            )
        )
    ]

    collapsed, excluded = apply_histone_collapse(data, histones)

    assert excluded == {"20"}
    assert set(collapsed["GeneID"]) == {"10", "30"}


def test_no_histone_collapse_retains_all_geneids():
    data = _synthetic_histone_data()

    uncollapsed, excluded = apply_histone_collapse(data, [], enabled=False)

    pd.testing.assert_frame_equal(uncollapsed, data)
    assert excluded == set()


def test_global_help_exposes_histone_collapse_override():
    result = CliRunner().invoke(main, ["--help"])

    assert result.exit_code == 0, result.output
    assert "--histone-collapse / --no-histone-collapse" in result.output
    assert "[default: histone-collapse]" in result.output
