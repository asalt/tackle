import json
from pathlib import Path

import pandas as pd

from tackle import pca_replay
from tackle.pca_replay_rmd import render_pca_replay_rmd


def test_write_pca2_replay_preserves_exact_pre_svd_orientation(tmp_path, monkeypatch):
    captured = {}

    def fake_write_gctx(matrix, path, **kwargs):
        captured.update(
            mat=matrix.copy(),
            cdesc=kwargs["col_metadata"].copy(),
            rdesc=kwargs["row_metadata"].copy(),
            matrix_dtype=kwargs["matrix_dtype"],
            content_addressed=kwargs["content_addressed"],
        )
        output = Path(path)
        output.write_text("fake gctx\n", encoding="utf-8")
        return output

    monkeypatch.setattr(pca_replay, "write_gctx", fake_write_gctx)
    monkeypatch.setattr(pca_replay, "read_gctx_content_hash", lambda _path: "abc123")
    matrix = pd.DataFrame(
        [[-1.0, 2.0, 0.5], [1.0, -2.0, -0.5]],
        index=["SampleA", "SampleB"],
        columns=["101", "202", "303"],
    )
    metadata = pd.DataFrame(
        {"group": ["control", "case"], "recno": [1, 2]},
        index=["SampleA", "SampleB"],
    )

    files = pca_replay.write_pca2_replay(
        pca_matrix=matrix,
        sample_metadata=metadata,
        feature_symbols={"101": "AAA", "202": "BBB"},
        pca2_outname=str(tmp_path / "pca" / "analysis_pca2"),
        analysis_outpath=str(tmp_path),
        preprocessing={"center": True, "scale": False, "fillna": "min"},
        plot_parameters={
            "color": "group",
            "marker": None,
            "annotate": False,
            "frame": False,
            "encircle": False,
            "show_loadings": False,
            "ntop_loadings": 10,
            "file_formats": [".pdf"],
        },
        data_parameters={"taxon": "all"},
        metadata_colors={"group": {"control": "#111111", "case": "#eeeeee"}},
        separation_testing={
            "enabled": True,
            "group_fields": ["group"],
            "resolved_scopes": [
                {"name": "PC1_PC2", "pcs": ["PC1", "PC2"], "plot_key": "pc1_vs_pc2"}
            ],
            "p_adjust_method": "holm",
        },
    )

    assert list(captured["mat"].index) == ["SampleA", "SampleB"]
    assert list(captured["mat"].columns) == ["101", "202", "303"]
    assert list(captured["rdesc"].index) == ["SampleA", "SampleB"]
    assert list(captured["cdesc"].index) == ["101", "202", "303"]
    assert captured["cdesc"].loc["101", "GeneSymbol"] == "AAA"
    assert captured["matrix_dtype"] == "float64"
    assert captured["content_addressed"] is True

    context = json.loads(files.context_path.read_text(encoding="utf-8"))
    assert context["replay_contract_version"] == 5
    assert context["authoritative_input"]["orientation"] == "rows are samples; columns are features"
    assert context["prcomp_arguments"] == {"center": False, "scale.": False}
    assert context["sample_ids"] == ["SampleA", "SampleB"]
    assert context["feature_ids"] == ["101", "202", "303"]
    assert context["replot_default_pc_pairs"] == [[1, 2], [1, 3], [2, 3]]
    assert context["authoritative_input"]["storage_format"] == "gctx"
    assert context["authoritative_input"]["matrix_dtype"] == "float64"
    assert context["authoritative_input"]["gctx_content_hash"] == "abc123"
    assert context["authoritative_input"]["gctx_content_hash_algorithm"] in {
        "blake3",
        "blake2b-256",
    }
    assert context["separation_testing"]["group_fields"] == ["group"]
    assert files.stats_r_path.exists()
    assert "pca_welch_james" in files.stats_r_path.read_text(encoding="utf-8")
    assert files.pointer_path == tmp_path / "context" / "last_pca2_replay.json"


def test_pca_replay_rmd_uses_stored_matrix_without_preprocessing():
    rmd = render_pca_replay_rmd(title="PCA replay")

    assert "pca_mat <- as.matrix(stored_ds@mat)" in rmd
    assert "prcomp(pca_mat, center = FALSE, scale. = FALSE)" in rmd
    assert "pc_pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))" in rmd
    assert 'source("pca_stats.R", local = TRUE)' in rmd
    assert "pca_analyze_separation" in rmd
    assert "pca_analyze_single_pc_separation" in rmd
    assert "pca_plot_pairwise_separation" in rmd
    assert "including the sole comparison for a two-level factor" in rmd
    assert 'file.path(out_dir, "pca_single_pc_pairwise.tsv")' in rmd
    assert "without variance pooling, moderation, or" in rmd
    assert 'plot.caption.position = "plot"' in rmd
    assert "fillna_func" not in rmd
    assert "safe_scale" not in rmd
