from pathlib import Path

from tackle.scriptgen import render_tacklerun_skeleton, summarize_config


def test_make_run_skeleton_includes_config_summary(tmp_path: Path):
    conf = tmp_path / "example.conf"
    conf.write_text(
        "\n".join(
            [
                "[s1]",
                "recno=1",
                "runno=1",
                "searchno=1",
                "label=LF",
                "group=A",
                "",
                "[s2]",
                "recno=2",
                "runno=1",
                "searchno=1",
                "label=LF",
                "group=B",
                "",
            ]
        )
    )

    summary = summarize_config(str(conf))
    assert summary.analysis_name == "example"
    assert summary.total_samples == 2
    assert "group" in summary.metadata_columns

    script = render_tacklerun_skeleton(conf_path=str(conf))
    assert f'CONF="{conf}"' in script
    assert "# Metadata columns (unique values):" in script
    assert "# - group: 2 (A, B)" in script
    assert "tackle \"${HEADMAIN[@]}\" \"$CONF\" \\" in script
    assert "VOLCANO_CONTRASTS=$(cat <<'EOF'" in script
    assert 'VOLCANO_LABEL_SCALE=1.0' in script
    assert 'VOLCANO_COMPARISON_LABEL_SCALE=1.0' in script
    assert 'VOLCANO_COMPARISON_WRAP_WIDTH=""' in script
    assert "read -r -d '' VOLCANO_CONTRASTS" not in script
    assert "plot_topdiff()" in script
    assert "run_pca() {" in script
    assert "run_cluster() {" in script
    assert 'pca2 --annotate --color "$DESIGN_COL" --max-pc 4' in script
    assert 'cluster2 --cut-by "$DESIGN_COL"' in script
    assert "EXPORT_NON_ZEROS=1" in script
    assert "HEAD_EXPORT_MSPC=(" in script
    assert "export --level area --level gct" in script
    assert "export --level MSPC --level evidence" in script
    assert "# run_pca && run_cluster" in script
    assert "mapfile -d '' -t files" in script
    assert "find \"$qstr\" -type f -name '*.tsv' -print0" in script
    assert "TOPDIFF_VARIANTS=(" in script
    assert "read -r -a v_arr" in script
    assert "local -n v_arr" not in script
