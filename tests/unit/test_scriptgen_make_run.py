import subprocess
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
    assert 'DESIGN_COL=""' in script
    assert script.index("# Metadata columns (unique values):") < script.index(
        'DESIGN_COL=""'
    )
    assert "tackle \"${HEADMAIN[@]}\" \"$CONF\" \\" in script
    assert "VOLCANO_CONTRASTS=$(cat <<'EOF'" in script
    assert 'VOLCANO_LABEL_SCALE=1.0' in script
    assert 'VOLCANO_COMPARISON_LABEL_SCALE=1.0' in script
    assert 'VOLCANO_COMPARISON_WRAP_WIDTH=""' in script
    assert "read -r -d '' VOLCANO_CONTRASTS" not in script
    assert "plot_topdiff()" in script
    assert "run_pca() {" in script
    assert "run_cluster() {" in script
    assert "run_correlation() {" in script
    assert (
        'pca2 --color "$DESIGN_COL" --max-pc 3 --annotate --encircle '
        '--test-by "$DESIGN_COL"'
    ) in script
    assert 'cluster2 --cut-by "$DESIGN_COL"' in script
    assert 'correlation --cluster --metric euclidean --annotate' in script
    assert 'correlation --cluster --metric pearson --annotate' in script
    assert 'correlation --cluster --metric pearson --z-score --annotate' in script
    assert (
        'correlation --cut-by "$DESIGN_COL" --cluster --metric euclidean --annotate'
    ) in script
    assert (
        'correlation --cut-by "$DESIGN_COL" --cluster --metric pearson --annotate'
    ) in script
    assert (
        'correlation --cut-by "$DESIGN_COL" --cluster --metric pearson '
        '--z-score --annotate'
    ) in script
    assert "EXPORT_NON_ZEROS=1" in script
    assert "HEAD_EXPORT_MSPC=(" in script
    assert "export --level area --level gct" in script
    assert "export --level MSPC --level evidence" in script
    assert "# run_correlation" in script
    assert "# run_pca && run_cluster" in script
    assert "mapfile -d '' -t files" in script
    assert "find \"$qstr\" -type f -name '*.tsv' -print0" in script
    assert "TOPDIFF_VARIANTS=(" in script
    assert "read -r -a v_arr" in script
    assert "local -n v_arr" not in script
    assert 'thebasename="$(basename "${CONF%.conf}")"' in script

    syntax_check = subprocess.run(
        ["bash", "-n"],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )
    assert syntax_check.returncode == 0, syntax_check.stderr


def test_make_run_design_uses_highest_nonunique_metadata_cardinality(
    tmp_path: Path,
) -> None:
    conf = tmp_path / "ranked.conf"
    sections = []
    for index in range(1, 7):
        sections.extend(
            [
                f"[s{index}]",
                f"recno={index}",
                "runno=1",
                "searchno=1",
                "label=LF",
                f"tube=tube{index}",
                f"batch=B{((index - 1) % 3) + 1}",
                f"group={'A' if index <= 3 else 'B'}",
                "",
            ]
        )
    conf.write_text("\n".join(sections))

    summary = summarize_config(str(conf))
    script = render_tacklerun_skeleton(conf_path=str(conf))

    assert summary.recommended_design_columns == ["batch", "group"]
    assert 'DESIGN_COL="batch"' in script
    assert script.index("# Metadata columns (unique values):") < script.index(
        'DESIGN_COL="batch"'
    )
