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
