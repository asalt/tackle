from pathlib import Path

from tackle.scriptgen import render_cluster2_sweep_skeleton


def test_make_cluster_skeleton_includes_cluster2_sweep(tmp_path: Path):
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

    script = render_cluster2_sweep_skeleton(conf_path=str(conf), k_start=3, k_end=5)
    assert f'CONF="{conf}"' in script
    assert "export TACKLE_CACHE=1" in script
    assert "export TACKLE_CLUSTER_DB=1" in script
    assert "run_cluster2_sweep" in script
    assert "cluster2 \\" in script
    assert "--cluster-db \\" in script
    assert "K_START=3" in script
    assert "K_END=5" in script
    assert "tackle cluster-summary" in script
