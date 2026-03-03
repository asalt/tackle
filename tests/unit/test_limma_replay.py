from __future__ import annotations

import json
from pathlib import Path

import pytest

from tackle.limma_replay import resolve_replay_dir, validate_replay_dir
from tackle.rmd_replay import write_limma_replay_bundle


def _write_replay_dir(replay_dir: Path) -> None:
    replay_dir.mkdir(parents=True, exist_ok=True)
    (replay_dir / "limma_input.gct").write_text("# dummy gct\n", encoding="utf-8")
    (replay_dir / "limma_replay_context.json").write_text(
        json.dumps(
            {
                "run_id": replay_dir.name,
                "result_labels": [],
                "contrast_expressions": [],
                "direct_coefs": [],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def test_resolve_replay_dir_via_pointer(tmp_path: Path) -> None:
    analysis = tmp_path / "analysis"
    replay = analysis / "volcano" / "mouse" / "replay" / "abc123"
    _write_replay_dir(replay)

    pointer = analysis / "context" / "last_volcano_replay.json"
    pointer.parent.mkdir(parents=True, exist_ok=True)
    pointer.write_text(
        json.dumps({"replay_dir": str(replay)}, indent=2) + "\n", encoding="utf-8"
    )

    resolved = resolve_replay_dir(analysis_dir=analysis)
    assert resolved == replay

    gct, ctx = validate_replay_dir(resolved)
    assert gct.name.endswith(".gct")
    assert ctx.name.endswith(".json")


def test_resolve_replay_dir_single_discovery(tmp_path: Path) -> None:
    analysis = tmp_path / "analysis"
    replay = analysis / "volcano" / "mouse" / "replay" / "onlyone"
    _write_replay_dir(replay)

    resolved = resolve_replay_dir(analysis_dir=analysis)
    assert resolved == replay


def test_resolve_replay_dir_multiple_requires_selection(tmp_path: Path) -> None:
    analysis = tmp_path / "analysis"
    _write_replay_dir(analysis / "volcano" / "mouse" / "replay" / "r1")
    _write_replay_dir(analysis / "volcano" / "mouse" / "replay" / "r2")

    with pytest.raises(RuntimeError):
        resolve_replay_dir(analysis_dir=analysis)


def test_resolve_replay_dir_from_volcano_dir_and_run_id(tmp_path: Path) -> None:
    analysis = tmp_path / "analysis"
    volcano = analysis / "volcano" / "mouse"
    _write_replay_dir(volcano / "replay" / "r1")
    _write_replay_dir(volcano / "replay" / "r2")

    resolved = resolve_replay_dir(analysis_dir=analysis, volcano_dir=volcano, run_id="r2")
    assert resolved == volcano / "replay" / "r2"


def test_write_limma_replay_bundle(tmp_path: Path) -> None:
    replay = tmp_path / "replay" / "abc"
    _write_replay_dir(replay)

    out = tmp_path / "report" / "rmd"
    bundle = write_limma_replay_bundle(
        report_dir=str(out),
        replay_dir=str(replay),
        title="Test Report",
        copy_inputs=True,
        force=True,
        render=False,
    )

    assert bundle.rmd_path.exists()
    assert bundle.render_sh.exists()
    assert bundle.gct_path.exists()
    assert bundle.context_path.exists()

