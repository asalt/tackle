from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from tackle.limma_replay import (
    resolve_replay_dir,
    validate_replay_dir,
    write_limma_replay_files,
)
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


def test_write_limma_replay_files_prefers_explicit_expression_matrix(tmp_path: Path, monkeypatch) -> None:
    captured = {}

    def fake_write_gct(*, out_path, mat, cdesc, rdesc, precision=4):
        captured["mat"] = mat.copy()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("# dummy gct\n", encoding="utf-8")
        return out_path

    monkeypatch.setattr("tackle.limma_replay._write_gct", fake_write_gct)

    analysis = tmp_path / "analysis"
    volcano = analysis / "volcano" / "mouse"
    sample_metadata = pd.DataFrame({"group": ["A", "B"]}, index=["S1", "S2"])

    explicit_edata = pd.DataFrame(
        {"S1": [10.5, 20.5], "S2": [30.5, 40.5]},
        index=["101", "202"],
    )
    result_df = pd.DataFrame(
        {"S1": [1.0, 2.0], "S2": [3.0, 4.0], "pAdj": [0.01, 0.02]},
        index=["101", "202"],
    )

    files = write_limma_replay_files(
        analysis_dir=str(analysis),
        volcano_dir=str(volcano),
        results={"A-B=A-B": result_df},
        sample_metadata=sample_metadata,
        expression_matrix=explicit_edata,
        impute_missing_values=True,
        imputation_backend="gaussian",
        gaussian_method="legacy",
        force=True,
    )

    assert captured["mat"].equals(explicit_edata)

    context = json.loads(files.context_path.read_text(encoding="utf-8"))
    assert context["impute_missing_values"] is True
    assert context["imputation_backend"] == "gaussian"
    assert context["gaussian_method"] == "legacy"
    assert context["path_base"] == "replay_dir"
    assert context["analysis_dir"] == "../../.."
    assert context["volcano_dir"] == "../.."
    assert context["replay_dir"] == "."
    assert context["gct_path"] == "limma_input.gct"
    assert context["stored_matrix_is_authoritative"] is True
    assert context["stored_matrix_role"] == "limma_input_imputed"
    assert context["has_pre_impute_matrix"] is False
    assert (files.replay_dir / "replay_explore.Rmd").exists()
    assert (files.replay_dir / "render_replay_explore.sh").exists()
    rmd = (files.replay_dir / "replay_explore.Rmd").read_text(encoding="utf-8")
    assert "limma_input.gct" in rmd
    assert "Stored Matrix Replay" in rmd
    assert "Optional Gaussian Recompute" in rmd


def test_write_limma_replay_files_writes_pre_impute_matrix_for_recompute(tmp_path: Path, monkeypatch) -> None:
    captured = {}

    def fake_write_gct(*, out_path, mat, cdesc, rdesc, precision=4):
        captured[out_path.name] = mat.copy()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("# dummy gct\n", encoding="utf-8")
        return out_path

    monkeypatch.setattr("tackle.limma_replay._write_gct", fake_write_gct)

    analysis = tmp_path / "analysis"
    volcano = analysis / "volcano" / "mouse"
    sample_metadata = pd.DataFrame({"group": ["A", "B"]}, index=["S1", "S2"])

    imputed_edata = pd.DataFrame(
        {"S1": [10.5, 20.5], "S2": [30.5, 40.5]},
        index=["101", "202"],
    )
    pre_impute_edata = pd.DataFrame(
        {"S1": [10.5, float("nan")], "S2": [float("nan"), 40.5]},
        index=["101", "202"],
    )
    result_df = pd.DataFrame(
        {"S1": [1.0, 2.0], "S2": [3.0, 4.0], "pAdj": [0.01, 0.02]},
        index=["101", "202"],
    )

    files = write_limma_replay_files(
        analysis_dir=str(analysis),
        volcano_dir=str(volcano),
        results={"A-B=A-B": result_df},
        sample_metadata=sample_metadata,
        expression_matrix=imputed_edata,
        pre_impute_expression_matrix=pre_impute_edata,
        impute_missing_values=True,
        imputation_backend="gaussian",
        gaussian_method="mqish",
        force=True,
    )

    assert captured["limma_input.gct"].equals(imputed_edata)
    assert captured["limma_input_pre_impute.gct"].equals(pre_impute_edata)
    assert files.pre_impute_gct_path == files.replay_dir / "limma_input_pre_impute.gct"

    context = json.loads(files.context_path.read_text(encoding="utf-8"))
    assert context["has_pre_impute_matrix"] is True
    assert context["gct_path"] == "limma_input.gct"
    assert context["pre_impute_gct_path"] == "limma_input_pre_impute.gct"
    assert context["pre_impute_matrix_shape"] == [2, 2]
    assert context["recompute_imputation_supported"] is True

    pointer = json.loads(files.pointer_path.read_text(encoding="utf-8"))
    assert pointer["path_base"] == "analysis_dir"
    assert pointer["analysis_dir"] == "."
    assert pointer["replay_dir"].startswith("volcano/mouse/replay/")
    assert pointer["gct_path"].endswith("/limma_input.gct")
    assert pointer["pre_impute_gct_path"].endswith("/limma_input_pre_impute.gct")
    assert not Path(pointer["replay_dir"]).is_absolute()
    assert not Path(pointer["gct_path"]).is_absolute()
    assert not Path(pointer["pre_impute_gct_path"]).is_absolute()

    rmd = (files.replay_dir / "replay_explore.Rmd").read_text(encoding="utf-8")
    assert 'pre_impute_gct_path <- "limma_input_pre_impute.gct"' in rmd
    assert "mat_stored" in rmd
