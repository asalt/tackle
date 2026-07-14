from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from tackle.limma_replay import (
    _find_existing_gct,
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


def test_modeled_gct_selection_never_uses_pre_impute_audit_file(tmp_path: Path) -> None:
    replay = tmp_path / "replay" / "run1"
    replay.mkdir(parents=True)
    modeled = replay / "limma_input_n2x3.gct"
    pre_impute = replay / "limma_input_pre_impute_n2x3.gct"
    modeled.write_text("# modeled\n", encoding="utf-8")
    pre_impute.write_text("# audit only\n", encoding="utf-8")
    context = replay / "limma_replay_context.json"
    context.write_text(
        json.dumps({"gct_path": modeled.name}, indent=2) + "\n",
        encoding="utf-8",
    )

    selected, _ = validate_replay_dir(replay)
    assert selected == modeled.resolve()
    assert _find_existing_gct(replay / "limma_input.gct") == modeled

    context.write_text(
        json.dumps({"gct_path": pre_impute.name}, indent=2) + "\n",
        encoding="utf-8",
    )
    selected_bad_context, _ = validate_replay_dir(replay)
    assert selected_bad_context == modeled

    context.write_text("{}\n", encoding="utf-8")
    selected_fallback, _ = validate_replay_dir(replay)
    assert selected_fallback == modeled


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
    rmd = bundle.rmd_path.read_text(encoding="utf-8")
    assert "lmFit(" in rmd
    assert "contrasts.fit(" in rmd
    assert "eBayes(" in rmd
    assert "deterministic_gaussian_impute" not in rmd
    assert "prcomp(" not in rmd
    assert "sample-distributions" not in rmd
    assert "sample-correlation" not in rmd
    assert "boxplot(" not in rmd


def test_write_limma_replay_files_prefers_explicit_expression_matrix(tmp_path: Path, monkeypatch) -> None:
    captured = {}

    def fake_write_gct(*, out_path, mat, cdesc, rdesc, precision=4):
        captured["mat"] = mat.copy()
        captured["precision"] = precision
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("# dummy gct\n", encoding="utf-8")
        return out_path

    monkeypatch.setattr("tackle.limma_replay._write_gct", fake_write_gct)
    monkeypatch.setattr(
        "tackle.limma_replay.read_gctx_content_hash", lambda _path: "abc123"
    )

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
        expression_matrix_source="stat_model_modeled_matrix",
        impute_missing_values=True,
        imputation_backend="gaussian",
        gaussian_method="legacy",
        force=True,
    )

    assert captured["mat"].equals(explicit_edata)
    assert captured["precision"] == 17

    context = json.loads(files.context_path.read_text(encoding="utf-8"))
    assert context["impute_missing_values"] is True
    assert context["imputation_backend"] == "gaussian"
    assert context["gaussian_method"] == "legacy"
    assert context["path_base"] == "replay_dir"
    assert context["analysis_dir"] == "../../../.."
    assert context["volcano_dir"] == "../.."
    assert context["replay_dir"] == "."
    assert context["gct_path"] == "limma_input.gctx"
    assert context["stored_matrix_is_authoritative"] is True
    assert context["stored_matrix_role"] == "limma_input_imputed"
    assert context["expression_matrix_source"] == "stat_model_modeled_matrix"
    assert context["replay_contract_version"] == 4
    assert context["gctx_content_hash"] == "abc123"
    assert context["gctx_content_hash_algorithm"] in {"blake3", "blake2b-256"}
    assert context["replay_scope"] == "limma_only"
    assert context["matrix_storage_format"] == "gctx"
    assert context["matrix_storage_dtype"] == "float64"
    assert context["has_pre_impute_matrix"] is False
    assert context["recompute_imputation_supported"] is False
    assert context["pca_included"] is False
    assert (files.replay_dir / "replay_explore.Rmd").exists()
    assert (files.replay_dir / "render_replay_explore.sh").exists()
    rmd = (files.replay_dir / "replay_explore.Rmd").read_text(encoding="utf-8")
    assert "limma_input.gctx" in rmd
    assert "Load the Modeled Matrix" in rmd
    assert "lmFit(" in rmd
    assert "as.matrix(mat_stored)" in rmd
    assert "contrasts.fit(" in rmd
    assert "eBayes(" in rmd
    assert "deterministic_gaussian_impute" not in rmd
    assert "Optional Gaussian Recompute" not in rmd
    assert "prcomp(" not in rmd
    assert "sample-distributions" not in rmd
    assert "sample-correlation" not in rmd


def test_write_limma_replay_files_keeps_pre_impute_matrix_as_audit_only(
    tmp_path: Path, monkeypatch
) -> None:
    captured = {}

    def fake_write_gct(*, out_path, mat, cdesc, rdesc, precision=4):
        captured[out_path.name] = mat.copy()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("# dummy gct\n", encoding="utf-8")
        return out_path

    monkeypatch.setattr("tackle.limma_replay._write_gct", fake_write_gct)
    monkeypatch.setattr(
        "tackle.limma_replay.read_gctx_content_hash", lambda _path: "abc123"
    )

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

    assert captured["limma_input.gctx"].equals(imputed_edata)
    assert captured["limma_input_pre_impute.gctx"].equals(pre_impute_edata)
    assert files.pre_impute_gct_path == files.replay_dir / "limma_input_pre_impute.gctx"

    context = json.loads(files.context_path.read_text(encoding="utf-8"))
    assert context["has_pre_impute_matrix"] is True
    assert context["gct_path"] == "limma_input.gctx"
    assert context["pre_impute_gct_path"] == "limma_input_pre_impute.gctx"
    assert context["pre_impute_matrix_shape"] == [2, 2]
    assert context["pre_impute_matrix_role"] == "audit_only"
    assert context["recompute_imputation_supported"] is False

    pointer = json.loads(files.pointer_path.read_text(encoding="utf-8"))
    assert pointer["path_base"] == "analysis_dir"
    assert pointer["analysis_dir"] == "."
    assert pointer["replay_dir"].startswith("volcano/mouse/replay/")
    assert pointer["gct_path"].endswith("/limma_input.gctx")
    assert pointer["pre_impute_gct_path"].endswith("/limma_input_pre_impute.gctx")
    assert not Path(pointer["replay_dir"]).is_absolute()
    assert not Path(pointer["gct_path"]).is_absolute()
    assert not Path(pointer["pre_impute_gct_path"]).is_absolute()

    rmd = (files.replay_dir / "replay_explore.Rmd").read_text(encoding="utf-8")
    assert "mat_stored" in rmd
    assert "mat_pre_impute" not in rmd
    assert "deterministic_gaussian_impute" not in rmd


def test_write_limma_replay_files_reuses_existing_gctx_for_same_inputs(tmp_path: Path) -> None:
    analysis = tmp_path / "analysis"
    volcano = analysis / "volcano" / "mouse"
    sample_metadata = pd.DataFrame({"group": ["A", "B"]}, index=["S1", "S2"])
    edata = pd.DataFrame({"S1": [10.5, 20.5], "S2": [30.5, 40.5]}, index=["101", "202"])
    result_df = pd.DataFrame(
        {"S1": [1.0, 2.0], "S2": [3.0, 4.0], "pAdj": [0.01, 0.02]},
        index=["101", "202"],
    )

    first = write_limma_replay_files(
        analysis_dir=str(analysis),
        volcano_dir=str(volcano),
        results={"A-B=A-B": result_df},
        sample_metadata=sample_metadata,
        expression_matrix=edata,
        force=True,
    )
    first_mtime_ns = first.gct_path.stat().st_mtime_ns
    second = write_limma_replay_files(
        analysis_dir=str(analysis),
        volcano_dir=str(volcano),
        results={"A-B=A-B": result_df},
        sample_metadata=sample_metadata,
        expression_matrix=edata,
        force=True,
    )

    assert first.replay_dir == second.replay_dir
    assert first.gct_path == second.gct_path
    assert second.gct_path.stat().st_mtime_ns == first_mtime_ns
    assert first.gct_path.name.startswith("limma_input.")
    assert first.gct_path.suffix == ".gctx"
