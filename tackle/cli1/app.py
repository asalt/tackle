from __future__ import annotations

import os
import shutil
from pathlib import Path

import click

from .analysis_dir import AnalysisDir, BuildSpec, derive_analysis_dir
from .hashutil import gct_dims, sha256_file, tsv_dims
from .manifest import ManifestDB, write_manifest_json


def _resolve_analysis_dir(
    *,
    analysis_dir: str | None,
    conf: str | None,
    name: str | None,
    result_dir: str,
) -> Path:
    if analysis_dir:
        return Path(analysis_dir)
    if not conf or not name:
        raise click.UsageError("Provide either --analysis-dir OR (--conf AND --name).")
    return derive_analysis_dir(result_dir=result_dir, conf_path=conf, name=name)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--analysis-dir", type=click.Path(), default=None, show_default=True)
@click.option("--conf", type=click.Path(exists=True, dir_okay=False), default=None, show_default=True)
@click.option("--name", type=str, default=None, show_default=True)
@click.option("--result-dir", type=click.Path(), default="./results", show_default=True)
@click.option("--cache/--no-cache", default=False, show_default=True, help="Enable TACKLE_CACHE for this run.")
@click.pass_context
def cli(ctx, analysis_dir, conf, name, result_dir, cache):
    """tackle1: clean-slate CLI (experimental, not wired into `tackle` entrypoint)."""
    ctx.ensure_object(dict)
    ctx.obj["analysis_dir"] = analysis_dir
    ctx.obj["conf"] = conf
    ctx.obj["name"] = name
    ctx.obj["result_dir"] = result_dir
    ctx.obj["cache"] = bool(cache)


@cli.command("init")
@click.option("--data-dir", type=click.Path(), default="./data/e2g/", show_default=True)
@click.option("--taxon", type=str, default="all", show_default=True)
@click.option("--non-zeros", type=int, default=1, show_default=True)
@click.option("--unique-pepts", type=int, default=0, show_default=True)
@click.option("--normalize-across-species/--no-normalize-across-species", default=False, show_default=True)
@click.option("--fill-na-zero/--no-fill-na-zero", default=True, show_default=True)
@click.option("--impute-missing-values/--no-impute-missing-values", default=False, show_default=True)
@click.option("--imputation-backend", type=click.Choice(["gaussian", "lupine"]), default="gaussian", show_default=True)
@click.option("--lupine-mode", type=click.Choice(["local", "joint"]), default="local", show_default=True)
@click.option(
    "--norm-method",
    type=click.Choice(["none", "median", "ifot", "ifot_ki", "ifot_tf", "quantile75", "quantile90", "genefile"]),
    default="none",
    show_default=True,
)
@click.option("--genefile-norm", type=click.Path(exists=True, dir_okay=False), default=None, show_default=True)
@click.option(
    "--ref-group",
    "ref_group_cols",
    multiple=True,
    default=(),
    show_default=True,
    help="One or more metadata columns to match samples when normalizing by a control sample (repeatable).",
)
@click.option(
    "--ref-label-col",
    default="label",
    show_default=True,
    help="Metadata column used to identify the control sample within each ref-group (default: label).",
)
@click.option(
    "--ref-control-value",
    default=None,
    show_default=True,
    help="Value in --ref-label-col used as the denominator/control within each ref-group.",
)
@click.option("--force", is_flag=True, default=False, show_default=True, help="Overwrite existing analysis.json if present.")
@click.pass_context
def init_cmd(
    ctx,
    data_dir,
    taxon,
    non_zeros,
    unique_pepts,
    normalize_across_species,
    fill_na_zero,
    impute_missing_values,
    imputation_backend,
    lupine_mode,
    norm_method,
    genefile_norm,
    ref_group_cols,
    ref_label_col,
    ref_control_value,
    force,
):
    """Initialize an analysis directory (write analysis.json + manifest)."""
    analysis_dir = _resolve_analysis_dir(
        analysis_dir=ctx.obj.get("analysis_dir"),
        conf=ctx.obj.get("conf"),
        name=ctx.obj.get("name"),
        result_dir=ctx.obj.get("result_dir") or "./results",
    )
    conf = ctx.obj.get("conf")
    name = ctx.obj.get("name")
    if conf is None or name is None:
        raise click.UsageError("init requires --conf and --name (or use derived --analysis-dir with those set).")

    if norm_method == "genefile" and not genefile_norm:
        raise click.BadParameter("--genefile-norm is required when --norm-method genefile is used.")
    if norm_method != "genefile":
        genefile_norm = None

    ref_group_cols = [str(x) for x in ref_group_cols if str(x).strip()] or None
    if ref_control_value and not ref_group_cols:
        raise click.BadParameter("--ref-group is required when --ref-control-value is set.")
    if ref_group_cols and not ref_control_value:
        raise click.BadParameter("--ref-control-value is required when --ref-group is set.")

    build_spec = BuildSpec(
        taxon=taxon,
        non_zeros=int(non_zeros),
        unique_pepts=int(unique_pepts),
        normalize_across_species=bool(normalize_across_species),
        fill_na_zero=bool(fill_na_zero),
        impute_missing_values=bool(impute_missing_values),
        imputation_backend=str(imputation_backend),
        lupine_mode=str(lupine_mode),
        norm_method=str(norm_method),
        genefile_norm=str(genefile_norm) if genefile_norm else None,
        ref_group_cols=ref_group_cols,
        ref_label_col=str(ref_label_col or "label"),
        ref_control_value=str(ref_control_value) if ref_control_value else None,
    )

    ad = AnalysisDir(analysis_dir)
    cfg = ad.write_analysis_config(
        conf_path=conf,
        name=name,
        result_dir=ctx.obj.get("result_dir") or "./results",
        data_dir=data_dir,
        build_spec=build_spec,
        force=bool(force),
    )

    ManifestDB(ad.manifest_sqlite_path).init_db()
    write_manifest_json(
        path=ad.manifest_json_path,
        analysis=cfg.to_json(),
        artifacts={},
    )

    click.echo(f"Initialized analysis dir: {ad.path}")
    click.echo(f"- analysis.json: {ad.analysis_json_path}")
    click.echo(f"- manifest.sqlite: {ad.manifest_sqlite_path}")
    click.echo(f"- inputs: {ad.inputs_dir}")
    click.echo(f"- matrices: {ad.matrices_dir}")


def _data_from_analysis(cfg, *, only_local: bool = False):
    from tackle.containers import Data

    norm = cfg.build_spec.norm_method
    flags = {
        "ifot": norm == "ifot",
        "ifot_ki": norm == "ifot_ki",
        "ifot_tf": norm == "ifot_tf",
        "median": norm == "median",
        "quantile75": norm == "quantile75",
        "quantile90": norm == "quantile90",
        "genefile_norm": cfg.build_spec.genefile_norm,
    }
    analysis_dir = Path(cfg.analysis_dir)
    legacy_out = analysis_dir / "_tackle"

    norm_info = None
    if cfg.build_spec.ref_group_cols and cfg.build_spec.ref_control_value:
        norm_info = {
            "control": cfg.build_spec.ref_control_value,
            "group": cfg.build_spec.ref_group_cols,
            "label": cfg.build_spec.ref_label_col or "label",
        }
    data_obj = Data(
        experiment_file=str(analysis_dir / cfg.conf_copied),
        base_dir=cfg.result_dir,  # for shared cache placement
        data_dir=cfg.data_dir,
        taxon=cfg.build_spec.taxon,
        non_zeros=cfg.build_spec.non_zeros,
        unique_pepts=cfg.build_spec.unique_pepts,
        normalize_across_species=cfg.build_spec.normalize_across_species,
        fill_na_zero=cfg.build_spec.fill_na_zero,
        impute_missing_values=cfg.build_spec.impute_missing_values,
        imputation_backend=cfg.build_spec.imputation_backend,
        lupine_mode=cfg.build_spec.lupine_mode,
        only_local=only_local,
        metrics=False,
        norm_info=norm_info,
        set_outpath=False,
        outpath=str(legacy_out),
        outpath_name=analysis_dir.name,
        **flags,
    )
    return data_obj


@cli.command("build")
@click.option("--only-local/--no-only-local", default=False, show_default=True)
@click.option("--force/--no-force", default=False, show_default=True, help="Recompute and overwrite canonical matrices.")
@click.pass_context
def build_cmd(ctx, only_local, force):
    """Ensure the base matrices exist (area TSV + GCT + MSPC TSV)."""
    analysis_dir = _resolve_analysis_dir(
        analysis_dir=ctx.obj.get("analysis_dir"),
        conf=ctx.obj.get("conf"),
        name=ctx.obj.get("name"),
        result_dir=ctx.obj.get("result_dir") or "./results",
    )
    ad = AnalysisDir(analysis_dir)
    cfg = ad.load_analysis_config()

    if ctx.obj.get("cache"):
        os.environ["TACKLE_CACHE"] = "1"

    manifest = ManifestDB(ad.manifest_sqlite_path)
    run_id = manifest.start_run(command="build", params={"force": bool(force), "only_local": bool(only_local)})

    canonical = ad.canonical_matrix_paths()
    ad.ensure_layout()

    try:
        data_obj = None
        for role, out_path in canonical.items():
            if out_path.exists() and not force:
                continue
            if data_obj is None:
                data_obj = _data_from_analysis(cfg, only_local=bool(only_local))

            if role == "area_tsv":
                src = data_obj.perform_data_export("area", force=bool(force))
            elif role == "gct":
                src = data_obj.perform_data_export("gct", force=bool(force))
            elif role == "mspc_tsv":
                src = data_obj.perform_data_export("MSPC", force=bool(force))
            else:
                raise ValueError(f"Unknown role: {role}")

            out_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, out_path)

        artifacts_summary = {}
        for role, out_path in canonical.items():
            if not out_path.exists():
                continue
            sha = sha256_file(out_path)
            if role == "gct":
                rows, cols = gct_dims(out_path)
            else:
                rows, cols = tsv_dims(out_path)
            rel = ad.relpath(out_path)
            manifest.upsert_artifact(
                run_id=run_id,
                type="matrix",
                role=role,
                path=rel,
                sha256=sha,
                rows=rows,
                cols=cols,
                params={"canonical": True},
            )
            artifacts_summary[role] = {"path": rel, "sha256": sha, "rows": rows, "cols": cols}

        write_manifest_json(
            path=ad.manifest_json_path,
            analysis=cfg.to_json(),
            artifacts={"matrices": artifacts_summary},
        )

        manifest.finish_run(run_id=run_id, status="ok")
    except Exception as e:
        manifest.finish_run(run_id=run_id, status="error", message=str(e))
        raise

    click.echo(f"Built base matrices for: {ad.path}")
    for role, out_path in canonical.items():
        click.echo(f"- {role}: {out_path}")


@cli.command("ls")
@click.pass_context
def ls_cmd(ctx):
    """Show analysis metadata and known artifact paths."""
    analysis_dir = _resolve_analysis_dir(
        analysis_dir=ctx.obj.get("analysis_dir"),
        conf=ctx.obj.get("conf"),
        name=ctx.obj.get("name"),
        result_dir=ctx.obj.get("result_dir") or "./results",
    )
    ad = AnalysisDir(analysis_dir)
    cfg = ad.load_analysis_config()
    manifest = ManifestDB(ad.manifest_sqlite_path)

    click.echo(f"Analysis dir: {ad.path}")
    click.echo(f"- Config: {ad.resolve_copied_conf_path(cfg)}")
    click.echo(f"- Build spec: taxon={cfg.build_spec.taxon} non_zeros={cfg.build_spec.non_zeros} norm={cfg.build_spec.norm_method}")
    click.echo("Canonical matrices:")
    for role, out_path in ad.canonical_matrix_paths().items():
        status = "OK" if out_path.exists() else "MISSING"
        click.echo(f"- {role}: {out_path} [{status}]")

    click.echo("Manifest artifacts (most recent first):")
    for art in manifest.list_artifacts(type="matrix")[:20]:
        click.echo(f"- {art.role}: {art.path} ({art.rows}x{art.cols})")
