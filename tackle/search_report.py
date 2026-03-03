#!/usr/bin/env python3
"""Build methods-focused MS search reports from FragPipe metadata files."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import html
import pathlib
import re
import sys
from dataclasses import dataclass

WORKSPACE_MOUNT_PREFIX = "/mnt/e/MSPC001494/"
DEFAULT_RUN_DIRS = (
    pathlib.Path("data/processed/prof/fragpipe"),
    pathlib.Path("data/processed/IMAC/fragpipe"),
)
DEFAULT_OUT_HTML = pathlib.Path("results/fragpipe_methods_report.html")
DEFAULT_OUT_MZML_TSV = pathlib.Path("results/fragpipe_mzml_tracking.tsv")


@dataclass
class FraggerParam:
    key: str
    value: str
    comment: str
    active: bool


@dataclass
class ManifestEntry:
    run_name: str
    source_path: str
    record: str
    group: str
    acquisition: str
    file_name: str
    workspace_candidate: str
    exists_source_path: bool
    exists_workspace_candidate: bool

    @property
    def status(self) -> str:
        return "present" if self.exists_source_path or self.exists_workspace_candidate else "missing"


@dataclass
class RunData:
    run_name: str
    run_dir: pathlib.Path
    workflow_versions: dict[str, str]
    workflow_params: dict[str, str]
    fragger_params: dict[str, str]
    fragger_rows: list[FraggerParam]
    manifest_entries: list[ManifestEntry]
    parse_warnings: list[str]


@dataclass
class SearchReportOutputs:
    out_html: pathlib.Path
    out_mzml_tsv: pathlib.Path
    manifest_total: int
    missing_total: int


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        action="append",
        dest="run_dirs",
        help=(
            "FragPipe output directory. "
            "Can be provided multiple times. "
            "Default: data/processed/prof/fragpipe and data/processed/IMAC/fragpipe"
        ),
    )
    parser.add_argument(
        "--out-html",
        default=str(DEFAULT_OUT_HTML),
        help="Output HTML file path.",
    )
    parser.add_argument(
        "--out-mzml-tsv",
        default=str(DEFAULT_OUT_MZML_TSV),
        help="Output TSV with manifest-based mzML tracking.",
    )
    parser.add_argument(
        "--workspace-root",
        default=".",
        help="Workspace root used for local path remapping checks.",
    )
    return parser.parse_args(argv)


def infer_run_name(run_dir: pathlib.Path) -> str:
    if run_dir.name.lower() == "fragpipe" and run_dir.parent != run_dir:
        return run_dir.parent.name
    return run_dir.name


def parse_workflow(path: pathlib.Path) -> tuple[dict[str, str], dict[str, str]]:
    versions: dict[str, str] = {}
    params: dict[str, str] = {}
    version_re = re.compile(r"^#\s*(.+?)\s+version\s+(.+)$")
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                match = version_re.match(line)
                if match:
                    versions[match.group(1).strip()] = match.group(2).strip()
                continue
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            params[key.strip()] = value.strip()
    return versions, params


def parse_assignment(line: str) -> tuple[str, str, str] | None:
    comment = ""
    core = line
    if "#" in line:
        core, comment = line.split("#", 1)
        comment = comment.strip()
    if "=" not in core:
        return None
    key, value = core.split("=", 1)
    key = key.strip()
    value = value.strip()
    if not key:
        return None
    return key, value, comment


def parse_fragger_params(path: pathlib.Path) -> tuple[dict[str, str], list[FraggerParam]]:
    active_params: dict[str, str] = {}
    rows: list[FraggerParam] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                parsed = parse_assignment(stripped[1:].lstrip())
                if parsed is None:
                    continue
                key, value, comment = parsed
                rows.append(FraggerParam(key=key, value=value, comment=comment, active=False))
                continue
            parsed = parse_assignment(stripped)
            if parsed is None:
                continue
            key, value, comment = parsed
            active_params[key] = value
            rows.append(FraggerParam(key=key, value=value, comment=comment, active=True))
    return active_params, rows


def remap_manifest_path(source_path: str, workspace_root: pathlib.Path) -> pathlib.Path | None:
    text = source_path.strip()
    if not text:
        return None

    as_path = pathlib.Path(text)
    if as_path.is_absolute() and text.startswith(str(workspace_root)):
        return as_path

    if text.startswith(WORKSPACE_MOUNT_PREFIX):
        rel = text[len(WORKSPACE_MOUNT_PREFIX) :]
        return workspace_root / rel

    if not as_path.is_absolute():
        return workspace_root / as_path

    return None


def parse_manifest(
    path: pathlib.Path,
    run_name: str,
    workspace_root: pathlib.Path,
) -> list[ManifestEntry]:
    entries: list[ManifestEntry] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                parts = line.split()
            source_path = parts[0].strip() if len(parts) >= 1 else ""
            record = parts[1].strip() if len(parts) >= 2 else ""
            group = parts[2].strip() if len(parts) >= 3 else ""
            acquisition = parts[3].strip() if len(parts) >= 4 else ""

            candidate_path = remap_manifest_path(source_path, workspace_root)
            candidate_text = str(candidate_path) if candidate_path else ""

            exists_source = pathlib.Path(source_path).exists() if source_path else False
            exists_candidate = candidate_path.exists() if candidate_path else False

            entries.append(
                ManifestEntry(
                    run_name=run_name,
                    source_path=source_path,
                    record=record,
                    group=group,
                    acquisition=acquisition,
                    file_name=pathlib.Path(source_path).name if source_path else "",
                    workspace_candidate=candidate_text,
                    exists_source_path=exists_source,
                    exists_workspace_candidate=exists_candidate,
                )
            )
    return entries


def load_run(run_dir: pathlib.Path, workspace_root: pathlib.Path) -> RunData:
    run_name = infer_run_name(run_dir)
    warnings: list[str] = []

    workflow_path = run_dir / "fragpipe.workflow"
    fragger_path = run_dir / "fragger.params"
    manifest_path = run_dir / "fragpipe-files.fp-manifest"

    workflow_versions: dict[str, str] = {}
    workflow_params: dict[str, str] = {}
    fragger_params: dict[str, str] = {}
    fragger_rows: list[FraggerParam] = []
    manifest_entries: list[ManifestEntry] = []

    if workflow_path.exists():
        workflow_versions, workflow_params = parse_workflow(workflow_path)
    else:
        warnings.append(f"Missing file: {workflow_path}")

    if fragger_path.exists():
        fragger_params, fragger_rows = parse_fragger_params(fragger_path)
    else:
        warnings.append(f"Missing file: {fragger_path}")

    if manifest_path.exists():
        manifest_entries = parse_manifest(manifest_path, run_name, workspace_root)
    else:
        warnings.append(f"Missing file: {manifest_path}")

    return RunData(
        run_name=run_name,
        run_dir=run_dir,
        workflow_versions=workflow_versions,
        workflow_params=workflow_params,
        fragger_params=fragger_params,
        fragger_rows=fragger_rows,
        manifest_entries=manifest_entries,
        parse_warnings=warnings,
    )


def to_unit_text(value: str, unit_code: str) -> str:
    unit_map = {"0": "Da", "1": "ppm"}
    unit = unit_map.get(unit_code, unit_code)
    return f"{value} {unit}".strip()


def parse_float(text: str) -> float | None:
    token = text.strip().split()[0] if text.strip() else ""
    if not token:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def nonzero_text(text: str) -> bool:
    parsed = parse_float(text)
    if parsed is None:
        return bool(text.strip())
    return parsed != 0.0


def is_non_default_value(text: str) -> bool:
    value = text.strip()
    if not value:
        return False
    lowered = value.lower()
    if lowered in {"0", "0.0", "false", "none", "null"}:
        return False
    parsed = parse_float(value)
    if parsed is not None and parsed == 0.0:
        return False
    return True


def fragger_comment_lookup(run: RunData) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in run.fragger_rows:
        if row.comment and row.key not in lookup:
            lookup[row.key] = row.comment
    for row in run.fragger_rows:
        if row.active and row.comment:
            lookup[row.key] = row.comment
    return lookup


def tooltip_label(key: str, comments: dict[str, str]) -> str:
    label = html.escape(key)
    comment = comments.get(key, "").strip()
    if comment:
        return f"<span class='has-tip' title='{html.escape(comment)}'>{label}</span>"
    return label


def get_modification_lists(run: RunData) -> tuple[list[str], list[str], list[str]]:
    active_variable: list[str] = []
    inactive_variable: list[str] = []
    fixed_mods: list[str] = []

    for row in run.fragger_rows:
        if row.key.startswith("variable_mod_"):
            label = f"{row.key} = {row.value}"
            if row.active:
                active_variable.append(label)
            else:
                inactive_variable.append(label)

    for key, value in sorted(run.fragger_params.items()):
        if not key.startswith("add_"):
            continue
        if nonzero_text(value):
            fixed_mods.append(f"{key} = {value}")

    return active_variable, inactive_variable, fixed_mods


def render_html_table(headers: list[str], rows: list[list[str]]) -> str:
    pieces = ["<table>", "<thead><tr>"]
    for head in headers:
        pieces.append(f"<th>{html.escape(head)}</th>")
    pieces.append("</tr></thead><tbody>")
    for row in rows:
        pieces.append("<tr>")
        for value in row:
            pieces.append(f"<td>{value}</td>")
        pieces.append("</tr>")
    pieces.append("</tbody></table>")
    return "".join(pieces)


def render_j2_template(template_name: str, **context: object) -> str:
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape
    except ImportError as e:
        raise RuntimeError(
            "Jinja2 is required for search report rendering. Install with: pip install jinja2"
        ) from e

    template_dir = pathlib.Path(__file__).resolve().parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(enabled_extensions=("html", "xml", "j2")),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = env.get_template(template_name)
    return template.render(**context)


def highlighted_methods(run: RunData) -> list[str]:
    fr = run.fragger_params
    wf = run.workflow_params

    data_type_map = {
        "0": "DDA",
        "1": "DIA",
        "2": "GPF-DIA",
        "3": "DDA+",
    }
    data_type = data_type_map.get(fr.get("data_type", ""), fr.get("data_type", ""))

    precursor_text = to_unit_text(
        fr.get("precursor_true_tolerance", ""),
        fr.get("precursor_true_units", ""),
    ).strip()
    fragment_text = to_unit_text(
        fr.get("fragment_mass_tolerance", ""),
        fr.get("fragment_mass_units", ""),
    ).strip()

    active_var_mods, _, fixed_mods = get_modification_lists(run)

    highlights = [
        f"Database FASTA: {wf.get('database.db-path', fr.get('database_name', 'n/a'))}",
        f"Search type: {data_type or 'n/a'}; enzyme: {fr.get('search_enzyme_name_1', 'n/a')}; missed cleavages: {fr.get('allowed_missed_cleavage_1', 'n/a')}",
        f"Precursor tolerance: {precursor_text or 'n/a'}; fragment tolerance: {fragment_text or 'n/a'}",
        f"IonQuant enabled: {wf.get('ionquant.run-ionquant', 'n/a')}; TMT-Integrator enabled: {wf.get('tmtintegrator.run-tmtintegrator', 'n/a')}",
        f"MSBooster enabled: {wf.get('msbooster.run-msbooster', 'n/a')}",
        f"Second enzyme setting: {fr.get('search_enzyme_name_2', 'n/a') or 'n/a'}",
        f"Active variable mods: {len(active_var_mods)}; active fixed mods: {len(fixed_mods)}",
    ]
    return highlights


def key_msfragger_rows(run: RunData) -> list[list[str]]:
    fr = run.fragger_params
    comments = fragger_comment_lookup(run)
    rows: list[list[str]] = []

    def add_row(key: str) -> None:
        if key not in fr:
            return
        rows.append([tooltip_label(key, comments), html.escape(fr.get(key, ""))])

    base_keys = [
        "calibrate_mass",
        "clip_nTerm_M",
        "num_enzyme_termini",
        "search_enzyme_name_1",
        "search_enzyme_cut_1",
        "allowed_missed_cleavage_1",
        "clear_mz_range",
        "digest_min_length",
        "digest_max_length",
        "digest_mass_range",
        "max_fragment_charge",
        "deisotope",
        "deneutralloss",
        "isotope_error",
        "fragment_ion_series",
        "remove_precursor_peak",
        "write_calibrated_mzml",
    ]
    for key in base_keys:
        add_row(key)

    if fr.get("ion_series_definitions", "").strip():
        add_row("ion_series_definitions")

    second_name = fr.get("search_enzyme_name_2", "").strip().lower()
    second_cut = fr.get("search_enzyme_cut_2", "").strip()
    second_enzyme_used = second_name not in {"", "null", "none"} or bool(second_cut)
    if second_enzyme_used:
        for key in ["search_enzyme_name_2", "search_enzyme_cut_2", "allowed_missed_cleavage_2"]:
            add_row(key)

    mass_offset_related = any(
        [
            is_non_default_value(fr.get("mass_offsets", "")),
            is_non_default_value(fr.get("mass_offsets_detailed", "")),
            is_non_default_value(fr.get("use_detailed_offsets", "")),
            is_non_default_value(fr.get("mass_diff_to_variable_mod", "")),
            is_non_default_value(fr.get("localize_delta_mass", "")),
        ]
    )
    if mass_offset_related:
        if is_non_default_value(fr.get("mass_offsets", "")):
            add_row("mass_offsets")
        if fr.get("mass_offsets_detailed", "").strip():
            add_row("mass_offsets_detailed")
        if is_non_default_value(fr.get("use_detailed_offsets", "")):
            add_row("use_detailed_offsets")
        if is_non_default_value(fr.get("mass_diff_to_variable_mod", "")):
            add_row("mass_diff_to_variable_mod")
        if is_non_default_value(fr.get("localize_delta_mass", "")):
            add_row("localize_delta_mass")
        if fr.get("delta_mass_exclude_ranges", "").strip():
            add_row("delta_mass_exclude_ranges")

    return rows


def diff_keys(run_data: list[RunData], attr: str) -> list[str]:
    union: set[str] = set()
    for run in run_data:
        mapping = getattr(run, attr)
        union.update(mapping.keys())

    differing: list[str] = []
    for key in sorted(union):
        values = []
        for run in run_data:
            mapping = getattr(run, attr)
            values.append(mapping.get(key, ""))
        if len(set(values)) > 1:
            differing.append(key)
    return differing


def write_mzml_tsv(path: pathlib.Path, run_data: list[RunData]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "run_name",
        "source_path",
        "file_name",
        "record",
        "group",
        "acquisition",
        "workspace_candidate",
        "exists_source_path",
        "exists_workspace_candidate",
        "status",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        for run in run_data:
            for entry in run.manifest_entries:
                writer.writerow(
                    [
                        entry.run_name,
                        entry.source_path,
                        entry.file_name,
                        entry.record,
                        entry.group,
                        entry.acquisition,
                        entry.workspace_candidate,
                        str(entry.exists_source_path).lower(),
                        str(entry.exists_workspace_candidate).lower(),
                        entry.status,
                    ]
                )


def build_html_report(
    run_data: list[RunData],
    out_html: pathlib.Path,
    out_mzml_tsv: pathlib.Path,
) -> str:
    now = dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")

    tool_names: set[str] = set()
    for run in run_data:
        tool_names.update(run.workflow_versions.keys())
    tool_rows: list[list[str]] = []
    for tool in sorted(tool_names):
        row = [html.escape(tool)]
        for run in run_data:
            row.append(html.escape(run.workflow_versions.get(tool, "")))
        tool_rows.append(row)

    summary_rows: list[list[str]] = []
    for run in run_data:
        missing_count = sum(1 for item in run.manifest_entries if item.status == "missing")
        total = len(run.manifest_entries)
        summary_rows.append(
            [
                html.escape(run.run_name),
                html.escape(str(run.run_dir)),
                html.escape(run.workflow_versions.get("FragPipe", "")),
                html.escape(run.workflow_versions.get("MSFragger", "")),
                html.escape(run.workflow_params.get("database.db-path", run.fragger_params.get("database_name", ""))),
                html.escape(f"{total}"),
                html.escape(f"{missing_count}"),
            ]
        )

    workflow_diff = diff_keys(run_data, "workflow_params")
    fragger_diff = diff_keys(run_data, "fragger_params")

    workflow_diff_rows: list[list[str]] = []
    for key in workflow_diff:
        row = [html.escape(key)]
        for run in run_data:
            row.append(html.escape(run.workflow_params.get(key, "")))
        workflow_diff_rows.append(row)

    fragger_diff_rows: list[list[str]] = []
    for key in fragger_diff:
        row = [html.escape(key)]
        for run in run_data:
            row.append(html.escape(run.fragger_params.get(key, "")))
        fragger_diff_rows.append(row)

    manifest_rows: list[list[str]] = []
    for run in run_data:
        for entry in run.manifest_entries:
            status_class = "missing" if entry.status == "missing" else "present"
            status_html = f"<span class='{status_class}'>{html.escape(entry.status)}</span>"
            manifest_rows.append(
                [
                    html.escape(entry.run_name),
                    html.escape(entry.file_name),
                    html.escape(entry.record),
                    html.escape(entry.group),
                    html.escape(entry.acquisition),
                    status_html,
                    html.escape(entry.source_path),
                    html.escape(entry.workspace_candidate),
                ]
            )

    warnings: list[str] = []
    for run in run_data:
        for warning in run.parse_warnings:
            warnings.append(f"{run.run_name}: {warning}")

    run_sections: list[str] = []
    for run in run_data:
        active_var_mods, _, fixed_mods = get_modification_lists(run)
        highlights = highlighted_methods(run)
        key_rows = key_msfragger_rows(run)
        param_rows = [[html.escape(k), html.escape(v)] for k, v in sorted(run.fragger_params.items())]
        workflow_rows = [[html.escape(k), html.escape(v)] for k, v in sorted(run.workflow_params.items())]

        section_parts = [
            f"<section><h2>Run: {html.escape(run.run_name)}</h2>",
            "<h3>Method Highlights</h3>",
            "<ul>",
            "".join(f"<li>{html.escape(item)}</li>" for item in highlights),
            "</ul>",
            "<h3>Key MSFragger Settings</h3>",
            "<p class='meta'>Hover parameter names for descriptions from fragger.params comments.</p>",
            render_html_table(["Parameter", "Value"], key_rows),
            "<h3>Modifications</h3>",
            f"<p><strong>Active variable mods ({len(active_var_mods)}):</strong></p>",
        ]
        if active_var_mods:
            section_parts.extend(
                [
                    "<ul>",
                    "".join(f"<li>{html.escape(item)}</li>" for item in active_var_mods),
                    "</ul>",
                ]
            )
        else:
            section_parts.append("<p>None detected.</p>")

        section_parts.append(
            f"<p><strong>Active fixed mods (non-zero) ({len(fixed_mods)}):</strong></p>"
        )
        if fixed_mods:
            section_parts.extend(
                [
                    "<ul>",
                    "".join(f"<li>{html.escape(item)}</li>" for item in fixed_mods),
                    "</ul>",
                ]
            )
        else:
            section_parts.append("<p>None detected.</p>")

        section_parts.extend(
            [
                "<details><summary>All active fragger.params key/value pairs</summary>",
                render_html_table(["Parameter", "Value"], param_rows),
                "</details>",
                "<details><summary>All fragpipe.workflow key/value pairs</summary>",
                render_html_table(["Parameter", "Value"], workflow_rows),
                "</details>",
                "</section>",
            ]
        )
        run_sections.append("".join(section_parts))

    context = {
        "generated_at": now,
        "input_dirs_text": ", ".join(str(run.run_dir) for run in run_data),
        "mzml_tsv_path": str(out_mzml_tsv),
        "run_summary_table_html": render_html_table(
            [
                "Run",
                "Run Directory",
                "FragPipe Version",
                "MSFragger Version",
                "Database FASTA",
                "Manifest Entries",
                "Missing mzML",
            ],
            summary_rows,
        ),
        "workflow_diff_count": len(workflow_diff_rows),
        "workflow_diff_table_html": render_html_table(
            ["Parameter", *[run.run_name for run in run_data]],
            workflow_diff_rows,
        ),
        "fragger_diff_count": len(fragger_diff_rows),
        "fragger_diff_table_html": render_html_table(
            ["Parameter", *[run.run_name for run in run_data]],
            fragger_diff_rows,
        ),
        "manifest_table_html": render_html_table(
            [
                "Run",
                "File",
                "Record",
                "Group",
                "Acquisition",
                "Status",
                "Manifest Source Path",
                "Workspace Candidate Path",
            ],
            manifest_rows,
        ),
        "run_sections_html": "".join(run_sections),
        "tool_versions_table_html": render_html_table(
            ["Tool", *[run.run_name for run in run_data]],
            tool_rows,
        ),
        "warnings": warnings,
    }
    return render_j2_template("search_report.html.j2", **context)


def resolve_run_dirs(
    run_dirs: list[str | pathlib.Path] | None,
    workspace_root: pathlib.Path,
) -> list[pathlib.Path]:
    raw_run_dirs = run_dirs if run_dirs else [str(path) for path in DEFAULT_RUN_DIRS]
    resolved: list[pathlib.Path] = []
    for item in raw_run_dirs:
        path = pathlib.Path(item).expanduser()
        if not path.is_absolute():
            path = (workspace_root / path).resolve()
        resolved.append(path)
    return resolved


def build_fragpipe_methods_report(
    *,
    run_dirs: list[str | pathlib.Path] | None = None,
    out_html: str | pathlib.Path = DEFAULT_OUT_HTML,
    out_mzml_tsv: str | pathlib.Path = DEFAULT_OUT_MZML_TSV,
    workspace_root: str | pathlib.Path = ".",
) -> SearchReportOutputs:
    workspace_root_path = pathlib.Path(workspace_root).expanduser().resolve()
    out_html_path = pathlib.Path(out_html).expanduser().resolve()
    out_mzml_tsv_path = pathlib.Path(out_mzml_tsv).expanduser().resolve()
    resolved_run_dirs = resolve_run_dirs(run_dirs, workspace_root_path)

    run_data = [load_run(path, workspace_root_path) for path in resolved_run_dirs]

    write_mzml_tsv(out_mzml_tsv_path, run_data)
    html_report = build_html_report(run_data, out_html_path, out_mzml_tsv_path)
    out_html_path.parent.mkdir(parents=True, exist_ok=True)
    out_html_path.write_text(html_report, encoding="utf-8")

    missing_total = sum(
        1
        for run in run_data
        for entry in run.manifest_entries
        if entry.status == "missing"
    )
    manifest_total = sum(len(run.manifest_entries) for run in run_data)
    return SearchReportOutputs(
        out_html=out_html_path,
        out_mzml_tsv=out_mzml_tsv_path,
        manifest_total=manifest_total,
        missing_total=missing_total,
    )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = build_fragpipe_methods_report(
        run_dirs=args.run_dirs,
        out_html=args.out_html,
        out_mzml_tsv=args.out_mzml_tsv,
        workspace_root=args.workspace_root,
    )
    print(f"Wrote HTML report: {result.out_html}")
    print(f"Wrote mzML tracking TSV: {result.out_mzml_tsv}")
    print(
        "Manifest entries: "
        f"{result.manifest_total}; missing mzML entries: {result.missing_total}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
