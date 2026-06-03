from __future__ import annotations

import json
import os
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from jinja2 import Environment, FileSystemLoader

from .utils import get_gaussian_imputation_defaults


@dataclass(frozen=True)
class MethodsReferenceOutputs:
    markdown_path: Path
    context_path: Path


REFERENCE_REGISTRY: Mapping[str, Mapping[str, str]] = {
    "limma": {
        "label": "limma moderated linear models",
        "pmid": "25605792",
        "doi": "10.1093/nar/gkv007",
        "url": "https://pubmed.ncbi.nlm.nih.gov/25605792/",
    },
    "benjamini_hochberg": {
        "label": "Benjamini-Hochberg false discovery rate correction",
        "pmid": "",
        "doi": "10.1111/j.2517-6161.1995.tb02031.x",
        "url": "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
    },
    "perseus": {
        "label": "Perseus-style proteomics imputation/workflow",
        "pmid": "27348712",
        "doi": "10.1038/nmeth.3901",
        "url": "https://pubmed.ncbi.nlm.nih.gov/27348712/",
    },
    "gsea_mootha_2003": {
        "label": "Gene set enrichment analysis, original pathway analysis application",
        "pmid": "12808457",
        "doi": "10.1038/ng1180",
        "url": "https://pubmed.ncbi.nlm.nih.gov/12808457/",
    },
    "gsea_subramanian_2005": {
        "label": "Gene Set Enrichment Analysis",
        "pmid": "16199517",
        "doi": "10.1073/pnas.0506580102",
        "url": "https://pubmed.ncbi.nlm.nih.gov/16199517/",
    },
    "msigdb_hallmark": {
        "label": "Molecular Signatures Database hallmark gene set collection",
        "pmid": "26771021",
        "doi": "10.1016/j.cels.2015.12.004",
        "url": "https://pubmed.ncbi.nlm.nih.gov/26771021/",
    },
    "fgsea": {
        "label": "fgsea fast gene set enrichment analysis",
        "pmid": "",
        "doi": "10.1101/060012",
        "url": "https://bioconductor.org/packages/fgsea",
    },
    "mitocarta3": {
        "label": "MitoCarta3.0 mitochondrial protein inventory",
        "pmid": "33174596",
        "doi": "10.1093/nar/gkaa1011",
        "url": "https://pubmed.ncbi.nlm.nih.gov/33174596/",
    },
    "gene_ontology": {
        "label": "Gene Ontology resource",
        "pmid": "33290552",
        "doi": "10.1093/nar/gkaa1113",
        "url": "https://pubmed.ncbi.nlm.nih.gov/33290552/",
    },
    "kegg": {
        "label": "KEGG database resource",
        "pmid": "33125081",
        "doi": "10.1093/nar/gkaa970",
        "url": "https://pubmed.ncbi.nlm.nih.gov/33125081/",
    },
    "reactome": {
        "label": "Reactome pathway knowledgebase",
        "pmid": "34788843",
        "doi": "10.1093/nar/gkab1028",
        "url": "https://pubmed.ncbi.nlm.nih.gov/34788843/",
    },
    "uniprot": {
        "label": "UniProt Universal Protein Knowledgebase",
        "pmid": "36408920",
        "doi": "10.1093/nar/gkac1052",
        "url": "https://pubmed.ncbi.nlm.nih.gov/36408920/",
    },
}

REFERENCE_ORDER: Tuple[str, ...] = (
    "limma",
    "benjamini_hochberg",
    "perseus",
    "gsea_mootha_2003",
    "gsea_subramanian_2005",
    "msigdb_hallmark",
    "fgsea",
    "mitocarta3",
    "gene_ontology",
    "kegg",
    "reactome",
    "uniprot",
)

ANNOTATION_TOKENS: Mapping[str, Tuple[str, ...]] = {
    "mitocarta3": ("MitoCarta", "MitoCarta_Pathways"),
    "gene_ontology": ("Gene Ontology", "GeneOntology", "GO_", "GOBP", "GOCC", "GOMF"),
    "kegg": ("KEGG", "curated.CP.KEGG"),
    "reactome": ("Reactome", "REACTOME", "curated.CP.Reactome"),
    "uniprot": ("UniProt", "Uniprot", "uniprot"),
}


def _display_path(path: Path) -> str:
    try:
        rel = os.path.relpath(path.resolve(), Path.cwd().resolve()).replace("\\", "/")
    except Exception:
        return str(path)
    if rel == ".":
        return "./"
    return rel if rel.startswith(".") else f"./{rel}"


def _safe_json_load(path: Path) -> Optional[Dict[str, Any]]:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return obj if isinstance(obj, dict) else None


def _glob_files(root: Path, pattern: str, *, limit: int = 200) -> List[Path]:
    files: List[Path] = []
    report_html = root / "report" / "html"
    for path in root.glob(pattern):
        if len(files) >= limit:
            break
        if not path.is_file():
            continue
        try:
            path.relative_to(report_html)
            continue
        except ValueError:
            pass
        files.append(path)
    return sorted(files)


def _has_any_path(root: Path, patterns: Sequence[str]) -> bool:
    for pattern in patterns:
        for _ in root.glob(pattern):
            return True
    return False


def _read_header(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return handle.readline().strip()
    except Exception:
        return ""


def _reference_list(keys: Iterable[str]) -> List[Dict[str, str]]:
    requested = set(keys)
    ordered_keys = [key for key in REFERENCE_ORDER if key in requested]
    ordered_keys.extend(sorted(requested.difference(ordered_keys)))
    refs: List[Dict[str, str]] = []
    for key in ordered_keys:
        ref = REFERENCE_REGISTRY.get(key)
        if not ref:
            continue
        entry = {"key": key}
        entry.update({str(k): str(v) for k, v in ref.items()})
        refs.append(entry)
    return refs


def _summarize_imputation_contexts(
    contexts: Sequence[Mapping[str, Any]],
) -> Tuple[List[str], Set[str], List[Dict[str, Any]]]:
    summaries: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
    refs: Set[str] = set()

    for ctx in contexts:
        if not bool(ctx.get("impute_missing_values")):
            continue
        backend = str(ctx.get("imputation_backend") or "unknown").strip().lower()
        if backend == "gaussian":
            method = str(ctx.get("gaussian_method") or "legacy").strip().lower()
            try:
                defaults = get_gaussian_imputation_defaults(method)
            except Exception:
                defaults = {}
            label_bits = [f"backend=gaussian", f"method={method}"]
            if defaults.get("downshift") is not None:
                label_bits.append(f"downshift={defaults.get('downshift')}")
            if defaults.get("effective_width") is not None:
                label_bits.append(f"effective_width={defaults.get('effective_width')}")
            elif defaults.get("scale") is not None:
                label_bits.append(f"scale={defaults.get('scale')}")
            key = "|".join(label_bits)
            refs.add("perseus")
        elif backend == "lupine":
            mode = str(ctx.get("lupine_mode") or "local").strip().lower()
            label_bits = [f"backend=lupine", f"mode={mode}"]
            key = "|".join(label_bits)
        else:
            label_bits = [f"backend={backend or 'unknown'}"]
            key = "|".join(label_bits)

        item = summaries.setdefault(
            key,
            {
                "settings": label_bits,
                "count": 0,
                "run_ids": [],
            },
        )
        item["count"] += 1
        run_id = str(ctx.get("run_id") or "").strip()
        if run_id and run_id not in item["run_ids"]:
            item["run_ids"].append(run_id)

    lines: List[str] = []
    for item in summaries.values():
        settings = ", ".join(item["settings"])
        line = f"Missing-value imputation context detected: {settings}; replay contexts={item['count']}."
        if item["run_ids"]:
            shown = ", ".join(item["run_ids"][:5])
            if len(item["run_ids"]) > 5:
                shown = f"{shown}, ..."
            line = f"{line} Run IDs: {shown}."
        lines.append(line)

    return lines, refs, list(summaries.values())


def _detect_annotation_references(tsv_files: Sequence[Path]) -> Tuple[List[str], Set[str]]:
    hits: "OrderedDict[str, List[str]]" = OrderedDict()
    refs: Set[str] = set()
    for path in tsv_files:
        header = _read_header(path)
        if not header:
            continue
        for ref_key, tokens in ANNOTATION_TOKENS.items():
            if any(token in header for token in tokens):
                refs.add(ref_key)
                hits.setdefault(ref_key, []).append(path.name)

    lines: List[str] = []
    for ref_key, files in hits.items():
        label = REFERENCE_REGISTRY.get(ref_key, {}).get("label", ref_key)
        shown = ", ".join(files[:5])
        if len(files) > 5:
            shown = f"{shown}, ..."
        lines.append(f"{label} annotations detected in exported table headers: {shown}.")
    return lines, refs


def build_methods_references_context(
    *,
    base_dir: Path,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    root = Path(base_dir).expanduser().resolve()
    refs: Set[str] = set()
    sections: List[Dict[str, Any]] = []

    volcano_tsvs = _glob_files(root, "volcano/**/*.tsv")
    replay_context_paths = _glob_files(root, "volcano/**/limma_replay_context.json")
    replay_context_paths.extend(_glob_files(root, "report/rmd/**/limma_replay_context.json"))
    replay_contexts = [
        ctx for ctx in (_safe_json_load(path) for path in replay_context_paths) if ctx is not None
    ]

    if volcano_tsvs or replay_contexts:
        refs.update({"limma", "benjamini_hochberg"})
        lines = [
            "Differential-analysis outputs were detected from volcano TSVs and/or limma replay contexts.",
            "Use the exact contrast formula, imputation settings, and stored replay matrix from the corresponding replay context when preparing final methods text.",
        ]
        sections.append(
            {
                "key": "differential_analysis",
                "title": "Differential analysis",
                "lines": lines,
            }
        )

    imputation_lines, imputation_refs, imputation_summary = _summarize_imputation_contexts(
        replay_contexts
    )
    if imputation_lines:
        refs.update(imputation_refs)
        sections.append(
            {
                "key": "missing_value_imputation",
                "title": "Missing-value imputation",
                "lines": [
                    "Missing-value imputation was detected in replay context metadata.",
                    "The stored imputed matrix remains the authoritative replay input; the settings below document the method identity and parameters.",
                    *imputation_lines,
                ],
            }
        )

    has_gsea = _has_any_path(
        root,
        (
            "gsea/**",
            "GSEA/**",
            "ssGSEA/**",
            "ssgsea/**",
            "enrichment/**",
            "replot_gsea/**",
        ),
    )
    if has_gsea:
        refs.update(
            {
                "gsea_mootha_2003",
                "gsea_subramanian_2005",
                "msigdb_hallmark",
                "fgsea",
            }
        )
        sections.append(
            {
                "key": "gene_set_enrichment",
                "title": "Gene set enrichment",
                "lines": [
                    "Gene-set enrichment outputs were detected.",
                    "Final methods text should name the exact gene-set collection and ranking statistic used for the run.",
                ],
            }
        )

    annotation_tsvs = []
    annotation_tsvs.extend(_glob_files(root, "export/**/*.tsv", limit=150))
    annotation_tsvs.extend(volcano_tsvs[:150])
    annotation_lines, annotation_refs = _detect_annotation_references(annotation_tsvs)
    if annotation_lines:
        refs.update(annotation_refs)
        sections.append(
            {
                "key": "annotation_sources",
                "title": "Annotation sources",
                "lines": [
                    "Annotation-like columns were detected in exported tables.",
                    "Only cite annotation sources that are actually present in the final reported outputs.",
                    *annotation_lines,
                ],
            }
        )

    if not sections:
        sections.append(
            {
                "key": "general",
                "title": "General notes",
                "lines": [
                    "No volcano, imputation, enrichment, or known annotation-source outputs were detected from the scanned analysis directory.",
                    "This document can still be used as a methods/reference scaffold after additional outputs are generated.",
                ],
            }
        )

    return {
        "title": title or "Tackle methods and references",
        "generated_at": datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"),
        "analysis_base_dir": _display_path(root),
        "sections": sections,
        "references": _reference_list(refs),
        "reference_count": len(refs),
        "detected": {
            "volcano_tsv_count": len(volcano_tsvs),
            "replay_context_count": len(replay_contexts),
            "imputation_contexts": imputation_summary,
            "gsea_outputs_detected": bool(has_gsea),
            "annotation_tsv_count": len(annotation_tsvs),
        },
        "notes": [
            "Generated Markdown is intentionally separate from data tables to avoid adding provenance-only columns to analytical outputs.",
            "Render to PDF with pandoc if desired: pandoc methods_references.md -o methods_references.pdf",
        ],
    }


def render_methods_references_markdown(context: Mapping[str, Any]) -> str:
    template_dir = Path(__file__).resolve().parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=False,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = env.get_template("methods_references.md.j2")
    return template.render(**context)


def write_methods_references(
    *,
    base_dir: Path,
    out_dir: Path,
    title: Optional[str] = None,
) -> MethodsReferenceOutputs:
    out_root = Path(out_dir).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    context = build_methods_references_context(base_dir=Path(base_dir), title=title)

    markdown_path = out_root / "methods_references.md"
    context_path = out_root / "methods_references.json"
    markdown_path.write_text(render_methods_references_markdown(context), encoding="utf-8")
    context_path.write_text(
        json.dumps(context, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return MethodsReferenceOutputs(markdown_path=markdown_path, context_path=context_path)
