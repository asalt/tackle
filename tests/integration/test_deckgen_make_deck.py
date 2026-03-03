import os
from pathlib import Path

import pandas as pd
import pytest


def _write_png(path: Path, text: str = "img") -> None:
    import matplotlib

    try:
        matplotlib.use("Agg", force=True)
    except Exception:
        pass
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(2, 1))
    fig.text(0.5, 0.5, text, ha="center", va="center")
    fig.savefig(path, dpi=120)
    plt.close(fig)


def test_build_volcano_table_assets_writes_png_and_pdf(tmp_path):
    from tackle.deckgen import build_volcano_table_assets

    base = tmp_path
    (base / "volcano").mkdir(parents=True, exist_ok=True)

    pd.DataFrame(
        {
            "GeneID": ["g1", "g2", "g3", "g4"],
            "GeneSymbol": ["AAA", "BBB", "CCC", "DDD"],
            "log2_FC": [1.5, -2.2, 0.2, -0.1],
            "pAdj": [0.001, 0.02, 0.5, 0.9],
            "GeneDescription": ["alpha", "beta", "gamma", "delta"],
        }
    ).to_csv(base / "volcano" / "condA_by_condB.tsv", sep="\t", index=False)

    out_dir = base / "report" / "deck" / "tables"
    assets = build_volcano_table_assets(
        base_dir=str(base),
        out_dir=str(out_dir),
        topn=2,
        direction="both",
        p_cutoff=0.05,
        formats=(".png", ".pdf"),
        force=False,
    )

    assert len(assets) == 1
    asset = assets[0]
    assert asset.png_path and os.path.exists(asset.png_path)
    assert asset.pdf_path and os.path.exists(asset.pdf_path)

    # Basic sanity: files are non-empty
    assert Path(asset.png_path).stat().st_size > 0
    assert Path(asset.pdf_path).stat().st_size > 0


def test_collect_plot_image_assets_discovers_expected_categories(tmp_path):
    from tackle.deckgen import collect_plot_image_assets

    _write_png(tmp_path / "volcano" / "groupA_volcano_plot.png", "volcano")
    _write_png(tmp_path / "pca" / "pcaplot_groupA.png", "pca")
    _write_png(tmp_path / "metrics" / "metrics_dist_groupA.png", "metrics")
    _write_png(tmp_path / "cluster" / "clustermap_groupA.png", "cluster")
    _write_png(tmp_path / "topdiff" / "cluster" / "clustermap_topdiff_groupA.png", "topdiff")
    _write_png(tmp_path / "report" / "deck" / "tables" / "generated.png", "generated")
    _write_png(tmp_path / "misc" / "other_plot.png", "other")

    assets = collect_plot_image_assets(
        base_dir=str(tmp_path),
        include_kinds=("volcano", "pca", "metrics", "cluster", "topdiff-cluster"),
        exclude_dirs=(str(tmp_path / "report" / "deck"),),
    )

    categories = {asset.category for asset in assets}
    assert categories == {"volcano", "pca", "metrics", "cluster", "topdiff-cluster"}
    # Generated deck assets should not be re-ingested.
    assert all("report/deck" not in asset.source_relpath.replace("\\", "/") for asset in assets)
    # Unclassified misc image should not be included.
    assert all("misc/other_plot.png" != asset.source_relpath.replace("\\", "/") for asset in assets)


@pytest.mark.skipif(
    __import__("importlib.util").util.find_spec("pptx") is None,
    reason="python-pptx not installed",
)
def test_build_pptx_deck_writes_file(tmp_path):
    from tackle.deckgen import TableAsset, build_pptx_deck

    # Create a tiny placeholder PNG to embed.
    png = tmp_path / "table.png"
    import matplotlib

    try:
        matplotlib.use("Agg", force=True)
    except Exception:
        pass
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(2, 1))
    fig.text(0.5, 0.5, "hi", ha="center", va="center")
    fig.savefig(png, dpi=150)
    plt.close(fig)

    assets = [
        TableAsset(
            title="Slide 1",
            source_relpath="volcano/x.tsv",
            source_path=str(tmp_path / "x.tsv"),
            table_rows=1,
            table_cols=1,
            png_path=str(png),
            pdf_path=None,
        )
    ]

    out = tmp_path / "deck.pptx"
    res = build_pptx_deck(out_path=str(out), title="Test", assets=assets)
    assert res == str(out)
    assert out.exists()
    assert out.stat().st_size > 0
