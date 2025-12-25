from pathlib import Path
import os

import pytest

from tackle.volcanoplot import safe_path_with_ext


def test_safe_path_keeps_short_names(tmp_path):
    base = tmp_path / "volcano_output"
    expected = base.with_name(base.name + ".tsv")

    result = safe_path_with_ext(str(base), ".tsv")

    assert result == str(expected)


@pytest.mark.parametrize("ext", [".tsv", "tsv"])
def test_safe_path_appends_single_extension(tmp_path, ext):
    base = tmp_path / "volcano_output"

    result = safe_path_with_ext(str(base), ext)

    assert Path(result).suffixes == [".tsv"]


def test_safe_path_truncates_long_component(tmp_path):
    long_base = "gene_set_" + "x" * 270
    path = tmp_path / long_base

    result = safe_path_with_ext(str(path), ".tsv")

    final_component = Path(result).name
    limit = os.pathconf(str(tmp_path), "PC_NAME_MAX")

    assert final_component.endswith(".tsv")
    assert len(final_component.encode("utf-8")) <= limit


def test_safe_path_preserves_dots_in_base_name(tmp_path):
    path = tmp_path / "sig_0.123"

    result = safe_path_with_ext(str(path), ".tsv")

    assert Path(result).suffixes == [".123", ".tsv"]
