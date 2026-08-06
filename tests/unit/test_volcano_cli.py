from click.testing import CliRunner

from tackle.main import volcano


def _invoke_volcano(monkeypatch, options):
    captured = {}

    def fake_volcanoplot(*args, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr("tackle.volcanoplot.volcanoplot", fake_volcanoplot)
    result = CliRunner().invoke(volcano, options, obj={})
    assert result.exit_code == 0, result.output
    return captured


def test_volcano_side_label_counts_default_to_no_override(monkeypatch):
    captured = _invoke_volcano(monkeypatch, [])

    assert captured["number"] == 35
    assert captured["number_left"] is None
    assert captured["number_right"] is None
    assert captured["label_size_by"] == "fixed"
    assert captured["label_size_min"] == 2.4
    assert captured["label_size_max"] == 4.0


def test_volcano_side_label_count_aliases_are_forwarded(monkeypatch):
    captured = _invoke_volcano(
        monkeypatch,
        ["-n", "20", "-n-left", "3", "--number-right", "7"],
    )

    assert captured["number"] == 20
    assert captured["number_left"] == 3
    assert captured["number_right"] == 7


def test_volcano_side_label_counts_reject_negative_values():
    result = CliRunner().invoke(volcano, ["--number-left", "-1"], obj={})

    assert result.exit_code != 0
    assert "x>=0" in result.output


def test_volcano_density_label_size_options_are_forwarded(monkeypatch):
    captured = _invoke_volcano(
        monkeypatch,
        [
            "--label-size-by",
            "density",
            "--label-size-min",
            "2.0",
            "--label-size-max",
            "5.0",
            "--number-by",
            "hybrid",
        ],
    )

    assert captured["label_size_by"] == "density"
    assert captured["label_size_min"] == 2.0
    assert captured["label_size_max"] == 5.0
    assert captured["number_by"] == "hybrid"


def test_volcano_density_label_size_range_must_be_ordered():
    result = CliRunner().invoke(
        volcano,
        ["--label-size-min", "5", "--label-size-max", "2"],
        obj={},
    )

    assert result.exit_code != 0
    assert "must be less than or equal" in result.output
