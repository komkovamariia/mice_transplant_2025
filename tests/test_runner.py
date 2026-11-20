import sys
import zipfile

import matplotlib.pyplot as plt
import nbformat
import pytest

from src.figures import save_figure
from scripts import run_analysis
from scripts.run_analysis import (
    STRATUM_RUNS,
    _cell_description,
    archive_first_approach_figures,
    prepare_first_approach_output_directories,
    resolve_sequence,
    resolve_strata,
)


def test_resolve_sequence_supports_individual_combined_and_alias_values():
    assert resolve_sequence("1") == ["1"]
    assert resolve_sequence("1,3") == ["1", "3"]
    assert resolve_sequence("all") == ["1", "2", "3", "4"]
    assert resolve_sequence("compare") == ["4"]
    assert resolve_sequence("1,1,2") == ["1", "2"]


def test_resolve_sequence_rejects_unknown_values():
    with pytest.raises(SystemExit, match="Unknown approach"):
        resolve_sequence("5")


def test_resolve_strata_preserves_the_nine_stratum_order():
    expected = [key for _, key in STRATUM_RUNS]
    assert resolve_strata("all") == expected
    assert resolve_strata("cd4_thymus,thymus_combined") == [
        "cd4_thymus",
        "thymus_combined",
    ]


def test_resolve_strata_rejects_unknown_values():
    with pytest.raises(SystemExit, match="Unknown stratum"):
        resolve_strata("whole_mouse")


def test_cell_description_prefers_explicit_metadata():
    cell = nbformat.v4.new_code_cell(
        "print('fallback')",
        metadata={"analysis_step": "Fit stratum model"},
    )
    assert _cell_description(cell) == "Fit stratum model"


def test_cell_description_uses_first_executable_line():
    cell = nbformat.v4.new_code_cell("\n\nresult = run_stratum()\n")
    assert _cell_description(cell) == "result = run_stratum()"


def test_prepare_and_archive_first_approach_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr("scripts.run_analysis.repo_root", lambda: tmp_path)
    figure_root = tmp_path / "figures" / "01_set_count"
    selected = figure_root / "cd4_thymus"
    retained = figure_root / "cd8_thymus"
    selected.mkdir(parents=True)
    retained.mkdir(parents=True)
    (selected / "stale.png").write_bytes(b"stale")
    (retained / "retained.png").write_bytes(b"retained")
    stale_results = tmp_path / "results" / "01_set_count" / "cd4_thymus"
    retained_results = tmp_path / "results" / "01_set_count" / "cd8_thymus"
    stale_results.mkdir(parents=True)
    retained_results.mkdir(parents=True)
    (stale_results / "g2_stale.csv").write_text("stale\n", encoding="utf-8")
    (retained_results / "retained.csv").write_text("retained\n", encoding="utf-8")
    (tmp_path / "figures" / "01_set_count.zip").write_bytes(b"old")

    prepare_first_approach_output_directories(["cd4_thymus"])

    assert not selected.exists()
    assert retained.is_dir()
    assert not stale_results.exists()
    assert retained_results.is_dir()
    assert not (tmp_path / "figures" / "01_set_count.zip").exists()

    archive_path = archive_first_approach_figures()
    assert archive_path.is_file()
    with zipfile.ZipFile(archive_path) as archive:
        assert archive.namelist() == ["01_set_count/cd8_thymus/retained.png"]


def test_adaptive_outputs_preserve_fixed_results(tmp_path, monkeypatch):
    monkeypatch.setattr("scripts.run_analysis.repo_root", lambda: tmp_path)
    fixed = tmp_path / "figures" / "01_set_count.zip"
    fixed.parent.mkdir()
    fixed.write_bytes(b"fixed")
    result = tmp_path / "results" / "01_set_count" / "all_combined" / "summary.csv"
    result.parent.mkdir(parents=True)
    result.write_text("fixed\n")
    adaptive = tmp_path / "figures" / "01_set_count_adaptive" / "all_combined"
    adaptive.mkdir(parents=True)

    prepare_first_approach_output_directories(["all_combined"], output_name="01_set_count_adaptive")
    assert not adaptive.exists()
    adaptive.mkdir()
    (adaptive / "plot.png").write_bytes(b"adaptive")
    archive = archive_first_approach_figures(output_name="01_set_count_adaptive")
    assert archive.name == "01_set_count_adaptive.zip"
    assert fixed.read_bytes() == b"fixed"
    assert result.read_text() == "fixed\n"


@pytest.mark.parametrize(("options", "mode"), [([], "adaptive"), (["--v-thresholds", "fixed"], "fixed")])
def test_standard_run_uses_dynamic_cutoffs_unless_fixed_is_requested(tmp_path, monkeypatch, options, mode):
    monkeypatch.setattr(sys, "argv", ["run_analysis.py", "--strata", "all_combined", *options])
    monkeypatch.setattr(run_analysis, "repo_root", lambda: tmp_path)
    executed = []
    monkeypatch.setattr(run_analysis, "execute_notebook", lambda *args, **kwargs: executed.append(kwargs))
    verified = []
    monkeypatch.setattr(run_analysis, "verify_first_approach_outputs", lambda *args, **kwargs: verified.append(kwargs))
    monkeypatch.setattr(run_analysis, "archive_first_approach_figures", lambda **_: tmp_path / "figures.zip")
    assert run_analysis.main() == 0
    assert [run["environment"]["MICE_TCR_V_THRESHOLDS"] for run in executed] == [mode]
    assert [run["threshold_mode"] for run in verified] == [mode]


def test_publication_figure_output_is_png_only_and_chain_specific(tmp_path, monkeypatch):
    monkeypatch.setattr("src.figures.repository_root", lambda: tmp_path)
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot([0, 1], [0, 1])
    ax.set_title("TRA test figure")

    path = save_figure(
        fig,
        "01_set_count",
        "all_combined",
        "test_figure",
        chain="TRA",
    )

    expected = (
        tmp_path
        / "figures"
        / "01_set_count"
        / "all_combined"
        / "TRA"
        / "test_figure.png"
    )
    assert path == expected
    assert path.is_file()
    assert not path.with_suffix(".pdf").exists()
