import zipfile

import nbformat
import pytest

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


def test_resolve_strata_preserves_the_prespecified_eight_run_order():
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
