import nbformat
import pytest

from scripts.run_analysis import (
    STRATUM_RUNS,
    _cell_description,
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
