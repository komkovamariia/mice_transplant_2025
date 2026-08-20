import json
import re
from pathlib import Path

import pandas as pd
import pytest

from scripts import run_analysis


ROOT = Path(__file__).resolve().parents[1]


def test_restored_notebook_retains_historical_scope_and_is_clean_source():
    notebook = json.loads((ROOT / "venn_original.ipynb").read_text(encoding="utf-8"))
    cells = notebook["cells"]
    code_cells = [cell for cell in cells if cell["cell_type"] == "code"]
    joined = "\n".join("".join(cell.get("source", [])) for cell in cells)

    assert len(cells) >= 40
    assert len(code_cells) >= 30
    assert not re.search(r"[\u0400-\u04FF]", joined)
    assert "deseq" not in joined.lower()
    assert all(not cell.get("outputs") for cell in code_cells)
    assert all(cell.get("execution_count") is None for cell in code_cells)
    assert joined.count("plt.show(") >= 20
    assert "plt.savefig = _article_savefig" in joined
    assert "plt.show = _article_show" in joined


def test_output_verification_rejects_missing_stratum_artifacts(tmp_path, monkeypatch):
    monkeypatch.setattr(run_analysis, "repo_root", lambda: tmp_path)
    with pytest.raises(RuntimeError, match="missing"):
        run_analysis.verify_first_approach_outputs(["cd4_thymus"])


def test_output_verification_accepts_a_complete_single_stratum(tmp_path, monkeypatch):
    monkeypatch.setattr(run_analysis, "repo_root", lambda: tmp_path)
    stem = run_analysis.STRATUM_BY_KEY["cd4_thymus"]
    result_dir = tmp_path / "results" / "01_set_count" / "cd4_thymus"
    figure_dir = tmp_path / "figures" / "01_set_count" / "cd4_thymus"
    audit_dir = tmp_path / "audit_runs"
    log_dir = tmp_path / "logs"
    for directory in (result_dir, figure_dir, audit_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    (audit_dir / f"{stem}.executed.ipynb").write_text("{}", encoding="utf-8")
    (log_dir / f"{stem}.log").write_text("complete\n", encoding="utf-8")
    (result_dir / "run_summary.json").write_text("{}\n", encoding="utf-8")
    (result_dir / "trav_ranking.csv").write_text("v_gene\nTRAV1\n", encoding="utf-8")
    (result_dir / "distinctive_trav.csv").write_text(
        "v_gene\nTRAV1\n", encoding="utf-8"
    )
    png = figure_dir / "figure.png"
    pdf = figure_dir / "figure.pdf"
    png.write_bytes(b"png")
    pdf.write_bytes(b"pdf")
    pd.DataFrame(
        [
            {
                "png": str(png.relative_to(tmp_path)),
                "pdf": str(pdf.relative_to(tmp_path)),
            }
        ]
    ).to_csv(result_dir / "figure_manifest.csv", index=False)

    run_analysis.verify_first_approach_outputs(["cd4_thymus"])
