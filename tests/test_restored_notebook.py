import ast
import hashlib
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
    assert joined.count("plt.show(") >= 10
    assert "plt.savefig = _article_savefig" in joined
    assert "plt.show = _article_show" in joined
    assert "sns.set_theme(style=\"white\"" in joined
    assert 'DOWNSAMPLE_DEPTH = int(os.environ.get("MICE_TCR_DOWNSAMPLE_UMI", "15000"))' in joined
    assert 'functionality="f"' in joined
    assert "by_umi=True" in joined
    assert "downsample_group_tables(" in joined
    assert 'globals()[f"raw_ct_{group}_aaV_tra"]' in joined
    assert 'globals()[f"ct_{group}_aaV_tra"]' in joined
    assert 'globals()[f"raw_ct_{group}_aaV_tra"] for group in ("g1", "g2", "g5", "g6")' in joined
    assert "RAREFACTION_SEEDS" in joined
    assert "BALANCE_SEEDS" in joined
    assert "collect_pairwise_similarity_boxplot_table(results_v_js)" in joined
    assert "top_n_labels=10" in joined
    assert "_adaptive_annotation_color" in joined
    assert "def plot_v_region_heatmap(" not in joined
    assert joined.count(
        'save_path=str(FIGURE_DIR / "TRA_g1_v_segment_localization.png")'
    ) == 1
    assert joined.count(
        'save_path=str(FIGURE_DIR / "TRB_g1_v_segment_localization.png")'
    ) == 1
    assert 'region_order = ["100", "110", "101", "111"]' in joined
    assert "np.log2(len(target_region_order))" in joined
    assert "min_full_freq_pct=0.6" in joined
    assert "min_retained_pct=66" in joined
    assert "impact_score" not in joined
    assert "specific_score" not in joined
    assert 'groups = ["g1", "g5"]' in joined
    assert 'GROUPS_TO_RUN = ["g1", "g5", "g6"]' in joined
    assert "cmp_g2_tra_all" not in joined
    assert "cmp_g2_all" not in joined
    assert "g2_remaining_aaV_tra" not in joined
    assert "g2_remaining_aaV_trb" not in joined
    assert "double_triple_g2" not in joined
    assert "_heatmaps.png" not in joined
    assert "_scatter.png" not in joined
    assert "_boxplot.png" not in joined
    assert "mouse_venn_panel" not in joined
    assert "distinctive_trav_across_eight_strata" not in joined
    assert "plt.subplots(1, 2" not in joined
    assert "plt.subplots(2, 2" not in joined
    assert "TRA g1: V-segment representation before and after background subtraction" in joined
    assert "TRA g1: V-segment localization across g1-containing overlap regions" in joined
    assert "TRA g1: Fisher and edgeR effect concordance" in joined


def _load_notebook_function(function_name):
    notebook = json.loads((ROOT / "venn_original.ipynb").read_text(encoding="utf-8"))
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell.get("source", []))
        tree = ast.parse(source)
        matches = [
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == function_name
        ]
        if matches:
            module = ast.Module(body=[matches[0]], type_ignores=[])
            namespace = {"np": __import__("numpy"), "pd": pd}
            exec(compile(module, "venn_original.ipynb", "exec"), namespace)
            return namespace[function_name]
    raise AssertionError(f"Function {function_name} was not found")


def test_article_score_and_entropy_use_retention_and_four_g1_regions():
    make_v_region_summary = _load_notebook_function("make_v_region_summary")
    make_v_profile = _load_notebook_function("make_v_profile")
    select_bubble_genes = _load_notebook_function("get_bubble_labeled_v_genes")
    make_v_profile.__globals__["make_v_region_summary"] = make_v_region_summary
    intersection = pd.DataFrame(
        {
            "v": ["TRAV1"] * 7,
            "region_id": ["100", "110", "101", "111", "010", "011", "001"],
            "region_name": ["region"] * 7,
            "clone_id": list("abcdefg"),
            "freq_vs_g1": [0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0],
        }
    )
    comparison = pd.DataFrame(
        {
            "v": ["TRAV1"],
            "frequency_in_full_g1": [0.2],
            "retained_cdr3_share": [0.75],
            "loss_share": [0.25],
        }
    )

    region_summary = make_v_region_summary(intersection, "freq_vs_g1")
    profile, _ = make_v_profile(comparison, intersection, "g1", "freq_vs_g1")

    assert set(region_summary["region_id"]) == {"100", "110", "101", "111"}
    assert profile.loc[0, "region_entropy"] == pytest.approx(1.0)
    assert profile.loc[0, "retention_score"] == pytest.approx(0.15)
    assert profile.loc[0, "loss_score"] == pytest.approx(0.05)

    candidates = pd.DataFrame(
        {
            "v": ["TRAV-low-share", "TRAV-middle", "TRAV-top"],
            "retention_score": [0.99, 0.20, 0.30],
            "retained_cdr3_share_pct": [100.0, 70.0, 80.0],
            "cdr3_in_full_g1": [10, 10, 10],
            "frequency_in_full_g1_pct": [0.5, 1.0, 2.0],
        }
    )
    selected_genes, _ = select_bubble_genes(candidates, "g1", top_n=2)
    assert selected_genes == ["TRAV-top", "TRAV-middle"]




def test_structural_downsampling_is_deterministic_and_excludes_low_depth_samples():
    notebook = json.loads((ROOT / "venn_original.ipynb").read_text(encoding="utf-8"))
    source = "".join(next(cell for cell in notebook["cells"] if cell["id"] == "venn-inputs")["source"])
    tree = ast.parse(source)
    selected = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name in {"_stable_sample_seed", "downsample_count_table"}
    ]
    namespace = {
        "np": __import__("numpy"),
        "pd": pd,
        "hashlib": hashlib,
        "DOWNSAMPLE_DEPTH": 15000,
        "DOWNSAMPLE_SEED": 1031,
    }
    exec(compile(ast.Module(body=selected, type_ignores=[]), "downsample_helpers", "exec"), namespace)
    downsample = namespace["downsample_count_table"]

    index = pd.Index([("AAA", "TRAV1"), ("BBB", "TRAV2")], tupleize_cols=False)
    table = pd.DataFrame(
        {
            "deep": [12000, 8000],
            "shallow": [10000, 4999],
        },
        index=index,
    )
    first, audit_first = downsample(table, depth=15000, seed=1031)
    second, audit_second = downsample(table, depth=15000, seed=1031)

    pd.testing.assert_frame_equal(first, second)
    pd.testing.assert_frame_equal(audit_first, audit_second)
    assert first.columns.tolist() == ["deep"]
    assert int(first["deep"].sum()) == 15000
    assert audit_first.set_index("sample_id").loc["shallow", "kept_for_structural_analysis"] == False

    with pytest.raises(ValueError, match="No sample reached"):
        downsample(
            table,
            depth=15000,
            seed=1031,
            eligible_samples={"shallow"},
        )

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
    (result_dir / "run_summary.json").write_text(
        json.dumps(
            {"edger_status": "edgeR quasi-likelihood model fitted successfully."}
        )
        + "\n",
        encoding="utf-8",
    )
    (result_dir / "trav_ranking.csv").write_text("v_gene\nTRAV1\n", encoding="utf-8")
    (result_dir / "distinctive_trav.csv").write_text(
        "v_gene\nTRAV1\n", encoding="utf-8"
    )
    manifest_rows = []
    for chain in ("TRA", "TRB"):
        chain_dir = figure_dir / chain
        chain_dir.mkdir(parents=True, exist_ok=True)
        png = chain_dir / f"{chain}_g1_v_segment_localization.png"
        png.write_bytes(b"png")
        manifest_rows.append(
            {
                "chain": chain,
                "png": str(png.relative_to(tmp_path)),
            }
        )
    pd.DataFrame(manifest_rows).to_csv(
        result_dir / "figure_manifest.csv", index=False
    )

    run_analysis.verify_first_approach_outputs(["cd4_thymus"])


def test_output_verification_rejects_an_unfitted_edger_model(tmp_path, monkeypatch):
    monkeypatch.setattr(run_analysis, "repo_root", lambda: tmp_path)
    stem = run_analysis.STRATUM_BY_KEY["cd4_thymus"]
    result_dir = tmp_path / "results" / "01_set_count" / "cd4_thymus"
    (tmp_path / "audit_runs").mkdir(parents=True)
    (tmp_path / "logs").mkdir(parents=True)
    result_dir.mkdir(parents=True)
    (tmp_path / "audit_runs" / f"{stem}.executed.ipynb").write_text(
        "{}", encoding="utf-8"
    )
    (tmp_path / "logs" / f"{stem}.log").write_text("complete\n", encoding="utf-8")
    (result_dir / "run_summary.json").write_text(
        json.dumps({"edger_status": "not fitted"}) + "\n", encoding="utf-8"
    )

    with pytest.raises(RuntimeError, match="edgeR did not complete"):
        run_analysis.verify_first_approach_outputs(["cd4_thymus"])
