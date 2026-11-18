import ast
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
import pytest

from scripts import run_analysis
from src.article_figures import ARTICLE_FIGURES, article_name, bubble_candidates


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
    assert "_render_article(fig, chain, panel)" in joined
    assert "plt.savefig = _article_savefig" in joined
    assert "plt.show = _article_show" in joined
    assert "sns.set_theme(style=\"white\"" in joined
    assert 'DOWNSAMPLE_DEPTH = int(os.environ.get("MICE_TCR_DOWNSAMPLE_UMI", "15000"))' in joined
    assert 'functionality="f"' in joined
    assert "by_umi=True" in joined
    assert "downsample_group_tables(" in joined
    assert 'globals()[f"raw_ct_{group}_aaV_tra"]' in joined
    assert 'globals()[f"ct_{group}_aaV_tra"]' in joined
    assert "prepare_counts(_tables_for_chain(chain, raw=True)" in joined
    assert "RAREFACTION_SEEDS" in joined
    assert "BALANCE_SEEDS" in joined
    assert "RAREFACTION_SEEDS = tuple(range(1201, 1301))" in joined
    assert 'REPRESENTATIVE_SEED_STRATA = {"all_combined", "cd4_combined", "cd8_combined"}' in joined
    assert "choose_representative_rarefaction_seed" in joined
    assert "top10_profile_distance_to_median" in joined
    assert "spike_in_used_for_selection" in joined
    assert "TRA_g1_top10_only_g1_stability.png" in joined
    assert 'profile["score_s"] = (' in joined
    assert 'profile["retention_score"] = profile["score_s"]' in joined
    assert "collect_pairwise_similarity_boxplot_table(results_v_js)" in joined
    assert 'output_dir = output_dir / "surplus"' in joined
    assert '_render_structural_figures(g1_profile, "TRA")' in joined
    assert '_render_structural_figures(g1_profile, "TRB")' in joined
    assert '_render_article(plot_subsets(profiles, chain), chain, "cd4_cd8")' in joined
    assert '_render_article(plot_depth_check(depth, chain, DOWNSAMPLE_DEPTH), chain, "depth")' in joined
    assert '_render_article(plot_mds(points, chain, DOWNSAMPLE_DEPTH), chain, "mds")' in joined
    assert "np.log2(len(target_region_order))" in joined
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
    for chain in ("tra", "trb"):
        for group in ("g1", "g2", "g3", "g4", "g5", "g6"):
            assert f"set(ct_{group}_aaV_{chain}.index)" in joined
            assert f"set(ct_{group}_aaV_{chain})" not in joined
    assert 'differential_chains = ("TRA", "TRB") if STRATUM == "all_combined"' in joined


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
    assert profile.loc[0, "score_s"] == pytest.approx(0.15)
    assert profile.loc[0, "retention_score"] == pytest.approx(0.15)
    assert profile.loc[0, "loss_score"] == pytest.approx(0.05)

    candidates = pd.DataFrame(
        {
            "v": ["TRAV-low-share", "TRAV-low-retention", "TRAV-middle", "TRAV-top"],
            "score_s": [0.99, 0.90, 0.20, 0.30],
            "retained_cdr3_share_pct": [100.0, 65.0, 66.0, 80.0],
            "cdr3_in_full_g1": [10, 10, 10, 10],
            "frequency_in_full_g1_pct": [0.5, 5.0, 0.6, 2.0],
            "region_entropy": [0.5, 0.4, 0.3, 0.2],
        }
    )
    displayed, labeled = bubble_candidates(candidates, limit=2)
    assert len(displayed) == 4
    assert labeled["v"].tolist() == ["TRAV-top", "TRAV-middle"]




def test_rarefaction_figures_use_count_distribution_and_boxplot_summary():
    notebook = json.loads((ROOT / "venn_original.ipynb").read_text(encoding="utf-8"))
    source = "".join(
        next(
            cell
            for cell in notebook["cells"]
            if cell["id"] == "venn-structural-robustness"
        )["source"]
    )

    assert 'ax.bar(' in source
    assert 'Reference top-10 V segments recovered' in source
    assert 'ax.boxplot(' in source
    assert '"mean | seed | IQR"' in source
    assert '100 reruns: IQR and median' in source
    assert 'fig.text(' not in source
    assert 'note=' not in source


def test_structural_score_and_only_g1_share_match_article_definitions():
    notebook = json.loads((ROOT / "venn_original.ipynb").read_text(encoding="utf-8"))
    source = "".join(
        next(cell for cell in notebook["cells"] if cell["id"] == "venn-inputs")["source"]
    )
    tree = ast.parse(source)
    selected = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name in {
            "_feature_v",
            "_structural_v_summary",
            "_ranked_structural_candidates",
        }
    ]
    namespace = {
        "np": __import__("numpy"),
        "pd": pd,
        "ast": ast,
    }
    exec(
        compile(ast.Module(body=selected, type_ignores=[]), "structural_summary", "exec"),
        namespace,
    )

    features = pd.Index(
        [
            ("AAA", "TRAV1"),
            ("AAB", "TRAV1"),
            ("AAC", "TRAV1"),
            ("BBB", "TRAV2"),
            ("BBC", "TRAV2"),
        ],
        tupleize_cols=False,
    )
    g1 = pd.DataFrame({"s1": [10, 10, 10, 10, 10]}, index=features)
    empty = pd.DataFrame({"s": []}, index=pd.Index([], dtype="object"))
    g5 = pd.DataFrame(
        {"s5": [4]},
        index=pd.Index([("AAA", "TRAV1")], tupleize_cols=False),
    )
    tables = {
        "g1": g1,
        "g2": empty,
        "g4": empty,
        "g5": g5,
        "g6": empty,
    }

    summary = namespace["_structural_v_summary"](tables).set_index("v")
    assert summary.loc["TRAV1", "score_s"] == pytest.approx(2 / 5)
    assert summary.loc["TRAV1", "retention"] == pytest.approx(2 / 3)
    assert summary.loc["TRAV1", "only_g1_share"] == pytest.approx(2 / 3)
    assert summary.loc["TRAV2", "score_s"] == pytest.approx(2 / 5)
    ranked = namespace["_ranked_structural_candidates"](summary.reset_index())
    assert ranked["v"].tolist()[:2] == ["TRAV1", "TRAV2"]


def test_figure_chain_detection_ignores_repository_name_and_routes_trb():
    notebook = json.loads((ROOT / "venn_original.ipynb").read_text(encoding="utf-8"))
    source = "".join(
        next(cell for cell in notebook["cells"] if cell["id"] == "venn-figure-capture")["source"]
    )
    tree = ast.parse(source)
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "_figure_chain"
    )
    namespace = {
        "Path": Path,
        "re": re,
        "ACTIVE_CHAIN": "TRA and TRB",
        "_figure_title": lambda value: value,
    }
    exec(
        compile(ast.Module(body=[function], type_ignores=[]), "figure_chain", "exec"),
        namespace,
    )
    detect = namespace["_figure_chain"]

    assert detect(
        "TRB g1: 15,000-UMI rarefaction stability",
        "/projects/mice_transplant_2025/figures/TRB_g1_rarefaction_stability.png",
    ) == "TRB"
    assert detect(
        "TRA g1: 15,000-UMI rarefaction stability",
        "/projects/mice_transplant_2025/figures/TRA_g1_rarefaction_stability.png",
    ) == "TRA"


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


@pytest.mark.parametrize("stratum", ["cd4_thymus", "all_combined"])
def test_output_verification_accepts_a_complete_single_stratum(tmp_path, monkeypatch, stratum):
    monkeypatch.setattr(run_analysis, "repo_root", lambda: tmp_path)
    stem = run_analysis.STRATUM_BY_KEY[stratum]
    result_dir = tmp_path / "results" / "01_set_count" / stratum
    figure_dir = tmp_path / "figures" / "01_set_count" / stratum
    audit_dir = tmp_path / "audit_runs"
    log_dir = tmp_path / "logs"
    for directory in (result_dir, figure_dir, audit_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    (audit_dir / f"{stem}.executed.ipynb").write_text("{}", encoding="utf-8")
    (log_dir / f"{stem}.log").write_text("complete\n", encoding="utf-8")
    (result_dir / "run_summary.json").write_text(
        json.dumps(
            {"edger_status": "edgeR quasi-likelihood model fitted successfully.",
             "edger_status_by_chain": {"TRB": "edgeR quasi-likelihood model fitted successfully."}}
        )
        + "\n",
        encoding="utf-8",
    )
    (result_dir / "trav_ranking.csv").write_text("v_gene\nTRAV1\n", encoding="utf-8")
    (result_dir / "distinctive_trav.csv").write_text(
        "v_gene\nTRAV1\n", encoding="utf-8"
    )
    manifest_rows = []
    panels = ["localization"]
    if stratum == "all_combined":
        panels = [panel for panel in ARTICLE_FIGURES if panel != "spike_in"]
        for name in (
            "representative_downsampling_seed.json", "tra_rarefaction_seed_summary.csv",
            "tra_rarefaction_v_profiles.csv", "tra_top10_only_g1_by_seed.csv",
            "tra_top10_only_g1_summary.csv", "trbv_ranking.csv", "distinctive_trbv.csv",
            "tra_depth_check.csv", "trb_depth_check.csv", "tra_mds_coordinates.csv", "trb_mds_coordinates.csv",
        ):
            (result_dir / name).write_text("test\n", encoding="utf-8")
    for chain in ("TRA", "TRB"):
        names = [figure_dir / chain / article_name(chain, panel) for panel in panels]
        if stratum == "all_combined":
            extra = figure_dir / chain / "surplus"
            names += [extra / f"{chain}_g1_{stem}.png" for stem in (
                "rarefaction_stability", "balanced_mouse_stability",
            )]
            if chain == "TRA":
                names.append(extra / "TRA_g1_top10_only_g1_stability.png")
        for png in names:
            png.parent.mkdir(parents=True, exist_ok=True)
            png.write_bytes(b"png")
            manifest_rows.append({"chain": chain, "png": str(png.relative_to(tmp_path))})
    pd.DataFrame(manifest_rows).to_csv(
        result_dir / "figure_manifest.csv", index=False
    )

    run_analysis.verify_first_approach_outputs([stratum])


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
