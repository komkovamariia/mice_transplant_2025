import argparse
import json
import os

import numpy as np
import pandas as pd
import pytest
from scipy.stats import fisher_exact

from scripts import run_analysis, slurm_spike
from src import runtime, spike_in
from src.count_statistics import pooled_fisher_rows


@pytest.fixture
def cpu_environment(monkeypatch):
    for key in ("ARTICLE_N_JOBS", "SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_ID"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(runtime, "_DEFAULT_WORKERS", None)
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: set(range(48)))
    monkeypatch.setattr(os, "cpu_count", lambda: 64)


def test_explicit_worker_requests_cannot_escape_slurm_or_affinity(cpu_environment, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "101")
    monkeypatch.setenv("ARTICLE_N_JOBS", "112")
    assert runtime.available_cpus() == 8
    assert runtime.configure_runtime(112, verbose=False) == 8
    assert runtime._resolve_jobs(112) == 8
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: {0, 1})
    assert runtime._resolve_jobs(112) == 2


def test_login_defaults_to_one_and_missing_task_allocation_is_conservative(cpu_environment, monkeypatch):
    assert runtime.available_cpus() == 1
    monkeypatch.setenv("SLURM_JOB_ID", "101")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "48")
    monkeypatch.setenv("ARTICLE_N_JOBS", "24")
    assert runtime.available_cpus() == 1
    with pytest.raises(ValueError):
        runtime.configure_runtime(0, verbose=False)


def test_parallel_fisher_matches_historical_exact_tests_and_input_order(monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    rng = np.random.default_rng(1031)
    g1 = pd.Series(np.r_[0, 0, 100, rng.integers(0, 50, 61)], dtype=float)
    ctrl = pd.Series(np.r_[0, 100, 0, rng.integers(0, 50, 61)], dtype=float)
    ids = [f"aaV_{i}" for i in range(len(g1))]
    expected = []
    for feature, a, b in zip(ids, g1, ctrl):
        odds, p = fisher_exact([[a, b], [g1.sum() - a, ctrl.sum() - b]], alternative="two-sided")
        expected.append([feature, a, b, np.log2((a + 1) / (b + 1)), odds, p])
    sequential = pd.DataFrame(pooled_fisher_rows(ids, g1, ctrl, n_jobs=1, batch_size=7))
    parallel = pd.DataFrame(pooled_fisher_rows(ids, g1, ctrl, n_jobs=2, batch_size=7))
    pd.testing.assert_frame_equal(parallel, sequential, check_exact=True)
    pd.testing.assert_frame_equal(parallel, pd.DataFrame(expected, columns=parallel.columns), check_exact=True)


def resources(**overrides):
    values = dict(cpus=8, baseline_cpus=24, max_parallel=3, cpu_budget=24,
                  mem="64G", time="02:00:00", partition="short", constraint="hpc", account=None)
    return argparse.Namespace(**(values | overrides))


def test_resource_budget_rejects_nested_total_above_24_and_unbounded_memory():
    slurm_spike.validate_resources(resources())
    for overrides in ({"max_parallel": 4}, {"baseline_cpus": 25}, {"mem": "0"},
                      {"time": "03:00:00"}, {"cpus": 0}):
        with pytest.raises(ValueError):
            slurm_spike.validate_resources(resources(**overrides))


def test_submission_dependencies_throttle_and_repeat_protection(tmp_path, monkeypatch):
    path = tmp_path / "plan.json"
    plan = {"resources": vars(resources()), "root": str(tmp_path), "run_id": "test",
            "fractions": list(range(16))}
    monkeypatch.setattr(slurm_spike, "load_plan", lambda _: plan)
    submitted = []

    def fake_submit(command, **kwargs):
        submitted.append(command)
        return argparse.Namespace(stdout=str(100 + len(submitted)))

    monkeypatch.setattr(slurm_spike.subprocess, "run", fake_submit)
    slurm_spike.submit(path, "all")
    slurm_spike.submit(path, "all")
    assert len(submitted) == 3
    assert "--cpus-per-task=24" in submitted[0]
    assert "--cpus-per-task=8" in submitted[1]
    assert "--array=1-16%3" in submitted[1]
    assert "--dependency=afterok:101" in submitted[1]
    assert "--dependency=afterok:102" in submitted[2]
    assert "--cpus-per-task=1" in submitted[2]


def test_worker_refuses_login_node_execution(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(RuntimeError, match="Slurm allocation"):
        slurm_spike.worker(tmp_path / "plan.json", "baseline")


@pytest.fixture
def staged_run(tmp_path, monkeypatch):
    monkeypatch.setattr(run_analysis, "repo_root", lambda: tmp_path)
    monkeypatch.setattr(run_analysis, "analysis_fingerprint", lambda _: "fixed-source")
    monkeypatch.setattr(run_analysis, "verify_first_approach_outputs", lambda *a, **kw: None)
    monkeypatch.setattr(spike_in, "load_snapshot_tables", lambda *a: {})
    monkeypatch.setattr(spike_in, "select_targets", lambda *a, **kw: {"clones": ["one", "two"]})
    monkeypatch.setattr(spike_in, "write_baseline_metrics", lambda *a: None)
    calls = []

    def execute(path, stem, *, environment, resume):
        calls.append(environment.copy())
        case = environment["MICE_TCR_SPIKE_CASE"]
        result = tmp_path / "results/01_spike_in/test" / case
        result.mkdir(parents=True)
        (result / "tra_v_retention.csv").write_text("gene\nTRAV1\n")
        (result / "trav_ranking.csv").write_text("gene\nTRAV1\n")

    monkeypatch.setattr(run_analysis, "execute_notebook", execute)
    summaries = []

    def summarize(result, figures, cases):
        summaries.append(cases)
        figures.mkdir(parents=True)
        (figures / "summary.png").write_bytes(b"test")

    monkeypatch.setattr(spike_in, "summarize_experiment", summarize)
    arguments = argparse.Namespace(approach="1", strata="all_combined", resume=False,
        spike_stage="baseline", spike_run_id="test", spike_fractions="0.0000001,0.01",
        spike_clones=2, spike_max_g1_fraction=1e-5, spike_seed=1031, spike_v_gene=None)
    logger = argparse.Namespace(write=lambda message: None)
    return tmp_path, arguments, logger, calls, summaries


def test_independent_doses_keep_identity_and_baseline_runs_once(staged_run):
    root, args, logger, calls, summaries = staged_run
    run_analysis.run_spike_experiment(args, {}, logger)
    args.spike_stage = "dose"
    # Out-of-order completion is normal for Slurm arrays.
    for index in (2, 1):
        args.spike_dose_index = index
        run_analysis.run_spike_experiment(args, {}, logger)
    with pytest.raises(FileExistsError):
        run_analysis.run_spike_experiment(args, {}, logger)
    args.spike_stage = "finalize"
    run_analysis.run_spike_experiment(args, {}, logger)
    assert [(x["MICE_TCR_SPIKE_CASE"], float(x["MICE_TCR_SPIKE_FRACTION"])) for x in calls] == [
        ("baseline", 0), ("dose_02", 0.01), ("dose_01", 1e-7)]
    assert all(x["MICE_TCR_STRATUM"] == "all_combined" for x in calls)
    assert summaries == [["baseline", "dose_01", "dose_02"]]
    assert (root / "figures/01_spike_in/test.zip").exists()


def test_doses_require_completed_unchanged_baseline(staged_run):
    root, args, logger, calls, _ = staged_run
    args.spike_stage = "dose"
    args.spike_dose_index = 1
    with pytest.raises(FileNotFoundError):
        run_analysis.run_spike_experiment(args, {}, logger)
    args.spike_stage = "baseline"
    run_analysis.run_spike_experiment(args, {}, logger)
    args.spike_stage = "dose"
    args.spike_seed += 1
    with pytest.raises(ValueError, match="settings differ"):
        run_analysis.run_spike_experiment(args, {}, logger)
    args.spike_seed -= 1
    (root / "results/01_spike_in/test/selection.json").write_text("{}")
    with pytest.raises(ValueError, match="selection changed"):
        run_analysis.run_spike_experiment(args, {}, logger)
    assert len(calls) == 1


def test_partial_array_cannot_publish_summary(staged_run):
    root, args, logger, _, summaries = staged_run
    run_analysis.run_spike_experiment(args, {}, logger)
    args.spike_stage = "finalize"
    with pytest.raises(FileNotFoundError):
        run_analysis.run_spike_experiment(args, {}, logger)
    assert not summaries
    assert not (root / "figures/01_spike_in/test.zip").exists()
