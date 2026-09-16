#!/usr/bin/env python3
"""Prepare and submit bounded Aldan-3 baseline -> dose array -> summary jobs.

Preparation and diagnostics never submit jobs. Submission must be requested explicitly.
Each job has a fresh Python process/kernel and one case-specific output directory.
"""
from __future__ import annotations

import argparse
import contextlib
import getpass
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DENSE_FRACTIONS = "0.0000001,0.0000002,0.0000005,0.000001,0.000002,0.000005,0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01"
INPUT_KEYS = ("MICE_TCR_DATA_DIR", "MICE_TCR_METADATA_CSV",
              "MICE_TCR_CLONOSET_INDEX", "MICE_TCR_WORKING_DIR")


def thread_environment(cpus):
    return {**{name: str(cpus) for name in (
        "ARTICLE_N_JOBS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "POLARS_MAX_THREADS")}, "OMP_DYNAMIC": "FALSE", "MKL_DYNAMIC": "FALSE",
        "OMP_MAX_ACTIVE_LEVELS": "1", "MPLBACKEND": "Agg"}


def validate_resources(args):
    if min(args.cpus, args.baseline_cpus, args.max_parallel, args.cpu_budget) < 1:
        raise ValueError("CPU counts, concurrency and CPU budget must be positive.")
    if max(args.baseline_cpus, args.cpus * args.max_parallel) > args.cpu_budget:
        raise ValueError("Baseline CPUs or concurrent dose CPUs exceed --cpu-budget.")
    if not re.fullmatch(r"[1-9][0-9]*[MGT]", args.mem):
        raise ValueError("Use an explicit positive memory reservation, e.g. --mem 64G.")
    for value in (args.partition, args.constraint, args.account):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
            raise ValueError("Partition, constraint and account must be simple Slurm names.")
    match = re.fullmatch(r"(\d+):(\d{2}):(\d{2})", args.time)
    if not match or int(match[2]) > 59 or int(match[3]) > 59:
        raise ValueError("Time must use HH:MM:SS.")
    seconds = int(match[1]) * 3600 + int(match[2]) * 60 + int(match[3])
    limits = {"short": 2 * 3600, "medium": 16 * 3600}
    if seconds < 1 or seconds > limits.get(args.partition, math.inf):
        raise ValueError("Requested time exceeds the documented partition limit.")


def prepare(args):
    from scripts.run_analysis import analysis_fingerprint
    validate_resources(args)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_id):
        raise ValueError("Invalid run ID.")
    fractions = sorted(set(float(x) for x in args.fractions.split(",")))
    if not fractions or any(not math.isfinite(x) or not 0 < x < 1 for x in fractions):
        raise ValueError("Fractions must be finite values between zero and one.")
    if args.clones < 1:
        raise ValueError("At least one clone is required.")
    for directory in (ROOT / "results/01_spike_in", ROOT / "figures/01_spike_in",
                      ROOT / "audit_runs/spike_in", ROOT / "logs/spike_in"):
        if (directory / args.run_id).exists():
            raise FileExistsError("Analysis run already exists; use a new --run-id.")
    plan_dir = ROOT / "logs/slurm" / args.run_id
    plan_dir.mkdir(parents=True, exist_ok=False)
    plan_path = plan_dir / "plan.json"
    runner = [str(ROOT / "scripts/run_analysis.py"), "--approach", "1", "--mode", "spike-in",
              "--strata", "all_combined", "--spike-run-id", args.run_id,
              "--spike-clones", str(args.clones), "--spike-seed", str(args.seed),
              "--spike-fractions", ",".join(str(x) for x in fractions)]
    if args.v_gene:
        runner += ["--spike-v-gene", args.v_gene]
    environment = {key: os.environ.get(key) for key in INPUT_KEYS}
    for argument, key in zip((args.data_dir, args.metadata_csv, args.clonoset_index,
                              args.working_dir), INPUT_KEYS):
        if argument is not None:
            environment[key] = str(argument.expanduser().resolve())
    plan = {"schema": 1, "root": str(ROOT), "run_id": args.run_id,
            "python": sys.executable, "runner": runner, "fractions": fractions,
            "environment": environment, "source_sha256": analysis_fingerprint(ROOT),
            "resources": {key: getattr(args, key) for key in (
                "cpus", "baseline_cpus", "max_parallel", "cpu_budget", "mem",
                "time", "partition", "constraint", "account")}}
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    for stage in ("baseline", "dose", "finalize"):
        command = [sys.executable, str(ROOT / "scripts/slurm_spike.py"), "worker",
                   str(plan_path), "--stage", stage]
        script = plan_dir / f"{stage}.sh"
        script.write_text("#!/bin/bash\nset -euo pipefail\n"
            + f"cd {shlex.quote(str(ROOT))}\n"
            + 'exec srun --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" '
              '--cpu-bind=cores ' + shlex.join(command) + "\n")
        script.chmod(0o750)
    print(f"Prepared only; no jobs submitted. Plan: {plan_path}")
    print(f"Baseline: {args.baseline_cpus} CPUs; doses: up to {args.max_parallel} x {args.cpus} CPUs.")
    print(f"Memory reservation: {args.mem} per analysis job; choose using measured MaxRSS.")
    print(f"Submit pilot: {shlex.join([sys.executable, str(ROOT / 'scripts/slurm_spike.py'), 'submit', str(plan_path), '--stage', 'baseline'])}")
    return plan_path


def load_plan(path):
    from scripts.run_analysis import analysis_fingerprint
    plan = json.loads(path.read_text())
    if plan["root"] != str(ROOT) or plan["source_sha256"] != analysis_fingerprint(ROOT):
        raise ValueError("Repository path or source differs from the prepared plan.")
    if not Path(plan["python"]).is_file():
        raise FileNotFoundError("The prepared Python environment is unavailable.")
    validate_resources(argparse.Namespace(**plan["resources"]))
    return plan


def submission_command(plan, path, stage, dependency=None):
    resources = plan["resources"]
    cpus = {"baseline": resources["baseline_cpus"], "dose": resources["cpus"], "finalize": 1}[stage]
    command = ["sbatch", "--parsable", "--nodes=1", "--ntasks=1", "--no-requeue",
               f"--cpus-per-task={cpus}", f"--mem={resources['mem']}",
               f"--time={resources['time']}", f"--partition={resources['partition']}",
               f"--constraint={resources['constraint']}",
               f"--job-name=spike_{stage}_{plan['run_id']}",
               f"--output={path.parent / (stage + '_%A_%a.log')}",
               f"--chdir={plan['root']}", "--export=ALL"]
    if resources["account"]:
        command += [f"--account={resources['account']}"]
    if stage == "dose":
        command += [f"--array=1-{len(plan['fractions'])}%{resources['max_parallel']}"]
    if dependency:
        command += [f"--dependency=afterok:{dependency}", "--kill-on-invalid-dep=yes"]
    command += [str(path.parent / f"{stage}.sh")]
    return command


@contextlib.contextmanager
def submission_lock(path):
    lock = path.parent / "submission.lock"
    with lock.open("x"):
        try:
            yield
        finally:
            lock.unlink()


def submit(path, stage):
    plan = load_plan(path)
    jobs_path = path.parent / "jobs.json"
    with submission_lock(path):
        jobs = json.loads(jobs_path.read_text()) if jobs_path.exists() else {}
        requested = {"baseline": ["baseline"], "remaining": ["dose", "finalize"],
                     "all": ["baseline", "dose", "finalize"]}[stage]
        if stage == "remaining" and "baseline" not in jobs:
            raise ValueError("Submit the baseline first, or use --stage all.")
        for name in requested:
            if name in jobs:
                print(f"Already submitted {name}: {jobs[name]}")
                continue
            predecessor = {"dose": "baseline", "finalize": "dose"}.get(name)
            command = submission_command(plan, path, name, jobs.get(predecessor))
            print(shlex.join(command), flush=True)
            result = subprocess.run(command, text=True, capture_output=True, check=True)
            job_id = result.stdout.strip().split(";")[0]
            if not job_id.isdigit():
                raise RuntimeError(f"Unexpected sbatch response; inspect queue before retrying: {result.stdout!r}")
            jobs[name] = job_id
            temporary = jobs_path.with_suffix(".tmp")
            temporary.write_text(json.dumps(jobs, indent=2) + "\n")
            temporary.replace(jobs_path)
            print(f"Submitted {name}: {job_id}", flush=True)
    print(f"Job IDs: {jobs_path}")


def worker(path, stage):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Workers require a Slurm allocation; use submit from the login node.")
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    if cpus < 1:
        raise RuntimeError("Missing SLURM_CPUS_PER_TASK.")
    # set limits before importing the numerical runner or launching a kernel
    os.environ.update(thread_environment(cpus))
    plan = load_plan(path)
    expected = {"baseline": plan["resources"]["baseline_cpus"],
                "dose": plan["resources"]["cpus"], "finalize": 1}[stage]
    if cpus != expected:
        raise ValueError("The worker allocation does not match its prepared resource plan.")
    command = [plan["python"], *plan["runner"], "--spike-stage", stage]
    if stage == "dose":
        index = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
        if not 1 <= index <= len(plan["fractions"]):
            raise ValueError("Invalid dose array index.")
        command += ["--spike-dose-index", str(index)]
    for key, value in plan["environment"].items():
        if value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = value
    print(f"Worker stage={stage}, CPUs={cpus}, Python={plan['python']}", flush=True)
    os.execve(plan["python"], command, dict(os.environ))


def diagnose():
    """A few read-only scheduler calls, with no polling or job submission."""
    commands = [
        ["sinfo", "-o", "%P %l %D %c %m %f"],
        ["scontrol", "show", "partition", "short"],
        ["scontrol", "show", "partition", "medium"],
        ["squeue", "-u", getpass.getuser(), "-o", "%.18i %.12P %.24j %.8T %.10M %.6C %.10m %R"],
    ]
    for command in commands:
        print(shlex.join(command), flush=True)
        try:
            subprocess.run(command, check=True, timeout=20)
        except (OSError, subprocess.SubprocessError) as error:
            print(f"Diagnostic unavailable: {error}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare", help="Write a plan/scripts without submitting anything.")
    p.add_argument("--run-id", required=True)
    p.add_argument("--mem", required=True, help="Memory per job, selected from measured MaxRSS; e.g. 64G.")
    p.add_argument("--cpus", type=int, default=8, help="CPUs per dose notebook.")
    p.add_argument("--baseline-cpus", type=int, default=24)
    p.add_argument("--max-parallel", type=int, default=3)
    p.add_argument("--cpu-budget", type=int, default=24)
    p.add_argument("--partition", default="short")
    p.add_argument("--constraint", default="hpc")
    p.add_argument("--time", default="02:00:00")
    p.add_argument("--account")
    p.add_argument("--fractions", default=DENSE_FRACTIONS, help="Unit fractions, not percentages.")
    p.add_argument("--clones", type=int, default=10)
    p.add_argument("--seed", type=int, default=1031)
    p.add_argument("--v-gene")
    for name in ("data-dir", "metadata-csv", "clonoset-index", "working-dir"):
        p.add_argument(f"--{name}", type=Path)
    s = sub.add_parser("submit", help="Submit prepared jobs to Slurm.")
    s.add_argument("plan", type=Path)
    s.add_argument("--stage", choices=("baseline", "remaining", "all"), default="baseline")
    w = sub.add_parser("worker", help=argparse.SUPPRESS)
    w.add_argument("plan", type=Path)
    w.add_argument("--stage", choices=("baseline", "dose", "finalize"), required=True)
    sub.add_parser("diagnose", help="Read-only partition/resource/queue information.")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args)
    elif args.action == "submit":
        submit(args.plan.resolve(), args.stage)
    elif args.action == "worker":
        worker(args.plan.resolve(), args.stage)
    else:
        diagnose()


if __name__ == "__main__":
    main()
