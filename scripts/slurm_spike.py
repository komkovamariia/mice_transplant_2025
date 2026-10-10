#!/usr/bin/env python3
"""Prepare and submit the bounded Slurm spike-in run."""

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

DENSE_FRACTIONS = (
    "0.0000001,0.0000002,0.0000005,0.000001,0.000002,0.000005,"
    "0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,"
    "0.005,0.01"
)
INPUT_KEYS = (
    "MICE_TCR_DATA_DIR",
    "MICE_TCR_METADATA_CSV",
    "MICE_TCR_CLONOSET_INDEX",
    "MICE_TCR_WORKING_DIR",
)


def thread_environment(cpus: int) -> dict[str, str]:
    thread_keys = (
        "ARTICLE_N_JOBS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "POLARS_MAX_THREADS",
    )
    environment = {name: str(cpus) for name in thread_keys}
    environment.update(
        OMP_DYNAMIC="FALSE",
        MKL_DYNAMIC="FALSE",
        OMP_MAX_ACTIVE_LEVELS="1",
        MPLBACKEND="Agg",
    )
    return environment


def validate_resources(args: argparse.Namespace) -> None:
    counts = (args.cpus, args.baseline_cpus, args.max_parallel, args.cpu_budget)
    if min(counts) < 1:
        raise ValueError("CPU counts, concurrency and CPU budget must be positive.")

    if max(args.baseline_cpus, args.cpus * args.max_parallel) > args.cpu_budget:
        raise ValueError("Baseline CPUs or concurrent dose CPUs exceed --cpu-budget.")

    if not re.fullmatch(r"[1-9][0-9]*[MGT]", args.mem):
        raise ValueError("Use an explicit positive memory reservation, e.g. --mem 64G.")

    for value in (args.partition, args.constraint, args.account):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
            raise ValueError(
                "Partition, constraint and account must be simple Slurm names."
            )

    match = re.fullmatch(r"(\d+):(\d{2}):(\d{2})", args.time)
    if not match or int(match[2]) > 59 or int(match[3]) > 59:
        raise ValueError("Time must use HH:MM:SS.")

    seconds = int(match[1]) * 3600 + int(match[2]) * 60 + int(match[3])
    limits = {"short": 2 * 3600, "medium": 16 * 3600}
    if seconds < 1 or seconds > limits.get(args.partition, math.inf):
        raise ValueError("Requested time exceeds the documented partition limit.")


def prepare(args: argparse.Namespace) -> Path:
    from scripts.run_analysis import analysis_fingerprint

    validate_resources(args)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_id):
        raise ValueError("Invalid run ID.")

    fractions = sorted(set(float(value) for value in args.fractions.split(",")))
    if not fractions or any(
        not math.isfinite(value) or not 0 < value < 1 for value in fractions
    ):
        raise ValueError("Fractions must be finite values between zero and one.")
    if args.clones < 1:
        raise ValueError("At least one clone is required.")

    output_roots = (
        ROOT / "results/01_spike_in",
        ROOT / "figures/01_spike_in",
        ROOT / "audit_runs/spike_in",
        ROOT / "logs/spike_in",
    )
    if any((directory / args.run_id).exists() for directory in output_roots):
        raise FileExistsError("Analysis run already exists; use a new --run-id.")

    plan_dir = ROOT / "logs/slurm" / args.run_id
    plan_dir.mkdir(parents=True, exist_ok=False)
    plan_path = plan_dir / "plan.json"

    runner = [
        str(ROOT / "scripts/run_analysis.py"),
        "--approach", "1",
        "--mode", "spike-in",
        "--strata", "all_combined",
        "--spike-run-id", args.run_id,
        "--spike-clones", str(args.clones),
        "--spike-seed", str(args.seed),
        "--spike-fractions", ",".join(str(value) for value in fractions),
    ]
    if args.v_gene:
        runner += ["--spike-v-gene", args.v_gene]

    environment = {key: os.environ.get(key) for key in INPUT_KEYS}
    path_args = (
        args.data_dir,
        args.metadata_csv,
        args.clonoset_index,
        args.working_dir,
    )
    for path_arg, key in zip(path_args, INPUT_KEYS):
        if path_arg is not None:
            environment[key] = str(path_arg.expanduser().resolve())

    resource_keys = (
        "cpus",
        "baseline_cpus",
        "max_parallel",
        "cpu_budget",
        "mem",
        "time",
        "partition",
        "constraint",
        "account",
    )
    plan = {
        "schema": 1,
        "root": str(ROOT),
        "run_id": args.run_id,
        "python": sys.executable,
        "runner": runner,
        "fractions": fractions,
        "environment": environment,
        "source_sha256": analysis_fingerprint(ROOT),
        "resources": {key: getattr(args, key) for key in resource_keys},
    }
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")

    for stage in ("baseline", "dose", "finalize"):
        command = [
            sys.executable,
            str(ROOT / "scripts/slurm_spike.py"),
            "worker",
            str(plan_path),
            "--stage",
            stage,
        ]
        script = plan_dir / f"{stage}.sh"
        script.write_text(
            "#!/bin/bash\n"
            "set -euo pipefail\n"
            f"cd {shlex.quote(str(ROOT))}\n"
            'exec srun --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" '
            "--cpu-bind=cores "
            + shlex.join(command)
            + "\n",
            encoding="utf-8",
        )
        script.chmod(0o750)

    print(f"Prepared only; no jobs submitted. Plan: {plan_path}")
    print(
        f"Baseline: {args.baseline_cpus} CPUs; "
        f"doses: up to {args.max_parallel} x {args.cpus} CPUs."
    )
    print(f"Memory reservation: {args.mem} per job; choose from measured MaxRSS.")
    submit = [
        sys.executable,
        str(ROOT / "scripts/slurm_spike.py"),
        "submit",
        str(plan_path),
        "--stage",
        "baseline",
    ]
    print(f"Submit pilot: {shlex.join(submit)}")
    return plan_path


def load_plan(path: Path) -> dict:
    from scripts.run_analysis import analysis_fingerprint

    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan["root"] != str(ROOT) or plan["source_sha256"] != analysis_fingerprint(ROOT):
        raise ValueError("Repository path or source differs from the prepared plan.")
    if not Path(plan["python"]).is_file():
        raise FileNotFoundError("The prepared Python environment is unavailable.")

    validate_resources(argparse.Namespace(**plan["resources"]))
    return plan


def submission_command(
    plan: dict,
    path: Path,
    stage: str,
    dependency: str | None = None,
) -> list[str]:
    resources = plan["resources"]
    cpus = {
        "baseline": resources["baseline_cpus"],
        "dose": resources["cpus"],
        "finalize": 1,
    }[stage]
    command = [
        "sbatch",
        "--parsable",
        "--nodes=1",
        "--ntasks=1",
        "--no-requeue",
        f"--cpus-per-task={cpus}",
        f"--mem={resources['mem']}",
        f"--time={resources['time']}",
        f"--partition={resources['partition']}",
        f"--constraint={resources['constraint']}",
        f"--job-name=spike_{stage}_{plan['run_id']}",
        f"--output={path.parent / (stage + '_%A_%a.log')}",
        f"--chdir={plan['root']}",
        "--export=ALL",
    ]
    if resources["account"]:
        command.append(f"--account={resources['account']}")
    if stage == "dose":
        command.append(
            f"--array=1-{len(plan['fractions'])}%{resources['max_parallel']}"
        )
    if dependency:
        command += [
            f"--dependency=afterok:{dependency}",
            "--kill-on-invalid-dep=yes",
        ]
    command.append(str(path.parent / f"{stage}.sh"))
    return command


@contextlib.contextmanager
def submission_lock(path: Path):
    lock_path = path.parent / "submission.lock"
    with lock_path.open("x"):
        try:
            yield
        finally:
            lock_path.unlink()


def submit(path: Path, stage: str) -> None:
    plan = load_plan(path)
    jobs_path = path.parent / "jobs.json"

    with submission_lock(path):
        jobs = (
            json.loads(jobs_path.read_text(encoding="utf-8"))
            if jobs_path.exists()
            else {}
        )
        requested = {
            "baseline": ["baseline"],
            "remaining": ["dose", "finalize"],
            "all": ["baseline", "dose", "finalize"],
        }[stage]
        if stage == "remaining" and "baseline" not in jobs:
            raise ValueError("Submit the baseline first, or use --stage all.")

        for job_name in requested:
            if job_name in jobs:
                print(f"Already submitted {job_name}: {jobs[job_name]}")
                continue

            predecessor = {"dose": "baseline", "finalize": "dose"}.get(job_name)
            command = submission_command(
                plan,
                path,
                job_name,
                jobs.get(predecessor),
            )
            print(shlex.join(command), flush=True)
            result = subprocess.run(
                command,
                text=True,
                capture_output=True,
                check=True,
            )
            job_id = result.stdout.strip().split(";")[0]
            if not job_id.isdigit():
                raise RuntimeError(
                    "Unexpected sbatch response; inspect the queue before retrying: "
                    f"{result.stdout!r}"
                )

            jobs[job_name] = job_id
            temporary = jobs_path.with_suffix(".tmp")
            temporary.write_text(
                json.dumps(jobs, indent=2) + "\n",
                encoding="utf-8",
            )
            temporary.replace(jobs_path)
            print(f"Submitted {job_name}: {job_id}", flush=True)

    print(f"Job IDs: {jobs_path}")


def worker(path: Path, stage: str) -> None:
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError(
            "Workers require a Slurm allocation; use submit from the login node."
        )

    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    if cpus < 1:
        raise RuntimeError("Missing SLURM_CPUS_PER_TASK.")

    # set thread limits before numerical libraries or kernels start
    os.environ.update(thread_environment(cpus))
    plan = load_plan(path)

    expected_cpus = {
        "baseline": plan["resources"]["baseline_cpus"],
        "dose": plan["resources"]["cpus"],
        "finalize": 1,
    }[stage]
    if cpus != expected_cpus:
        raise ValueError(
            "The worker allocation does not match its prepared resource plan."
        )

    command = [plan["python"], *plan["runner"], "--spike-stage", stage]
    if stage == "dose":
        dose_index = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
        if not 1 <= dose_index <= len(plan["fractions"]):
            raise ValueError("Invalid dose array index.")
        command += ["--spike-dose-index", str(dose_index)]

    for key, value in plan["environment"].items():
        if value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = value

    print(
        f"Worker stage={stage}, CPUs={cpus}, Python={plan['python']}",
        flush=True,
    )
    os.execve(plan["python"], command, dict(os.environ))


def diagnose() -> None:
    """Run the read-only scheduler checks used before submission."""
    commands = [
        ["sinfo", "-o", "%P %l %D %c %m %f"],
        ["scontrol", "show", "partition", "short"],
        ["scontrol", "show", "partition", "medium"],
        [
            "squeue",
            "-u",
            getpass.getuser(),
            "-o",
            "%.18i %.12P %.24j %.8T %.10M %.6C %.10m %R",
        ],
    ]
    for command in commands:
        print(shlex.join(command), flush=True)
        try:
            subprocess.run(command, check=True, timeout=20)
        except (OSError, subprocess.SubprocessError) as error:
            print(f"Diagnostic unavailable: {error}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare",
        help="Write a plan and stage scripts without submitting jobs.",
    )
    prepare_parser.add_argument("--run-id", required=True)
    prepare_parser.add_argument(
        "--mem",
        required=True,
        help="Memory per job, chosen from measured MaxRSS; e.g. 64G.",
    )
    prepare_parser.add_argument("--cpus", type=int, default=8)
    prepare_parser.add_argument("--baseline-cpus", type=int, default=24)
    prepare_parser.add_argument("--max-parallel", type=int, default=3)
    prepare_parser.add_argument("--cpu-budget", type=int, default=24)
    prepare_parser.add_argument("--partition", default="short")
    prepare_parser.add_argument("--constraint", default="hpc")
    prepare_parser.add_argument("--time", default="02:00:00")
    prepare_parser.add_argument("--account")
    prepare_parser.add_argument(
        "--fractions",
        default=DENSE_FRACTIONS,
        help="Unit fractions, not percentages.",
    )
    prepare_parser.add_argument("--clones", type=int, default=10)
    prepare_parser.add_argument("--seed", type=int, default=1031)
    prepare_parser.add_argument("--v-gene")
    for name in ("data-dir", "metadata-csv", "clonoset-index", "working-dir"):
        prepare_parser.add_argument(f"--{name}", type=Path)

    submit_parser = subparsers.add_parser(
        "submit",
        help="Submit prepared jobs to Slurm.",
    )
    submit_parser.add_argument("plan", type=Path)
    submit_parser.add_argument(
        "--stage",
        choices=("baseline", "remaining", "all"),
        default="baseline",
    )

    worker_parser = subparsers.add_parser("worker", help=argparse.SUPPRESS)
    worker_parser.add_argument("plan", type=Path)
    worker_parser.add_argument(
        "--stage",
        choices=("baseline", "dose", "finalize"),
        required=True,
    )

    subparsers.add_parser(
        "diagnose",
        help="Show partition, resource and queue information.",
    )

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
