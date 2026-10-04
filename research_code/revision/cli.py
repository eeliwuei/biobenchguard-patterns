#!/usr/bin/env python3
"""Explicit reproduction modes. Fresh full runs rebuild inputs and fit models.

No mode publishes externally. Replotting and summarizing are never labelled full
reproduction. Every command requires a new output directory and records failures.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
CODE_ALLOWLIST = (
    "cli.py", "core.py", "data.py", "experiments.py", "summarize.py", "figures.py", "test_core.py",
    "e1/run_e1.py", "e1/report_e1.py", "e1/test_e1.py",
    "e2/run_e2.py", "e2/validate_saved_e2.py", "e2/fetch_public_crc.py",
)
PROTOCOL_ALLOWLIST = (
    "e1/protocol.json", "e1/freeze.json", "e1/original_to_replay_mapping.csv", "e1/source_audit.json",
    "e2/analysis_plan.md", "estimand_registry.yaml",
)
LOCKED_FIGURE_INPUTS = (
    "e1/main/real_failure_recovery.csv", "e3e4/operating_rates.csv",
    "e3e4/calibration_transfer_matrix.csv", "e3e4/mixed_response.csv",
    "e2/per_study_contrasts.csv", "e2/e2_results.json",
)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def require_file(path):
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise ValueError(f"Required input file is absent: {p}")
    return p


def require_directory(path):
    p = Path(path).expanduser().resolve()
    if not p.is_dir():
        raise ValueError(f"Required input directory is absent: {p}")
    return p


def import_source(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Run:
    def __init__(self, path, mode, arguments):
        proposed = Path(path).expanduser().absolute()
        if proposed.exists() or proposed.is_symlink():
            raise ValueError(f"Refusing to overwrite existing run directory: {proposed}")
        self.root = proposed.resolve()
        self.root.mkdir(parents=True)
        (self.root / "logs").mkdir()
        self.manifest = dict(mode=mode, status="RUNNING", started_utc=stamp(),
            arguments=arguments, orchestration_python=sys.executable, python=sys.version,
            platform=platform.platform(), cli_sha256=digest(__file__), commands=[], copied_sources=[])
        self.save()

    def save(self):
        (self.root / "cli_run_manifest.json").write_text(json.dumps(self.manifest, indent=2) + "\n")

    def command(self, name, command):
        command = [str(x) for x in command]
        item = dict(name=name, argv=command, started_utc=stamp(), status="RUNNING",
                    log=f"logs/{len(self.manifest['commands']):02d}_{name}.log")
        self.manifest["commands"].append(item); self.save()
        print(f"[{name}] starting; log: {self.root / item['log']}", flush=True)
        try:
            with (self.root / item["log"]).open("w") as log:
                # Respect explicitly supplied threading controls; defaults avoid nested BLAS oversubscription.
                env = os.environ.copy()
                for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
                    env.setdefault(key, "1")
                subprocess.run(command, cwd=self.root, stdout=log, stderr=subprocess.STDOUT,
                               env=env, check=True)
            item.update(status="COMPLETE", returncode=0, finished_utc=stamp())
        except subprocess.CalledProcessError as error:
            item.update(status="FAILED", returncode=error.returncode, finished_utc=stamp())
            raise
        except Exception as error:
            item.update(status="FAILED", returncode=None, finished_utc=stamp(), failure=repr(error))
            raise
        finally:
            self.save()

    def copy(self, source, relative, kind):
        source = require_file(source); target = self.root / relative
        if target.exists():
            raise ValueError(f"Refusing duplicate staged destination: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        self.manifest["copied_sources"].append(dict(source=str(source), destination=relative,
                                                    sha256=digest(source), kind=kind))
        self.save()

    def finish(self, status, **fields):
        self.manifest.update(status=status, finished_utc=stamp(), **fields); self.save()
        text = ["# Reproduction execution report", "", f"Mode: `{self.manifest['mode']}`.",
                f"Status: **{status}**.", "", self.manifest.get("scope", ""), "",
                "Commands, source hashes, interpreter paths, timestamps and logs are in `cli_run_manifest.json`.",
                "A completed numerical run does not establish clinical validity, a universal error rate, or a globally clean workflow."]
        (self.root / "REPRODUCTION_REPORT.md").write_text("\n".join(text) + "\n")


def stage_code(run):
    for relative in CODE_ALLOWLIST:
        run.copy(HERE / relative, "research_code/revision/" + relative, "analysis_source_code")
    return run.root / "research_code/revision"


def stage_protocols(run, source, experiments, config=None):
    source = require_directory(source)
    required = []
    if "e1" in experiments:
        required += ["e1/protocol.json", "e1/freeze.json", "e1/original_to_replay_mapping.csv", "e1/source_audit.json"]
    if "e2" in experiments:
        required += ["e2/analysis_plan.md"]
    for relative in required:
        run.copy(source / relative, "protocol/" + relative, "frozen_scientific_protocol")
    if (source / "estimand_registry.yaml").is_file():
        run.copy(source / "estimand_registry.yaml", "protocol/estimand_registry.yaml", "scientific_estimand_registry")
    if "e3e4" in experiments:
        if config is None:
            raise ValueError("E3/E4 requires an explicit --config, normally protocol/e3_e4_final.json")
        run.copy(config, "protocol/e3_e4_run.json", "requested_reference_configuration")


def stage_raw_links(run, raw_root=None, e1_data=None):
    raw = run.root / "data_manifest/raw"; raw.mkdir(parents=True)
    links = {}
    if raw_root:
        source = require_directory(raw_root)
        for name in ["crop", "microbiome", "rice"]:
            if (source / name).is_dir():
                (raw / name).symlink_to((source / name).resolve(), target_is_directory=True)
                links[name] = str((source / name).resolve())
    if e1_data:
        source = require_file(e1_data); (raw / "e1").mkdir()
        frozen = run.root / "protocol/e1/freeze.json"
        if frozen.is_file() and digest(source) != json.loads(frozen.read_text())["data_sha256"]:
            raise ValueError("E1 raw data hash differs from the supplied frozen protocol")
        (raw / "e1/parkinsons.data").symlink_to(source)
        links["e1/parkinsons.data"] = str(source)
    run.manifest["raw_input_links"] = links; run.save()
    return raw


def validate_inputs(raw_root, e1_data, config):
    import numpy as np
    import pandas as pd
    raw_root = require_directory(raw_root); config_path = require_file(config)
    config = json.loads(config_path.read_text())
    loader = import_source(HERE / "data.py", "revision_loader_validation")
    requested = {c["base"] for c in config["conditions"]}
    if config.get("mixed"):
        requested.update(["wheat", "microbiome_crc"])
    records = {}
    for base in sorted(requested):
        X, y, group, metadata = loader.load_base(base, raw_root)
        if X.ndim != 2 or len(X) != len(y) or len(X) != len(group) or not np.isfinite(X).all():
            raise ValueError(f"Invalid matrix/label/group shape or features: {base}")
        if len(set(metadata["sample_ids"])) != len(y):
            raise ValueError(f"Duplicate sample IDs: {base}")
        if base != "rice" and not np.isfinite(y).all():
            raise ValueError(f"Missing genuine outcomes: {base}")
        if base == "microbiome_crc" and set(np.unique(y)) != {0, 1}:
            raise ValueError("CRC original outcome labels are not binary")
        for condition in [c for c in config["conditions"] if c["base"] == base]:
            if condition.get("sample_limit", len(X)) > len(X) or condition.get("feature_limit", X.shape[1]) > X.shape[1]:
                raise ValueError(f"Configuration limit exceeds input dimensions: {condition['id']}")
        records[base] = dict(shape=list(X.shape), n_groups=len(set(group)), raw_files=metadata["raw_files"],
                             outcome_scope=metadata["phenotype"], preprocessing_scope=metadata.get("processing", metadata.get("group_source")))
    e1 = require_file(e1_data); d = pd.read_csv(e1)
    features = [c for c in d.columns if c not in ["name", "status"]]
    groups = d.name.str.split("_").str[2]
    if len(d) != 195 or len(features) != 22 or groups.nunique() != 32 or not d.name.is_unique:
        raise ValueError("UCI Parkinson input does not match the frozen E1 schema")
    if not np.isfinite(d[features].to_numpy(float)).all() or d.status.sum() != 147:
        raise ValueError("E1 features/outcomes do not match the frozen schema")
    if digest(e1) != "a36116c2deac07895b79ce97a2f99a1b3720ac4892995ff5c3b9ccb397ed2bc5":
        raise ValueError("E1 raw content hash differs from the audited UCI Parkinson file")
    if d.assign(subject=groups).groupby("subject").status.nunique().max() != 1:
        raise ValueError("E1 status is inconsistent within subject")
    records["e1"] = dict(shape=[195,22], n_subjects=32, n_case_records=147, raw_sha256=digest(e1))
    return dict(status="PASS", validated_utc=stamp(), config_sha256=digest(config_path),
                input_scope="actual raw parsing and finite/schema checks; no predictive-model fitting; rice outcomes are intentionally not used", datasets=records)


def python_for(args, name):
    value = getattr(args, name + "_python", None)
    if not value:
        return sys.executable
    # Preserve the venv executable path: resolving its symlink can discard the venv.
    path = Path(value).expanduser().absolute()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"Python executable is absent or not executable: {path}")
    return str(path)


def fit_e1(run, code, python):
    run.command("e1_input_contract", [python, code / "e1/test_e1.py"])
    run.command("e1_main_actual_fitting", [python, code / "e1/run_e1.py", "--stage", "main"])
    run.command("e1_report_from_new_predictions", [python, code / "e1/report_e1.py"])


def fit_e2(run, code, python):
    run.command("e2_integrity_pilot", [python, code / "e2/run_e2.py", "--root", run.root, "--pilot"])
    run.command("e2_main_actual_fitting", [python, code / "e2/run_e2.py", "--root", run.root])
    run.command("e2_independent_prediction_checks", [python, code / "e2/validate_saved_e2.py", "--root", run.root])


def fit_reference(run, code, python, workers):
    derived = run.root / "data_manifest/derived"
    run.command("rebuild_derived_from_raw", [python, code / "data.py", "--raw", run.root / "data_manifest/raw", "--out", derived])
    raw_run = run.root / "results_raw/e3e4"
    run.command("e3e4_main_actual_fitting", [python, code / "experiments.py", "--config", run.root / "protocol/e3_e4_run.json",
                                            "--data", derived, "--out", raw_run, "--workers", workers])
    run.command("summarize_new_reference_instances", [python, code / "summarize.py", "--run", raw_run,
                                                     "--out", run.root / "results_summary/e3e4"])


def parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version="patterns-revision-cli 14.0")
    sub = ap.add_subparsers(dest="mode", required=True)
    def new(name, help_text):
        p = sub.add_parser(name, help=help_text)
        p.add_argument("--run-root", required=True, help="New output directory; any existing path is rejected")
        return p
    p = new("validate-inputs", "Parse supplied real raw data and verify schemas, dimensions, labels and hashes")
    p.add_argument("--raw-root", required=True); p.add_argument("--e1-data", required=True); p.add_argument("--config", required=True)
    new("smoke-test", "Fit five tiny synthetic engineering reference instances; not scientific or full reproduction")
    new("demo", "Alias for smoke-test: five synthetic engineering computations; no scientific operating-rate inference")
    p = new("replot-locked", "Render figures from explicitly supplied existing source tables; no model refitting")
    p.add_argument("--source-root", required=True)
    p = new("summarize", "Recompute E3/E4 summaries from a complete actual run; no model refitting")
    p.add_argument("--reference-run", required=True)
    for name, help_text in [("run-experiment", "Execute one real analysis family in a fresh workspace"),
                            ("reproduce-full", "Rebuild matrices, fit all E1/E2/E3/E4 analyses, summarize and replot in a fresh workspace")]:
        p = new(name, help_text)
        if name == "run-experiment": p.add_argument("--experiment", choices=["e1", "e2", "e3e4"], required=True)
        p.add_argument("--raw-root", required=name == "reproduce-full")
        p.add_argument("--e1-data", required=name == "reproduce-full")
        p.add_argument("--protocol-root", required=True)
        p.add_argument("--config", required=name == "reproduce-full", help="Explicit reference configuration; use frozen e3_e4_final.json for the full revision")
        p.add_argument("--workers", type=int, default=4)
        for family in ["e1", "e2", "reference"]:
            p.add_argument("--" + family + "-python", help="Optional explicit Python executable for this analysis family; default is invoking Python")
    return ap


def main():
    args = parser().parse_args()
    if getattr(args, "workers", 1) < 1:
        raise SystemExit("--workers must be positive")
    run = None
    try:
        run = Run(args.run_root, args.mode, vars(args))
        if args.mode == "validate-inputs":
            result = validate_inputs(args.raw_root, args.e1_data, args.config)
            (run.root / "input_validation.json").write_text(json.dumps(result, indent=2) + "\n")
            run.finish("COMPLETE", scope="Raw-input parsing and schema validation only; no model fitting or full reproduction.")
        elif args.mode in ("smoke-test", "demo"):
            import numpy as np
            code = stage_code(run); data = run.root / "engineering_inputs"; data.mkdir()
            rng = np.random.default_rng(881); np.savez(data / "wheat.npz", X=rng.normal(size=(80,25)), g=np.repeat(np.arange(4),20))
            config = {"status":"ENGINEERING_SMOKE_ONLY", "conditions":[dict(id="engineering_smoke",base="wheat",task="regression",family="linear",signal=.5,k=5,folds=2,alternative_L=[4]) ]}
            path = run.root / "smoke_config.json"; path.write_text(json.dumps(config, indent=2))
            target = run.root / "results_raw/smoke_reference"
            run.command("five_actual_smoke_fits", [sys.executable, code / "experiments.py", "--config", path, "--data", data, "--out", target, "--workers", "1", "--pilot"])
            manifest = json.loads((target / "run_manifest.json").read_text())
            if not (manifest["status"] == "COMPLETE" and manifest["pilot"] and manifest["n_attempted"] == 5 and manifest["n_failed"] == 0):
                raise ValueError("Engineering smoke run did not complete all five intended instances")
            run.finish("COMPLETE", scope="Five actual tiny synthetic engineering computations; pilot data cannot be used for manuscript claims or described as full reproduction.")
        elif args.mode == "replot-locked":
            source = require_directory(args.source_root); code = stage_code(run)
            for relative in LOCKED_FIGURE_INPUTS:
                run.copy(source / "results_summary" / relative, "results_summary/" + relative, "locked_figure_source_table_not_refitted")
            run.command("replot_existing_tables", [sys.executable, code / "figures.py", "--root", run.root])
            run.finish("COMPLETE", scope="Figures regenerated from copied locked source tables. No models or reference instances were refitted.")
        elif args.mode == "summarize":
            source = require_directory(args.reference_run); code = stage_code(run)
            run.command("summarize_existing_reference_run", [sys.executable, code / "summarize.py", "--run", source, "--out", run.root / "results_summary/e3e4"])
            run.finish("COMPLETE", scope="Reference summaries recomputed from an explicitly supplied existing complete run. No model refitting.")
        else:
            experiments = ["e1","e2","e3e4"] if args.mode == "reproduce-full" else [args.experiment]
            if args.mode == "reproduce-full":
                final_config = require_directory(args.protocol_root) / "e3_e4_final.json"
                if digest(require_file(args.config)) != digest(require_file(final_config)):
                    raise ValueError("reproduce-full requires the unchanged protocol-root/e3_e4_final.json; use run-experiment for a different reference configuration")
            if "e1" in experiments and not args.e1_data: raise ValueError("E1 requires --e1-data")
            if any(x in experiments for x in ["e2","e3e4"]) and not args.raw_root: raise ValueError("E2 or E3/E4 requires --raw-root")
            code = stage_code(run); stage_protocols(run, args.protocol_root, experiments, args.config)
            stage_raw_links(run, args.raw_root, args.e1_data)
            if "e1" in experiments: fit_e1(run, code, python_for(args,"e1"))
            if "e2" in experiments: fit_e2(run, code, python_for(args,"e2"))
            if "e3e4" in experiments: fit_reference(run, code, python_for(args,"reference"), args.workers)
            if args.mode == "reproduce-full":
                run.command("figures_from_new_results", [sys.executable, code / "figures.py", "--root", run.root])
            run.finish("COMPLETE", scope="Actual new model fitting from explicitly supplied raw inputs; reference matrices rebuilt where E3/E4 requested. No saved predictions or scientific summary tables were copied as results.", experiments_completed=experiments)
        print(f"COMPLETE: {run.root}", flush=True)
    except Exception as error:
        if run is not None:
            run.finish("FAILED", failure_type=type(error).__name__, failure=str(error))
        print(f"FAILED: {error}", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            raise SystemExit(error.returncode or 1)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
