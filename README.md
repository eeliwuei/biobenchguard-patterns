# Same score, different meaning: data-leakage checks in biological machine learning

Code S1 of Liu and Jang, *Same score, different meaning: interpreting data-leakage checks in biological machine learning*, Patterns (2026; submitted under the title "Controlled reference tasks characterize benchmark diagnostics in biological prediction"). The package holds the analysis code, fixed configurations, environment specifications, execution records and data-access instructions. The derived source tables of Data S1 are included under `results_summary/`.

## Archive and versions

The package is archived at Zenodo.

| Version | DOI | Notes |
|---|---|---|
| All versions | https://doi.org/10.5281/zenodo.23124515 | Resolves to the latest version |
| 1.0.1 (this version) | https://doi.org/10.5281/zenodo.23133354 | Documentation and license update |
| 1.0.0 | https://doi.org/10.5281/zenodo.23124516 | First archived version |

Both versions contain the same code, protocols and numerical results. Version 1.0.1 differs from 1.0.0 in five files: `README.md`, `NOTICE.md`, `LICENSE-DATA` (new), `reproduction/SOURCE_COMPARISON.csv` and `PACKAGE_MANIFEST.csv`. In `SOURCE_COMPARISON.csv`, the rows for `README.md`, `NOTICE.md` and `research_code/revision/figures.py` carry new hashes and extended notes, and a row for `LICENSE-DATA` was added; the recorded hash of `figures.py` now matches the file shipped in both versions; version 1.0.0 still recorded the hashes of `figures.py` and `README.md` from before the figure renumbering. Every file under `research_code/`, `protocol/` and `results_summary/`, and every other file under `reproduction/`, is byte-identical in the two versions. Inside the zip archive the top-level folder is now `Code_S1/` (1.0.0: `Code_S1_v16_figmap/`). No analysis was rerun for version 1.0.1.

## Reproduction

Start with [reproduction instructions](reproduction/README_REPRODUCTION.md) for installation, exact E1–E4 commands, input layout, expected outputs and portability limits. The instructions separate three reproduction tiers: Tier 1 redraws main Figures 1–7 from the stored source tables (Figures 1 and 2 are schematics with fixed text), Tier 2 refits the E3/E4 reference models from processed matrices, and Tier 3 recomputes E1–E4 from the original source files. E5 and Figure S1 have their own commands (see below), and the E1 exploratory comparison behind Tables 2–3 and Table S7 requires the original locked E1 prediction files. [Data access](DATA_ACCESS.md) lists acquisition routes and expected source hashes.

## Correspondence with the published article

| Article item | Location in this package |
|---|---|
| Figures 1–7 | `research_code/revision/figures.py`; Tier 1 of the reproduction instructions |
| Figure 3 and the E1 rows of Table S5 | `results_summary/e1/main` |
| Tables 2–3 and Table S7 | `results_summary/e1/exploratory` |
| Figure 4 and the E2 rows of Table S5 | `results_summary/e2` |
| Figures 5–7 and Tables S2–S4 | `results_summary/e3e4` |
| Table 4, Table S6 and Figure S1 | `results_summary/e5`; `research_code/revision/e5/summarize_e5.py` |
| Table S8 | Cross-environment comparisons recorded in `reproduction/SCIENTIFIC_REPRODUCTION_REPORT.md` and `reproduction/PRIOR_EXECUTION_EVIDENCE.json`; dependency pins in `reproduction/requirements.lock` |

Table 1 (notation) and Table S1 (related methodological work) have no source file. Tier 1 reproduces the published Figures 4–7 and gives Figures 1–3 in their earlier layout. For the published article, Figures 1 and 2 were redrawn with wider spacing between boxes, full arrow shafts and heads and slightly larger type, and the legend of Figure 3A was moved into the empty area between the two groups of points. These layout edits were made outside this package; the text of Figures 1 and 2 and all plotted values are those produced by `figures.py`.

The files keep the field names under which they were written. The article uses the notation of its Table 1:

| In the files | In the article |
|---|---|
| `n_calibration`, `calibration_n`; `B` inside `quantile_rule` strings and in the `quantile` entry of `protocol/estimand_registry.yaml` | Number of calibration instances, N_cal (Equation 3) |
| `cutoff`, `threshold`; `calibration_rank` | Calibrated cutoff τ; rank j of the calibration order statistic |
| `fixed_cutoff`, `fixed_threshold`, `illustrative_threshold` | Analyst-defined comparators (0.05 score units in E3 and E5; 0.02 in E2), not the calibrated cutoff τ |
| `exceedance_k`, `fixed_k`, `n_positive`, `source_eval_k`, `observed_k`, `fixed_n_positive` | Exceedance counts: the "Calibrated count/100" and "Fixed count/100" columns of Tables 4 and S6, the numerators of the rates in Tables S2 and S3, and the fixed-comparator counts behind Figure 5A |
| `n_reference_ge_observed` (per draw: `reference_ge_observed`) in `results_summary/e1/exploratory` | "Count/99" in Tables 3 and S7: the number of reference AUROCs at least as large as the observed AUROC (Supplemental Methods S14 calls it column `k`) |
| `k` in `results_summary/e3e4/reference_results.csv` | Feature budget k |
| `score_B`, `score_C`, `R2_B`, `R2_C` | Metric under workflows B and C |
| `arm` in `results_summary/e2` (`A_variance20`, `B_train_supervised20`, `C_full_supervised20`, `D_all_train_variable`) | Feature arms A–D of Supplemental Methods S4 (Document S1) |
| `bootstrap_B`, `n_bootstrap` | Number of bootstrap resamples; unrelated to workflow B |
| `phenotype` (1, 2, 4, 5) | Observed wheat outcome labels |

Protocol files, execution records and historical notes keep the figure and panel labels that were current when they were written, for example "old Figure 3A" in `protocol/e3_e4_final.json` and "Figure 5" for the calibration-transfer figure in the v15 numbering. Their hashes are recorded and they were not relabeled.

## Licenses

All files outside `results_summary/` other than the license texts `LICENSE` and `LICENSE-DATA` (the analysis code and its test fixtures, the protocols, the execution and reproduction records, the manifests and the documentation) are under the Apache License 2.0 (`LICENSE`). The derived numerical tables under `results_summary/`, which are the tables of Data S1, are under the Creative Commons Attribution-NonCommercial 4.0 International license (CC BY-NC 4.0; `LICENSE-DATA`; Copyright 2026 DeepFarm Inc.): they may be reused with credit to the article for noncommercial purposes only. Third-party source datasets are not included and keep their own terms; see `NOTICE.md` and `DATA_ACCESS.md`.

## Package history

Labels v14–v16 denote internal revisions made before archiving; archive version 1.0.0 is the v16 package. Numerical results and analysis code are unchanged from the v15 package; for the revised manuscript only the main-figure numbering was synchronized (see FIGURE_NUMBERING_UPDATE.md: Figure 4 = CRC decision, Figure 5 = reference operating rates, Figure 6 = calibration transfer, Figure 7 = mixed mechanisms; Figures 1–3 and Figure S1 kept their numbers) and the package manifest was regenerated. Historical notes below keep the figure numbers used at the time of their execution.

The v15 package added E5 input validation and safe publication, an explicitly post hoc E1 reference comparison, and a layout correction to the calibration-transfer figure (Figure 5 in the v15 numbering; Figure 6 in the published article). The original fitting programs, seeds and protocol freezes are unchanged. The CLI still reports interface version 14.0; that interface label does not identify the revision of every source file. See [source comparison](reproduction/SOURCE_COMPARISON.md).

## Authentic wheat outcomes (E5)

The complete E1–E4 CLI run and E5 are separate commands. E5 retains all four authentic standardized yield outcomes. Obtain the exact BGLR source object listed in the data guide and run from this directory after environment installation:

```sh
.venv/bin/python research_code/revision/e5/run_e5.py --raw /path/to/wheat.RData --protocol protocol/e5/scientific_protocol.md --freeze protocol/e5/freeze.json --out runs/e5 --workers 4
.venv/bin/python research_code/revision/e5/summarize_e5.py --root runs/e5
```

The fitting command's `--out` directory must not already exist; the subsequent summarizer reads that completed run. The fitting command performs 20,000 Ridge fits and writes 2,000 scores plus per-reference OOF predictions, fold IDs and selected predictors. The summary writes all 16 outcome/level rows, conditional intervals and Figure S1. Distributed E5 results include derived scores; raw outcomes and sample-level reconstructed matrices/predictions are excluded. The protocol file is a reader edition of the recorded execution freeze with unchanged scientific parameters; its distinct document hash is explicit. E5 randomizes partitions/candidates on fixed authentic data. It is not a clean-workflow biological null and did not establish additional practitioner utility over correct provenance inspection.

A second process in the same pinned environment reproduced all 1,000 saved arrays exactly across 200 reference files. This is computational repeatability, separate from the E1–E4 cross-platform comparisons.

That full E5 fitting evidence is inherited from v14. The v15 correction reran the summarizer, including all three Figure S1 formats: its 16-row operating CSV is byte-identical and all 12 bootstrap arrays (72,000 values) equal v14 exactly. It did not repeat the biological model fits. The hardened summarizer rejects incomplete manifests, failed rows, bad checksums, invalid grids and non-finite scores before statistical computation. Its staged rendering and backup/rollback publication preserve previous outputs on the tested catchable failures. Individual replacements are atomic; the six-file set is not a concurrent-reader or power-loss atomic snapshot. See [the validation and transaction contract](research_code/revision/e5/README_SUMMARY_VALIDATION.md).

## Tests and the E1 exploratory comparison

The v15 checks passed 64 E5 regression/fault-injection cases, seven existing core tests and five tiny CLI reference computations. The 64 cases also pass using only bundled derived fixtures. These are engineering checks, not another full E1–E5 raw-data refit. The small NPZ under `e5/tests/fixtures` contains only the 12 derived bootstrap arrays; it contains no sample-level outcomes, feature matrices or model predictions.

The [E1 exploratory instructions](research_code/revision/e1/README_EXPLORATORY_REFERENCES.md) explain the added reader script and all eight model/split/reference comparisons. It uses the already fixed first split (repeat 0, seed 110003), computes conditional Monte Carlo tail ranks, and keeps the two reference questions distinct. The comparison is post hoc, uses a strict rank < 0.05 descriptive line, applies no multiplicity correction and does not supply calibrated p values. The two distributed CSVs contain aggregate reference statistics only. Running the script requires the original locked E1 prediction/control files, which are not distributed here.
