# Reproduce the analyses

Run commands from the package directory. The main CLI entry point is `research_code/revision/cli.py`. Its package version identifies this command interface, not a research model or an AI model. Each CLI mode requires a new explicit `--run-root`; existing directories and symlinks are rejected. The standalone E5 summarizer instead reads its completed `--root` and uses the publication transaction documented below. Failures propagate as nonzero exits.

## Environment

Python 3.12 is required for the tested numerical pins. Create a new environment:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r reproduction/independent/requirements_main_pinned.txt
.venv/bin/python -m pip install pytest==9.0.2
.venv/bin/python research_code/revision/cli.py --help
.venv/bin/python research_code/revision/cli.py --version
```

The reference environment pins NumPy 1.26.4, pandas 2.3.3, SciPy 1.17.1 and scikit-learn 1.8.0. Original E1/E2 macOS runs used Python 3.12.14, NumPy 2.5.3, pandas 3.0.5 and scikit-learn 1.9.0, captured in `reproduction/requirements.lock`. That broad lock also contains document-production packages. E1/E2 RF values differ in the reference environment; installing its pins is not a promise of equality with the primary real-label tables. See the scientific reproduction report.

## Tier 1: regenerate figures from recorded numerical results

The included `results_summary/` contains the canonical figure inputs. E1/E2 tables use the original real-label runs; E3/E4 tables use the primary Linux reference run. No raw data or model fitting is needed:

```bash
.venv/bin/python research_code/revision/cli.py replot-locked \
  --source-root . --run-root runs/replot
```

This creates seven figures in PDF, SVG and PNG and records hashes of every copied source table. Replotting is a presentation check, not numerical reproduction by refitting.

## Engineering checks and a tiny fitted demonstration

```bash
.venv/bin/python -m pytest -q -p no:cacheprovider research_code/revision/test_core.py
.venv/bin/python -m pytest -q -p no:cacheprovider research_code/revision/e5/tests/test_summarize_e5.py
.venv/bin/python research_code/revision/cli.py demo --run-root runs/demo
```

`demo` is an alias of `smoke-test`: five tiny seeded synthetic reference computations test the execution path. They do not estimate biological operating rates. `summarize` recomputes tables from a supplied complete E3/E4 raw run and rejects a pilot run; it does not refit models.

The v15 run passed seven core tests, 64 E5 regression/fault-injection cases and the five-instance demo. The E5 suite also passed with only the package's canonical derived tables and the small recorded manifest/bootstrap fixtures; no private input path is required. Its NPZ fixture holds 12 derived threshold/rate/rate-difference bootstrap arrays, not raw outcomes or sample-level predictions. The corrected E5 summarizer's normal table and all 72,000 bootstrap values equal v14 exactly. [E5 validation and transaction semantics](../research_code/revision/e5/README_SUMMARY_VALIDATION.md) describe explicit failures, staging, backup/rollback and the limits of multi-file publication.

## Post hoc E1 reference interpretation from locked predictions

This optional summary reads the original completed E1 files `oof_predictions.csv`, `permutation_predictions.csv.gz` and `permutation_controls.csv`; these sample-level prediction files are not bundled. Supply the locked original run to reproduce the distributed primary comparison. A new E1 refit may differ numerically across environments and must retain its own identity.

```bash
E1_LOCKED=/path/to/original/results_raw/e1/main
.venv/bin/python research_code/revision/e1/reanalyze_e1_references.py \
  --input-dir "$E1_LOCKED" --out runs/e1-exploratory
```

The new output directory receives eight comparison rows and 792 reference-AUROC rows. No models are fitted and no input predictions are copied into those tables. The fixed first split is repeat 0, seed 110003; the statistic is recording-weighted pooled OOF AUROC. The 30-repeat mean is not substituted for this single-split observation. All eight combinations are retained. The conditional Monte Carlo rank is `(1 + count(reference AUROC >= observed AUROC)) / 100`. Record-label and identity-preserving subject-label references address different questions. The descriptive line is strict rank < 0.05; equality is shown separately. There is no multiplicity correction or calibrated-p-value claim. See [full interpretation limits](../research_code/revision/e1/README_EXPLORATORY_REFERENCES.md).

## Tier 2: refit reference models from processed matrices

For E3/E4, supply the processed `wheat.npz`, `microbiome_crc.npz` and `rice.npz` matrices produced by the included `data.py` loader. These matrices are not bundled because they preserve the underlying data. This mode fits the complete recorded reference configuration but does not repeat raw-file parsing or the original-label E1/E2 analyses:

```bash
PROCESSED=/path/to/derived-matrices
.venv/bin/python research_code/revision/experiments.py \
  --config protocol/e3_e4_final.json --data "$PROCESSED" \
  --out runs/from-processed --workers 4
.venv/bin/python research_code/revision/summarize.py \
  --run runs/from-processed --out runs/processed-summary
```

The recorded full execution includes this same processed-matrix fitting stage. The current packaging checks do not repeat its 14,960 evaluations. Use Tier 3 to rebuild matrices and run the original-label families from the specified source files.

## Tier 3: actual raw-to-result recomputation

Obtain authorized source inputs with this directory structure. The raw files are not distributed in this package:

```text
raw/
  crop/data/wheat.RData
  microbiome/data/crc_zeller_results/crc_zeller.metadata.txt
  microbiome/data/crc_zeller_results/RDP/crc_zeller.otu_table.100.denovo.rdp_assigned
  microbiome/data/crc_baxter_results/crc_baxter.metadata.txt
  microbiome/data/crc_baxter_results/RDP/crc_baxter.otu_table.100.denovo.rdp_assigned
  microbiome/data/crc_zackular_results/crc_zackular.metadata.txt
  microbiome/data/crc_zackular_results/RDP/crc_zackular.otu_table.100.denovo.rdp_assigned
  rice/data/RiceDiversity_44K_Genotypes_PLINK/sativas413.ped
  rice/data/germplasm.csv
```

Supply UCI `parkinsons.data` separately. Its expected SHA256 is `a36116c2deac07895b79ce97a2f99a1b3720ac4892995ff5c3b9ccb397ed2bc5`. Set `RAW` and `E1_DATA` to those input locations; each output directory below must be new:

```bash
RAW=/path/to/raw
E1_DATA=/path/to/parkinsons.data
.venv/bin/python research_code/revision/cli.py validate-inputs \
  --raw-root "$RAW" --e1-data "$E1_DATA" \
  --config protocol/e3_e4_final.json --run-root runs/input-check
.venv/bin/python research_code/revision/cli.py reproduce-full \
  --raw-root "$RAW" --e1-data "$E1_DATA" --protocol-root protocol \
  --config protocol/e3_e4_final.json --workers 4 --run-root runs/full
```

`reproduce-full` rebuilds the three reference matrices, actually fits E1, E2 and all 14,960 E3/E4 reference evaluations, recomputes summaries, checks E2 stored predictions and draws figures from the new outputs. It never substitutes included summary tables for new fitted results. The configuration must match `protocol/e3_e4_final.json` exactly. Raw inputs are linked into the fresh run and must remain available; sources are read without modification. The E1 frozen code/protocol hashes are checked.

The default interpreter runs every family. To use separate environments, set `--e1-python`, `--e2-python` and `--reference-python` to the corresponding executable paths. The original real-label and reference pins describe different environments and must not be silently conflated. Each analysis records its actual package versions.

To run only one family:

```bash
.venv/bin/python research_code/revision/cli.py run-experiment --experiment e1 \
  --e1-data "$E1_DATA" --protocol-root protocol --run-root runs/e1
.venv/bin/python research_code/revision/cli.py run-experiment --experiment e2 \
  --raw-root "$RAW" --protocol-root protocol --run-root runs/e2
.venv/bin/python research_code/revision/cli.py run-experiment --experiment e3e4 \
  --raw-root "$RAW" --protocol-root protocol --config protocol/e3_e4_final.json \
  --workers 4 --run-root runs/reference
.venv/bin/python research_code/revision/cli.py summarize \
  --reference-run runs/reference/results_raw/e3e4 --run-root runs/reference-summary
```

The standalone E1 schema test requires the supplied real input and is invoked by the E1 run mode. It cannot run from this raw-free archive alone.

## Data access and rights

| Input | Access and terms | Use in these analyses |
|---|---|---|
| UCI Parkinsons | [UCI dataset 174](https://archive.ics.uci.edu/dataset/174/parkinsons), CC BY 4.0 | Original disease labels; retain attribution. The analyzed file has 195 recordings and 32 parsed subjects, differing from the current catalogue prose. |
| microbiomeHD CRC | [Zenodo 840333](https://zenodo.org/records/840333), CC BY-NC 4.0 | Three processed cohorts; preserve noncommercial terms and original cohort attribution. The fetcher below checks published archive MD5 and refuses inconsistent existing files. |
| BGLR wheat | [BGLR manual](https://cran.r-project.org/web/packages/BGLR/BGLR.pdf) | Obtain an authorized package/source copy; the analysis-code license does not replace source-data terms. |
| Rice diversity panel | RiceDiversity genotype and population metadata underlying Zhao et al. (2011), DOI 10.1038/ncomms1467 | An independent redistribution license was not established by the recorded reproduction audit. Obtain authorized inputs; these files are not bundled. |

```bash
.venv/bin/python research_code/revision/e2/fetch_public_crc.py --destination "$RAW"
```

The three source CRC archives and their six analyzed members were independently checksum-matched in the original data audit. The fetcher obtains only those public CRC inputs; it does not infer rights for rice or wheat. No publisher full texts or third-party author R code are redistributed.

## Interpretation and execution status

Numerical output is conditional on the recorded task, model, feature geometry and protocol. E3 L=0 baseline-reference exceedance is not clean-pipeline FPR. The real E2 342-genus representation differs from the 195-feature CRC reference geometry, so historical calibration is not transported into the real task. E1 reference annotations are generated separately from the numerical/provenance auditor and are not supplied to its calculation.

See `SCIENTIFIC_REPRODUCTION_REPORT.md` for completed historical model-refitting evidence and portability limits. `INTERFACE_VALIDATION.json`, `PRIOR_EXECUTION_EVIDENCE.json` and `E5_REPRODUCTION.json` retain v14 execution identities; their earlier summary-code hash is historical and does not identify the corrected v15 E5 summarizer. `V15_VALIDATION.json` records the narrower current checks. The v15 packaging pass does not rerun the full numerical program. Fitting implementations, scientific parameters and seeds remain unchanged; v15 changes to summary validation, post hoc interpretation and figure presentation are distinguished in `SOURCE_COMPARISON.md` and the final source comparison.
