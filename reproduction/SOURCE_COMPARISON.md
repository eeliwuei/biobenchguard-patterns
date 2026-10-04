# Source comparison and execution identity

The source comparison distinguishes v15 engineering and post hoc changes from the
historical model-fitting evidence. Source files and the final package manifest
identify exact artifacts; a successful old execution does not mean a subsequently
edited summary script was used in that execution.

| v15 change | Scope and evidence |
|---|---|
| `e5/summarize_e5.py` | New complete manifest/score validation; actual denominators; finite interval checks; staged rendering; per-file replacement with backup/rollback, writer lock and recovery journal. This is a post-v14 engineering correction. Normal 16-row CSV is byte-identical to v14 and all 12 bootstrap arrays (72,000 values) agree exactly. |
| `e5/tests/` and its README | 64 regression/fault cases passed against original recorded inputs and separately with bundled fixtures only. Seven existing core tests and the five-instance CLI demo passed. These tests do not refit the biological experiments. |
| `e5/tests/fixtures/bootstrap_draws_v14.npz` | An explicit derived-data exception to excluding sample-level NPZ files: only threshold, rate and rate-difference bootstrap arrays are included. No feature matrix, biological outcomes, OOF predictions, sample IDs or selected predictors are present. The companion manifest is the recorded input identity for those regression fixtures. |
| `e1/reanalyze_e1_references.py` and two exploratory CSVs | New post hoc descriptive recomputation from the fixed first E1 split, retaining all eight combinations and all 792 reference AUROCs. The script and tables are copied from the executed portable reader analysis. No model fitting, seed search or revised original freeze. |
| `figures.py`, calibration-transfer figure (Figure 5 in the v15 numbering; Figure 6 in the current revision) | A presentation change to keep the complete colorbar label inside the canvas. The matrix values, cell ordering, color scale and scientific conclusions are retained. |
| README and data-boundary explanations | Identify each execution period, describe usable commands, preserve the reference interpretation limits and distinguish aggregate test fixtures from excluded sample-level inputs. |

The E1, E2 and E5 fitting scripts, the shared core/data/reference-experiment code,
scientific configuration, seed rules and protocol/freeze files are inherited
unchanged. The E3/E4 summary implementation is unchanged. This statement does
**not** include the corrected E5 summary implementation or the newly added E1
post hoc interpreter.

`PRIOR_EXECUTION_EVIDENCE.json` and `E5_REPRODUCTION.json` preserve historical full
execution evidence. `INTERFACE_VALIDATION.json` preserves the v14 environment and
interface checks. Their hashes and dates identify those past sources. In
particular, the earlier E5 summary-code hash in the historical report must not be
replaced with the v15 summary hash to imply retrospective execution.
`V15_VALIDATION.json` states the actual narrower current checks. No additional
full E1–E5 raw-input refit is claimed for v15.

The publication transaction is not an atomic snapshot across six filenames.
All statistics and rendering finish before formal outputs change. Catchable
commit failures are rolled back using backups, while a persistent rollback
failure leaves a journal and prevents another run until recovery. Concurrent
readers and abrupt process/power loss require the documented precautions in the
E5 validation README; file-by-file `os.replace` is not presented as group atomic.

Raw biological source files, sample-level reconstructed matrices, prediction
files, private working histories, audit instructions and software caches are
excluded. The original E1 prediction inputs are required separately to run the
post hoc interpreter. The numerical source tables and the explicitly described
bootstrap regression fixture are retained for reader checks.
