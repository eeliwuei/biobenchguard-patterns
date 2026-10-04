# E5 validated summary publication (v15)

Run `python research_code/revision/e5/summarize_e5.py --root /path/to/completed/e5/run`.
The root must contain the original `run_manifest.json` and `scores.csv`. No raw-data
model refitting occurs. The unchanged E5 protocol/fitting source remains bundled.
The v15 summary code hash identifies an engineering correction made after v14.

Before statistical computation or output writes, the summarizer checks the
manifest identity/status/failure list and score checksum; the complete four-outcome,
two-role, 100-replicate grid; IDs and their declared seeds; row statuses; finite
R2_B/R2_C/T_scope; and T_scope = R2_C − R2_B with absolute tolerance 1e-12 and relative
tolerance 1e-10. It checks recorded model counts against the completed score grid
(two arms × five folds per score). This is consistency verification of recorded
fits, not independent reconstruction of those fits. Actual successful counts
supply all rate denominators. Incomplete runs are rejected rather than summarized.

The original and reader-exposition protocol hashes are explicitly distinguished
by the unchanged `protocol/e5/freeze.json`; either recorded run identity is accepted.
The bundled reader exposition itself must match its frozen hash. Score bytes are
read once, hashed and parsed from that same snapshot. Input and freeze files are
never edited.

All statistics finish in memory before tables and figures are rendered into a
temporary directory. Finite interval bounds and observed-rate containment are
validated before plotting. Only yerr roundoff within 1e-12 is normalized to zero;
real negative errors, reversed intervals, non-finite values and invalid bounds
raise an error. Normal-input CSV values and the 12 bootstrap arrays match v14.

## Publication transaction and its limits

Existing output filenames remain `operating_rates.csv`, `bootstrap_draws.npz`, the
three `Figure_S1_authentic_wheat` formats, and `summary.json`. A writer acquires
`.e5-summary.lock`. After rendering succeeds, it backs up existing outputs in its
temporary directory and writes `.e5-summary-transaction.json` before replacing any
formal output. Files are replaced individually with `os.replace`; `summary.json`
is last. Catchable commit errors restore all old files or remove new files when
no previous outputs existed. A rendering error leaves formal outputs untouched.
Tests inject first, middle and last replacement failures with and without old files.

**This is not an atomic six-file snapshot.** Individual replacements are atomic;
concurrent readers may observe a mixture during a commit. Readers must not treat
outputs as a completed set while the writer lock or transaction journal exists.
A process kill, power loss, persistent filesystem failure or rollback failure is
not claimed automatically recoverable. The journal and backup directory are kept
on rollback failure, and any later invocation refuses an unresolved journal/lock.
For manual recovery, first ensure no writer is active, inspect the journal's
`previous_files` and recorded hashes, restore every previously existing file from
its `stage_directory/backup`, remove names originally absent, and verify all hashes
before removing the journal/lock. Keep that recovery evidence for the audit.
Do not simply delete a stale journal and call the interrupted output successful.

The lock, staging files and recovery journal are operational metadata, not scientific
success outputs. Standard readers should wait for a successful exit and absence of
both markers. The transaction does not rewrite protocol freezes or any input files.

## Tests

From the code-package root, use a Python environment with the pinned reference
NumPy/pandas versions and pytest:

```
python -m pytest research_code/revision/e5/tests/test_summarize_e5.py -q
python -m pytest research_code/revision/test_core.py -q
python research_code/revision/cli.py demo --run-root /path/to/new/demo_run
```

Tests use bundled canonical derived scores/operating tables and the recorded v14
manifest/bootstrap fixtures. Set `E5_REFERENCE_ROOT` to a complete recorded E5 run
to test against its original summary, bootstrap arrays and rendered outputs.
Fixtures and faulted inputs are always copied into pytest temporary directories.
The fixture NPZ contains derived bootstrap summary arrays, not sample-level raw
inputs or fitted predictions. Tests and the five-instance CLI demo are engineering
validation, not an additional full E1–E5 reproduction.
