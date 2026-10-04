# E1 post hoc reference interpretation

`reanalyze_e1_references.py` recomputes an exploratory interpretation of already
locked predictions. It was added after the original E1 execution; it is not a
new prespecified hypothesis, a new permutation method or additional model fitting.
The original first split is fixed at repeat 0 and seed 110003. No split search is
performed. The observed statistic is recording-weighted pooled OOF AUROC for that
single split, not the mean over the 30 original splits.

## Inputs and command

Supply a completed original E1 directory containing all three files:

- `oof_predictions.csv`
- `permutation_predictions.csv.gz`
- `permutation_controls.csv`

Run from the code-package root, using the numerical environment described in
`reproduction/README_REPRODUCTION.md` and an output directory that does not exist:

```bash
.venv/bin/python research_code/revision/e1/reanalyze_e1_references.py \
  --input-dir /path/to/original/results_raw/e1/main \
  --out runs/e1-exploratory
```

The primary distributed results come from the locked original predictions. A new
E1 run can supply compatible files, but numerical differences across environments
must retain that new run's identity rather than be labeled the original results.
The script validates the complete combinations, row IDs, original labels,
permutation seed rule, label-vector hashes and saved AUROCs. It independently
checks AUROC with exact pairwise ranks and scikit-learn before counting tails.
It writes eight comparison rows and 792 reference-statistic rows; it copies no
sample-level predictions into these outputs.

## Reference questions and descriptive boundary

For each of the eight model × split × reference combinations, the reported rank is
`(1 + count(reference AUROC >= observed AUROC)) / (99 + 1)`. Equality contributes
to the upper tail. The optional descriptive line uses **rank < 0.05**, strictly.
A value equal to 0.05 is identified separately and is not classified below it.
There is no multiple-comparison correction and no calibrated-p-value claim.
All eight combinations are retained:

| Model and split | Observed first-split AUROC | Record-label rank | Subject-label rank |
|---|---:|---:|---:|
| RF, record-wise | 0.9064625850340137 | 0.01 | 0.12 |
| Logistic, record-wise | 0.8808106575963719 | 0.01 | 0.02 |
| RF, subject-wise | 0.7949263038548752 | 0.01 | 0.01 |
| Logistic, subject-wise | 0.7073412698412698 | 0.01 | 0.05 |

Record-label randomization breaks within-subject label consistency. Subject-label
randomization retains that consistency while assigning outcomes across subjects.
They address different reference questions; differing ranks are not contradictory
results from the same validated null hypothesis. The original procedure redrew
labels when a fixed training or test fold had only one class. The ranks therefore
condition on that acceptance rule. Without a demonstrated conditional
exchangeability argument, they remain conditional Monte Carlo tail ranks rather
than formally calibrated significance-test p values.

The RF record-wise result is extreme relative to record-label randomization but
less extreme relative to identity-preserving randomization. This limits the
inference that exceeding the first reference alone establishes disease-specific
or new-subject generalization evidence. It motivates appropriate independent-
subject validation; it does not justify a rule that rank < 0.05 authorizes
clinical deployment. Rank 0.12 does not show absence of disease signal or prove
that subject identity explains every part of performance. Neither reference
certifies a clean benchmark or estimates universal sensitivity, specificity,
leakage FPR or a detection limit.

The interpretation follows established identity-confounding and subject-level
randomization principles, including Chaibub Neto et al. (2019),
DOI `10.1038/s41746-019-0178-x`; it is not claimed as a newly invented test.

## Distributed data boundary

`results_summary/e1/exploratory/` contains the eight-row comparison table and the
792-row reference-AUROC table. These are aggregate statistics over recordings,
not individual recordings, labels, subject IDs or OOF probabilities. The original
locked prediction/control inputs must be supplied separately and are not included
in the code archive. The two summary tables alone cannot reconstruct or rerun the
sample-level AUROC calculation.
