#!/usr/bin/env python3
"""Post hoc descriptive first-split E1 reference comparison; no model fitting.

Input directory: original E1 oof_predictions.csv and
permutation_predictions.csv.gz and permutation_controls.csv, as produced by the E1 reproduction pipeline.
Outputs: eight summary rows and 792 reference-statistic rows. No raw predictions
are copied. This is not a calibrated significance test or clinical validation.

Example:
  python reanalyze_e1_references.py --input-dir /path/to/results_raw/e1/main \
      --out /path/to/new/e1_reference_summary
"""
import argparse
import hashlib
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SEED = 110003
N_REFERENCE = 99
MODELS = {'random_forest', 'logistic_regression'}
SPLITS = {'record_wise', 'subject_wise'}
REFERENCES = {'record_label', 'subject_label'}
SCOPE = 'POST_HOC_DESCRIPTIVE_NOT_CALIBRATED_P_VALUE'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact_auc(labels, probabilities):
    y = np.asarray(labels, dtype=int)
    p = np.asarray(probabilities, dtype=float)
    require(set(y) == {0, 1}, 'Both outcome classes are required')
    require(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(),
            'Predictions must be finite probabilities in [0, 1]')
    pos = p[y == 1, None]
    neg = p[None, y == 0]
    twice_wins = int(2 * (pos > neg).sum() + (pos == neg).sum())
    return Fraction(twice_wins, 2 * pos.shape[0] * neg.shape[1])


def run(input_dir, out):
    require(not out.exists(), f'Output already exists: {out}')
    real = pd.read_csv(input_dir / 'oof_predictions.csv', float_precision='round_trip')
    null = pd.read_csv(input_dir / 'permutation_predictions.csv.gz', float_precision='round_trip')
    controls = pd.read_csv(input_dir / 'permutation_controls.csv', float_precision='round_trip')
    real_columns = {'model', 'scheme', 'repeat', 'seed', 'row_index',
                    'record_id', 'subject_id', 'y_true', 'pred_probability'}
    null_columns = {'model', 'scheme', 'null_scheme', 'permutation',
                    'row_index', 'y_permuted', 'pred_probability'}
    control_columns = {'model', 'scheme', 'null_scheme', 'permutation', 'original_split_seed',
                       'permutation_seed', 'label_vector_sha256', 'oof_auc'}
    require(control_columns <= set(controls), 'Missing reference-control columns')
    require(len(controls) == 792 and set(controls.original_split_seed) == {SEED},
            'Reference controls must use the frozen first split')
    keys = ['model', 'scheme', 'null_scheme', 'permutation']
    require(not controls.duplicated(keys).any(), 'Duplicate reference-control key')
    controls = controls.set_index(keys)
    require(real_columns <= set(real), 'Missing observed-prediction columns')
    require(null_columns <= set(null), 'Missing reference-prediction columns')
    observed = real[(real['repeat'] == 0) & (real['seed'] == SEED)]
    require(len(observed) == 4 * 195, 'Expected four complete fixed-first-split predictions')
    require(len(null) == 8 * N_REFERENCE * 195, 'Expected 792 complete reference predictions')
    require(set(observed.model) == MODELS and set(observed.scheme) == SPLITS,
            'Unexpected model or split identifiers')
    require(set(null.model) == MODELS and set(null.scheme) == SPLITS
            and set(null.null_scheme) == REFERENCES, 'Unexpected reference identifiers')
    truth = {}
    baseline = None
    expected_rows = np.arange(195)
    for (model, split), frame in observed.groupby(['model', 'scheme']):
        frame = frame.sort_values('row_index')
        require(np.array_equal(frame.row_index, expected_rows), 'Observed rows duplicate or incomplete')
        require(frame.record_id.is_unique, 'Duplicate record identifier')
        require(set(frame.y_true) == {0, 1}, 'Observed labels must be binary')
        identity = frame[['record_id', 'subject_id', 'y_true']].reset_index(drop=True)
        if baseline is None:
            baseline = identity
        require(identity.equals(baseline), 'Identity or labels differ across observed combinations')
        require(frame.subject_id.nunique() == 32 and frame.y_true.sum() == 147,
                'Unexpected E1 subject or class counts')
        require(frame.groupby('subject_id').y_true.nunique().eq(1).all(),
                'Original label is inconsistent within subject')
        ex = exact_auc(frame.y_true, frame.pred_probability)
        sk = float(roc_auc_score(frame.y_true, frame.pred_probability))
        require(abs(float(ex) - sk) < 1e-12, 'Independent AUROC implementations disagree')
        truth[model, split] = (ex, sk)
    require(len(truth) == 4, 'Missing observed combination')
    reference_rows = []
    vectors = {}
    for (model, split, kind, index), frame in null.groupby(
            ['model', 'scheme', 'null_scheme', 'permutation']):
        frame = frame.sort_values('row_index')
        require(np.array_equal(frame.row_index, expected_rows), 'Reference rows duplicate or incomplete')
        require(index in range(N_REFERENCE), 'Unexpected permutation index')
        require((model, split, kind, index) in controls.index, 'Missing reference-control key')
        control = controls.loc[(model, split, kind, index)]
        require(control.permutation_seed == 510001 + 10000 * (kind == 'record_label') + index,
                'Unexpected reference permutation seed')
        labels = frame.y_permuted.to_numpy()
        require(set(labels) == {0, 1}, 'Reference labels must be binary')
        label_hash = hashlib.sha256(np.asarray(labels, dtype=np.int64).tobytes()).hexdigest()
        require(label_hash == control.label_vector_sha256, 'Reference label hash mismatch')
        vk = kind, index
        if vk not in vectors:
            vectors[vk] = labels
        require(np.array_equal(labels, vectors[vk]), 'Reference labels differ across models or splits')
        if kind == 'subject_label':
            labeled = baseline.assign(permuted=labels)
            require(labeled.groupby('subject_id').permuted.nunique().eq(1).all(),
                    'Subject-label reference breaks within-subject label consistency')
            require(labeled.groupby('subject_id').permuted.first().sum() == 24,
                    'Subject-label reference changes subject case count')
        else:
            require(labels.sum() == 147, 'Row-label reference changes recording case count')
        exact = exact_auc(labels, frame.pred_probability)
        sk = float(roc_auc_score(labels, frame.pred_probability))
        observed_exact, observed_sk = truth[model, split]
        require(abs(float(exact) - sk) < 1e-12, 'Reference AUROC implementations disagree')
        require(np.isfinite(control.oof_auc) and abs(control.oof_auc - sk) < 1e-12,
                'Reference predictions disagree with saved control AUROC')
        ge = exact >= observed_exact
        require(ge == (sk >= observed_sk), 'Upper-tail comparison is sensitive to numerical rounding')
        reference_rows.append(dict(
            model=model, scheme=split, null_scheme=kind, permutation=int(index),
            original_split_seed=SEED, repeat=0,
            reference_auc=sk, reference_auc_exact_pairwise=float(exact),
            reference_auc_calculation_difference=abs(sk - float(exact)),
            observed_auc=observed_sk, reference_ge_observed=bool(ge),
            reference_eq_observed=bool(exact == observed_exact),
            analysis_status=SCOPE))
    refs = pd.DataFrame(reference_rows)
    require(len(refs) == 792 and len(vectors) == 198, 'Missing reference combination')
    summaries = []
    for (model, split, kind), group in refs.groupby(['model', 'scheme', 'null_scheme']):
        require(len(group) == N_REFERENCE and set(group.permutation) == set(range(N_REFERENCE)),
                'Reference index set incomplete')
        exact, sk = truth[model, split]
        ge = int(group.reference_ge_observed.sum())
        eq = int(group.reference_eq_observed.sum())
        rank = Fraction(1 + ge, N_REFERENCE + 1)
        summaries.append(dict(
            model=model, scheme=split, null_scheme=kind, original_split_seed=SEED, repeat=0,
            observed_auc=sk, observed_auc_exact_pairwise=float(exact),
            observed_auc_calculation_difference=abs(sk - float(exact)),
            n_reference=N_REFERENCE, n_reference_ge_observed=ge,
            n_reference_eq_observed=eq, monte_carlo_tail_rank=float(rank),
            descriptive_line=0.05, strictly_below_0_05=rank < Fraction(5, 100),
            on_0_05_line=rank == Fraction(5, 100),
            reference_mean=float(group.reference_auc.mean()),
            reference_median=float(group.reference_auc.median()),
            analysis_status=SCOPE))
    require(len(summaries) == 8, 'Missing summary combination')
    out.mkdir(parents=True, exist_ok=False)
    pd.DataFrame(summaries).to_csv(out / 'e1_exploratory_reference_comparisons.csv',
                                 index=False, float_format='%.17g')
    refs.to_csv(out / 'e1_exploratory_reference_statistics.csv',
                index=False, float_format='%.17g')
    return len(summaries), len(refs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    try:
        n_summary, n_ref = run(args.input_dir, args.out)
    except (ValueError, FileNotFoundError, OSError) as exc:
        parser.error(str(exc))
    print(f'Wrote {n_summary} descriptive comparisons and {n_ref} reference AUROCs.')
    print('No models refitted; ranks are not calibrated p values.')


if __name__ == '__main__':
    main()
