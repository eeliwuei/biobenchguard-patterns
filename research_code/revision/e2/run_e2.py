#!/usr/bin/env python3
"""Original-label CRC study decision case; see frozen protocol/e2/analysis_plan.md."""
from __future__ import annotations
import argparse, hashlib, json, platform, time, warnings
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler

SEED = 2026090802
STUDIES = ["crc_zeller", "crc_baxter", "crc_zackular"]
MODELS = ["logistic_C1", "random_forest_120"]
ARMS = ["A_variance20", "B_train_supervised20", "C_full_supervised20", "D_all_train_variable"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def auc(y, scores):
    y = np.asarray(y); n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((rankdata(scores)[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def load_strict_crc(raw):
    frames, outcomes, group, sample_ids, audit, hashes = [], [], [], [], [], {}
    for study in STUDIES:
        base = Path(raw) / "microbiome/data" / f"{study}_results"
        otu = base / "RDP" / f"{study}.otu_table.100.denovo.rdp_assigned"
        metadata = base / f"{study}.metadata.txt"
        hashes[str(otu.relative_to(raw))] = sha(otu)
        hashes[str(metadata.relative_to(raw))] = sha(metadata)
        with otu.open() as f:
            samples = f.readline().rstrip("\n").split("\t")[1:]
            assert len(samples) == len(set(samples)), "Duplicate OTU sample IDs"
            counts = defaultdict(lambda: np.zeros(len(samples)))
            total = np.zeros(len(samples)); assigned = np.zeros(len(samples))
            discarded_otus = 0; all_otus = 0
            for line in f:
                cells = line.rstrip("\n").split("\t")
                assert len(cells) == len(samples) + 1
                values = np.asarray(cells[1:], dtype=float)
                assert np.isfinite(values).all() and (values >= 0).all()
                total += values; all_otus += 1
                genus = next((t[3:] for t in cells[0].split(";") if t.startswith("g__")), "")
                if genus:
                    counts["g__" + genus] += values; assigned += values
                else:
                    discarded_otus += 1
        assert (total > 0).all(), "Zero-read sample requires explicit handling"
        relative = pd.DataFrame(counts, index=samples).div(total, axis=0)
        md = pd.read_csv(metadata, sep="\t", dtype=str, encoding="latin-1")
        id_column = max(md.columns, key=lambda c: md[c].isin(set(samples)).sum())
        assert not md[id_column].duplicated().any()
        md = md.set_index(id_column)
        missing_metadata_n = len(set(samples) - set(md.index))
        label = md["DiseaseState"].reindex(samples)
        keep = label.isin(["CRC", "H"])
        relative = relative.loc[keep]
        y = (label.loc[keep] == "CRC").to_numpy(dtype=int)
        assert len(np.unique(y)) == 2
        frames.append(relative); outcomes.append(y)
        group += [study] * len(y)
        sample_ids += [study + ":" + s for s in relative.index]
        audit.append(dict(study=study, metadata_id_column=id_column,
                          n_otu_samples=len(samples), n_eligible=len(y), n_crc=int(y.sum()),
                          n_h=int((y == 0).sum()), n_excluded=int((~keep).sum()),
                          n_missing_metadata=missing_metadata_n,
                          n_metadata_rows_without_otu=len(set(md.index) - set(samples)),
                          excluded_labels=json.dumps(label.loc[~keep].value_counts(dropna=False).to_dict()),
                          total_otu_rows=all_otus, unassigned_otu_rows_discarded=discarded_otus,
                          named_genera_available=len(relative.columns),
                          assigned_read_fraction_median=float(np.median((assigned / total)[keep])),
                          min_total_reads=int(total[keep].min())))
    features = sorted(set.union(*(set(frame.columns) for frame in frames)))
    relative = np.vstack([frame.reindex(columns=features, fill_value=0).to_numpy() for frame in frames])
    X = np.log1p(relative * 1e6)
    y = np.concatenate(outcomes); g = np.asarray(group); ids = np.asarray(sample_ids)
    assert len(ids) == len(set(ids)) and np.isfinite(X).all()
    meta = dict(n=len(y), p=len(features), sample_ids=ids.tolist(), features=features,
                source_files_sha256=hashes, preprocessing="strict named-genus; original all-OTU per-sample denominator; log1p(1e6*x); union storage; train-variable restriction per fold",
                outcome="original DiseaseState H=0, CRC=1", study_counts=audit,
                source="Duvallet et al. 2017, DOI 10.1038/s41467-017-01973-8; microbiomeHD Zenodo 840333",
                differs_from_historical_reference_geometry=True)
    return X, y, g, ids, np.asarray(features), meta


def association(X, y):
    xc = X - X.mean(axis=0); yc = y - y.mean()
    den = np.sqrt((xc * xc).sum(axis=0) * (yc * yc).sum())
    return np.divide(np.abs(xc.T @ yc), den, out=np.zeros(X.shape[1]), where=den > 0)


def selected_features(X, y, train, test, arm):
    variable = np.flatnonzero(np.var(X[train], axis=0) > 1e-12)
    assert len(variable) >= 20
    if arm == "D_all_train_variable":
        return variable
    if arm == "A_variance20":
        scores = np.var(X[train][:, variable], axis=0)
    else:
        scope = np.concatenate([train, test]) if arm == "C_full_supervised20" else train
        scores = association(X[scope][:, variable], y[scope])
    # Array feature order is lexical; stable sort provides deterministic tie-break.
    return variable[np.argsort(-scores, kind="stable")[:20]]


def fit_predict(X, y, train, test, arm, model_name, seed):
    selected = selected_features(X, y, train, test, arm)
    scaler = StandardScaler().fit(X[train][:, selected])
    a = scaler.transform(X[train][:, selected]); b = scaler.transform(X[test][:, selected])
    if model_name == MODELS[0]:
        model = LogisticRegression(C=1.0, max_iter=3000, solver="lbfgs", random_state=seed)
    else:
        model = RandomForestClassifier(n_estimators=120, max_features="sqrt", min_samples_leaf=1,
                                       bootstrap=True, n_jobs=1, random_state=seed)
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        model.fit(a, y[train])
    runtime = time.perf_counter() - start
    pred = model.predict_proba(b)[:, 1]
    assert np.isfinite(pred).all() and ((pred >= 0) & (pred <= 1)).all()
    return pred, selected, runtime, [str(w.message) for w in ws]


def interval(y, difference_pairs, rng, B=2000):
    neg, pos = np.flatnonzero(y == 0), np.flatnonzero(y == 1)
    values = np.empty((B, len(difference_pairs)))
    for b in range(B):
        index = np.concatenate([rng.choice(neg, len(neg)), rng.choice(pos, len(pos))])
        for j, (p1, p0) in enumerate(difference_pairs):
            values[b, j] = auc(y[index], p1[index]) - auc(y[index], p0[index])
    return [(float(np.quantile(values[:, j], .025)), float(np.quantile(values[:, j], .975)))
            for j in range(len(difference_pairs))]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--pilot", action="store_true")
    args = ap.parse_args(); root = Path(args.root).resolve()
    rawout = root / "results_raw/e2"; summary = root / "results_summary/e2"; proto = root / "protocol/e2"
    for path in [rawout, summary, proto]: path.mkdir(parents=True, exist_ok=True)
    npz, metap = rawout / "crc_strict_inputs.npz", rawout / "input_manifest.json"
    if npz.exists():
        data = np.load(npz); X, y, g, ids, features = [data[k] for k in ["X", "y", "g", "ids", "features"]]
        meta = json.loads(metap.read_text())
    else:
        X, y, g, ids, features, meta = load_strict_crc(root / "data_manifest/raw")
        np.savez_compressed(npz, X=X, y=y, g=g, ids=ids, features=features)
        metap.write_text(json.dumps(meta, indent=2))
        pd.DataFrame(meta["study_counts"]).to_csv(summary / "cohort_input_audit.csv", index=False)
    folds = [(st, np.flatnonzero(g != st), np.flatnonzero(g == st)) for st in STUDIES]
    if args.pilot:
        st, train, test = folds[0]; checks = {}
        for model in MODELS:
            pred, selected, runtime, ws = fit_predict(X, y, train, test, ARMS[1], model, SEED)
            swapped = y.copy(); swapped[test] = 1 - swapped[test]
            pred2, selected2, _, _ = fit_predict(X, swapped, train, test, ARMS[1], model, SEED)
            checks[model] = dict(finite_predictions=bool(np.isfinite(pred).all()),
                no_train_test_overlap=len(set(train) & set(test)) == 0,
                B_invariant_to_heldout_labels=bool(np.array_equal(selected, selected2) and np.array_equal(pred, pred2)),
                selected_k=len(selected), fitting_warnings=ws, runtime_seconds=runtime)
            assert checks[model]["B_invariant_to_heldout_labels"] and not ws
        record = dict(phase="pilot", timestamp=datetime.now(timezone.utc).isoformat(), scores_inspected=False,
                      plan_sha256=sha(proto / "analysis_plan.md"), script_sha256=sha(__file__),
                      n=len(y), p=X.shape[1], checks=checks)
        (rawout / "pilot_integrity.json").write_text(json.dumps(record, indent=2))
        print(json.dumps(record, indent=2), flush=True)
        return
    pilot = json.loads((rawout / "pilot_integrity.json").read_text())
    assert pilot["plan_sha256"] == sha(proto / "analysis_plan.md")
    assert pilot["script_sha256"] == sha(__file__), "Code changed after pilot; rerun integrity pilot"
    freeze = dict(timestamp=datetime.now(timezone.utc).isoformat(), plan_sha256=sha(proto / "analysis_plan.md"),
                  script_sha256=sha(__file__), input_npz_sha256=sha(npz), phase="before_main_fits",
                  seed=SEED, python=platform.python_version(), numpy=np.__version__,
                  pandas=pd.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__)
    (proto / "execution_freeze.json").write_text(json.dumps(freeze, indent=2))
    predictions, memberships, fitted, fold_metrics, selections = [], [], [], [], []
    store = {}
    for j, (study, train, test) in enumerate(folds):
        for role, idx in [("train", train), ("test", test)]:
            memberships += [dict(protocol="LOSO", split_id=study, sample_id=ids[i], study=g[i], role=role) for i in idx]
        for model in MODELS:
            for arm in ARMS:
                pred, selected, elapsed, ws = fit_predict(X, y, train, test, arm, model, SEED + j)
                store[(study, model, arm)] = pred
                fitted.append(dict(protocol="LOSO", split_id=study, model=model, arm=arm, train_n=len(train),
                                   test_n=len(test), n_features=len(selected), runtime_seconds=elapsed, warnings=json.dumps(ws)))
                for idx in selected:
                    selections.append(dict(protocol="LOSO", split_id=study, model=model, arm=arm, feature=features[idx]))
                fold_metrics.append(dict(study=study, model=model, arm=arm, n_train=len(train), n_test=len(test),
                                         n_crc=int(y[test].sum()), n_h=int((y[test] == 0).sum()), auroc=auc(y[test], pred)))
                predictions += [dict(protocol="LOSO", split_id=study, sample_id=ids[i], study=g[i], y=int(y[i]),
                                     model=model, arm=arm, prediction=float(p)) for i, p in zip(test, pred)]
        print("Finished LOSO", study, flush=True)
    fm = pd.DataFrame(fold_metrics); fm.to_csv(summary / "loso_per_study.csv", index=False)
    pd.DataFrame(predictions).to_csv(rawout / "loso_predictions.csv", index=False)
    pd.DataFrame(selections).to_csv(rawout / "selected_features.csv", index=False)
    totals = []
    for model in MODELS:
        for arm in ARMS:
            rows = fm[(fm.model == model) & (fm.arm == arm)]
            yy = np.concatenate([y[test] for _, _, test in folds])
            pp = np.concatenate([store[(study, model, arm)] for study, _, _ in folds])
            totals.append(dict(model=model, arm=arm, study_count=3, equal_study_mean_auroc=rows.auroc.mean(),
                               min_study_auroc=rows.auroc.min(), max_study_auroc=rows.auroc.max(),
                               pooled_oof_auroc=auc(yy, pp), n=len(yy)))
    pd.DataFrame(totals).to_csv(summary / "loso_summary.csv", index=False)
    contrasts = []
    for j, (study, train, test) in enumerate(folds):
        definitions = [dict(contrast="RF_minus_logistic_B", model="pair", arm_plus=ARMS[1], arm_minus=ARMS[1],
                            plus=store[(study, MODELS[1], ARMS[1])], minus=store[(study, MODELS[0], ARMS[1])])]
        for model in MODELS:
            for a, b, name in [(ARMS[2], ARMS[1], "C_minus_B"), (ARMS[1], ARMS[0], "B_minus_A"), (ARMS[2], ARMS[0], "C_minus_A")]:
                definitions.append(dict(contrast=name, model=model, arm_plus=a, arm_minus=b,
                                        plus=store[(study, model, a)], minus=store[(study, model, b)]))
        cis = interval(y[test], [(d["plus"], d["minus"]) for d in definitions], np.random.default_rng(SEED + 200 + j))
        for d, (lo, hi) in zip(definitions, cis):
            contrasts.append(dict(study=study, contrast=d["contrast"], model=d["model"],
                estimate=auc(y[test], d["plus"]) - auc(y[test], d["minus"]), lower=lo, upper=hi,
                interval="paired stratified sample bootstrap; fixed fitted OOF predictions; conditional within this study",
                bootstrap_B=2000, group_sampling_uncertainty="not estimated; only 3 observed studies"))
    ct = pd.DataFrame(contrasts); ct.to_csv(summary / "per_study_contrasts.csv", index=False)
    random_preds, random_scores = [], []
    for j, (target_study, _, target_test) in enumerate(folds):
        splitter = StratifiedShuffleSplit(n_splits=10, test_size=len(target_test), random_state=SEED + 1000 + j)
        for r, (train, test) in enumerate(splitter.split(X, y)):
            sid = f"matched_{target_study}_repeat{r:02d}"
            assert len(train) == len(y) - len(target_test)
            for role, idx in [("train", train), ("test", test)]:
                memberships += [dict(protocol="random_matched_size", split_id=sid, sample_id=ids[i], study=g[i], role=role) for i in idx]
            for model in MODELS:
                for arm in [ARMS[1], ARMS[3]]:
                    pred, selected, elapsed, ws = fit_predict(X, y, train, test, arm, model, SEED + 1000 + 100*j + r)
                    fitted.append(dict(protocol="random_matched_size", split_id=sid, model=model, arm=arm,
                        train_n=len(train), test_n=len(test), n_features=len(selected), runtime_seconds=elapsed, warnings=json.dumps(ws)))
                    random_scores.append(dict(matched_loso_study=target_study, repeat=r, model=model, arm=arm,
                                              n_train=len(train), n_test=len(test), auroc=auc(y[test], pred)))
                    random_preds += [dict(protocol="random_matched_size", split_id=sid, sample_id=ids[i], study=g[i],
                        y=int(y[i]), model=model, arm=arm, prediction=float(p)) for i, p in zip(test, pred)]
        print("Finished matched-size control", target_study, flush=True)
    rs = pd.DataFrame(random_scores); rs.to_csv(summary / "random_matched_size_scores.csv", index=False)
    pd.DataFrame(random_preds).to_csv(rawout / "random_matched_size_predictions.csv", index=False)
    pd.DataFrame(memberships).to_csv(rawout / "fold_membership.csv", index=False)
    fits = pd.DataFrame(fitted); fits.to_csv(rawout / "fit_execution.csv", index=False)
    bridge = []
    for (study, model, arm), rows in rs.groupby(["matched_loso_study", "model", "arm"]):
        observed = fm[(fm.study == study) & (fm.model == model) & (fm.arm == arm)].iloc[0]
        bridge.append(dict(matched_loso_study=study, model=model, arm=arm, train_n=int(observed.n_train),
                           test_n=int(observed.n_test), loso_auroc=observed.auroc,
                           random_matched_mean_auroc=rows.auroc.mean(), random_split_sd=rows.auroc.std(),
                           random_minus_loso=rows.auroc.mean()-observed.auroc, random_repeats=10,
                           estimand_note="random mixture new-sample target differs from LOSO new-study target"))
    pd.DataFrame(bridge).to_csv(summary / "matched_train_size_comparison.csv", index=False)
    b = fm[fm.arm == ARMS[1]].pivot(index="study", columns="model", values="auroc")
    d = fm[fm.arm == ARMS[3]].pivot(index="study", columns="model", values="auroc")
    delta = b[MODELS[1]] - b[MODELS[0]]; delta_d = (d[MODELS[1]] - d[MODELS[0]]).mean()
    if delta.mean() >= .02 and (delta > 0).all() and delta_d > 0:
        action = "Prioritize RF120 for another independent-study validation; retain logistic comparator"
    elif delta.mean() <= -.02 and (delta < 0).all() and delta_d < 0:
        action = "Prioritize logistic for another independent-study validation; retain RF120 comparator"
    else:
        action = "Retain both candidates for another independent-study validation; no preference under frozen policy"
    limits = "Three observed studies; heterogeneous protocols/case definitions; fixed-model comparison; no prospective cohort; no clinical utility or generic FPR claim"
    decisions = [dict(case="CRC_LOSO_candidate_priority", fixed_target="AUROC in a new independent CRC study",
        model_pair="logistic C1 versus RF120; clean training-only k20 selection",
        conventional_evidence=f"LOSO equal-study B AUROC logistic={b[MODELS[0]].mean():.6f}, RF={b[MODELS[1]].mean():.6f}; RF-LR by study={delta.to_dict()}; all-feature mean RF-LR={delta_d:.6f}",
        conventional_decision=action,
        calibrated_diagnostic_evidence="No independently validated task-matched calibration: strict named-genus/log transform, original outcome law and LOSO differ from reference calibration. calibration_not_established.",
        revised_decision=action,
        confirmatory_evidence="All-feature D comparison and same-size random splits are internal sensitivity checks. No new external confirmatory study was available.",
        action=action, limitations=limits + "; no demonstrated additional candidate-selection benefit over proper LOSO validation")]
    diagnostics = []
    for model in MODELS:
        subset = ct[(ct.model == model) & (ct.contrast == "C_minus_B")]
        value = float(subset.estimate.mean()); status = "exceeds_illustrative_0.02" if value > .02 else "does_not_exceed_illustrative_0.02"
        diagnostics.append(dict(case="CRC_LOSO_"+model, statistic="equal-study mean AUROC(C)-AUROC(B)",
            observed=value, illustrative_threshold=.02, baseline_response=status,
            calibration_status="calibration_not_established", historical_geometry_compatible=False,
            original_pipeline_leakage_status="not identifiable from C-B; B/C scope known by construction",
            action="Keep fold-internal B; C is an invalid deployment workflow; obtain compatible calibration before rate claims"))
        decisions.append(dict(case="CRC_LOSO_scope_response_"+model, fixed_target="AUROC in a new independent CRC study",
            model_pair=model+": train-only B versus nondeployable full-label C",
            conventional_evidence=f"Code-scope audit verifies B train-only and C heldout-label access; observed mean C-B={value:.6f}",
            conventional_decision="Retain B scope and reject C scope based on known information access, regardless of numeric response",
            calibrated_diagnostic_evidence=f"Illustrative 0.02 rule: {status}; no compatible original-outcome calibration; per-study paired intervals are conditional performance uncertainty, not detector operating rates",
            revised_decision="Known scope violation persists even when C-B is small or negative; numerical nonalert supplies no clean verdict",
            confirmatory_evidence="Same fitted protocol, raw labels, saved OOF predictions and exact feature selections; no independent natural negative control",
            action="Use B or no-selection D; seek another independent cohort and target-matched calibration before calibrated absence claims",
            limitations=limits+"; controlled replay on real labels is not a reconstruction of unknown historical pipeline actions"))
    pd.DataFrame(decisions).to_csv(summary / "real_benchmark_decision_table.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(summary / "actual_case_diagnostic_reports.csv", index=False)
    report = dict(frozen_policy_action=action, primary_equal_study_auroc={m:float(b[m].mean()) for m in MODELS},
                  rf_minus_logistic_B_by_study=delta.to_dict(), rf_minus_logistic_B_mean=float(delta.mean()),
                  rf_minus_logistic_D_mean=float(delta_d), diagnostics=diagnostics,
                  n=len(y), p=X.shape[1], studies=3, fit_counts=fits.groupby("model").size().to_dict(),
                  runtime_seconds_by_model=fits.groupby("model").runtime_seconds.sum().to_dict(),
                  warnings_count=int((fits.warnings!="[]").sum()), timestamp=datetime.now(timezone.utc).isoformat(),
                  limitations=limits, completion="main models actually fitted and predictions saved")
    (summary / "e2_results.json").write_text(json.dumps(report, indent=2))
    (summary / "E2_REPORT.md").write_text("# Original-label CRC decision case\n\n" +
        f"Action under the protocol frozen before fitting: **{action}**.\n\n" +
        f"{len(y)} original-label samples from three studies; {X.shape[1]} named-genus storage columns, with train-only variable-column restriction. Equal-study LOSO B-arm AUROC: logistic {b[MODELS[0]].mean():.3f}; RF120 {b[MODELS[1]].mean():.3f}. The individual-study RF-minus-logistic differences are {delta.to_dict()}.\n\n" +
        "The data-defined real cases were actually fitted. Prediction files, fold membership, selected taxa, counts, runtimes and input hashes are saved. Both model families have one fixed configuration and equal fit counts; CPU cost is recorded separately. The all-feature D variants and exact training-size random controls are reported separately. Random validation addresses a different target from a new study.\n\n" +
        "The decision table acknowledges that sound LOSO validation already supports the reported action; this case does not establish an additional candidate-selection benefit from calibrated diagnostics. No new independent prospective cohort was obtained.\n\n" +
        "## Selection-scope replay\n\n" + "\n".join(f"- {x['case']}: mean C-B={x['observed']:.6f}; {x['baseline_response']}; calibration_not_established." for x in diagnostics) +
        "\n\nC deliberately accesses held-out outcomes during selection and is ineligible for deployment even if its advantage is small or negative. C-B measures response to that replay; it cannot identify the processing history of an unknown pipeline. The legacy 0.02 threshold is an illustrative numeric baseline, with no FPR calibration.\n\n" +
        "Historical synthetic-reference thresholds are not transported: feature handling, label law and LOSO protocol differ. This honest refusal is part of the operating report, not evidence of a true-negative detector. Neither three real studies nor repeated splits provide generic detector FPR or power.\n\n" +
        "## Uncertainty and source\n\nPaired stratified bootstrap intervals condition on each observed study's frozen predictions. They exclude refitting and new-study sampling uncertainty. Three studies are insufficient for a precise claim about the distribution of future cohorts.\n\n" +
        "The input source is [Duvallet et al. 2017](https://www.nature.com/articles/s41467-017-01973-8), microbiomeHD Zenodo record 840333. This revision follows its named-genus exclusion rule but uses a new fixed log transform and models; it does not claim to reproduce the original paper's exact numerical experiment. H is the source harmonized control label, not a guarantee that every control is free of other disorders or confounders.\n\n" +
        "Wheat is not promoted to an independent E2 deployment case: [the BGLR manual](https://cran.r-project.org/web/packages/BGLR/BGLR.pdf) specifies yield averages in four mega-environments but provides no per-line year/site here, and distributed markers were previously filtered/imputed.\n")
    (rawout / "completion.json").write_text(json.dumps(dict(status="complete", report_sha256=sha(summary / "e2_results.json"),
        main_fit_count=len(fits), prediction_rows=len(predictions)+len(random_preds), freeze_sha256=sha(proto / "execution_freeze.json")), indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
