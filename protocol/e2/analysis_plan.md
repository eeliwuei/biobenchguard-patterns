# E2: frozen original-label CRC decision case

Version: 1; written before E2 pilot or main model fitting, 2026-09-08.

## Target and source

Primary question: Which fixed candidate, L2 logistic regression or a 120-tree random forest, should be prioritized for further validation on a **new independent CRC study**? The population is adults represented by the three available microbiomeHD CRC case/control studies (Zeller, Baxter, Zackular), not clinical screening deployment. Outcome is original metadata `DiseaseState`: CRC=1, H=0; other labels are excluded and counted. The statistical unit for cross-study generalization is a study (only three). All available eligible samples are used, without synthetic labels.

Duvallet et al. 2017, DOI 10.1038/s41467-017-01973-8, Methods/data-processing description supports RDP-assigned OTU genus aggregation, division by total sample reads, and discarding OTUs not assigned a genus. It does not validate this new analysis or provide a numeric diagnostic cutoff. Input provenance is archived by SHA256. The local source files came from the historical public microbiomeHD Zenodo record 840333; hashes are recorded separately from that bibliographic origin.

## Predictor construction and scope

Aggregate named genera only. Denominator is each sample's total count across all OTUs, before excluding unassigned taxa. Taxonomy union is used only as a lossless array storage schema; every fitting fold restricts columns to named genera with nonzero variance in its training samples. Thus held-out-only genera cannot enter a fitted model. Apply fixed `log1p(1e6 * relative_abundance)` sample-wise. This transform has no fitted parameters. No study ID, metadata, or artificial group embedding is a predictor.

The primary candidates use absolute Pearson association with training labels to retain k=20 train-variable genera (stable lexical feature tie-break). Train-only StandardScaler precedes logistic regression; the same scaling is used for RF to keep preparation consistent. Logistic regression: C=1, L2, solver=lbfgs, max_iter=3000. RF: 120 trees, max_features=sqrt, min_samples_leaf=1, unlimited depth, bootstrap=True. Both have **one fixed hyperparameter configuration, no search, identical folds and fit counts**. Equal fit counts are not a claim of equal computational cost; fit wall time is saved. The all-train-variable-genus variant is a prespecified strong baseline sensitivity for both model classes, preventing a weak selected linear baseline from being the only comparator.

No hyperparameter is selected from held-out cohorts; there is no tuning requiring an inner CV. Any learned centering, scaling, variance restriction and supervised selection for legitimate candidates is fitted on training data only.

## Primary evaluation and replay diagnostic

Leave one study out (LOSO), train on both others. Save all sample-level out-of-fold predictions, per-study AUROC, and an equal-study mean AUROC. Pooled OOF AUROC is a separate secondary quantity: it mixes between-study score comparability and is not the primary new-study estimand.

Within exactly the same fold, training cohort, model and fit budget, run four arms:

- A: k=20 highest training feature variance, without labels.
- B: k=20 highest absolute feature/outcome correlation using training samples only; this is the legitimate primary selected candidate.
- C: select the same k from the same **training-variable feature universe**, but calculate supervised association using training plus held-out labels. This is a deliberately nondeployable controlled replay of scope violation.
- D: all training-variable named genera, no supervised selection; legitimate strong baseline sensitivity.

The observed diagnostic is C minus B AUROC. It measures response to this replay intervention, not whether an unknown original pipeline leaked. B minus A and C minus A are also saved. C is never a valid candidate for selection or future deployment. A small or negative C-B does not show that the scope violation was absent.

## Controls, uncertainty and actions

Matched training-size control: for each LOSO held-out-study test size, create 10 repeated stratified random splits of all eligible samples with **exactly the same train/test counts**; evaluate both candidates' B and D arms. This is explicitly a second estimand, new samples from the mixture of observed studies. Random-to-LOSO differences are not represented as model reversal for a fixed target. Each sample appears in one test set within a random split; the repeated split results are not independent biological replicates.

Within each observed held-out study, obtain 2,000 paired stratified bootstrap resamples of its frozen OOF predictions for RF-minus-logistic and C-minus-B AUROC differences. These percentile intervals condition on these samples and fitted models; they do not quantify new-study sampling or model-refitting uncertainty. For the three-study mean, report each study and its range, not a narrow sample-level interval masquerading as cohort replication. No claims of generic FPR, power, or sensitivity are estimated from these real cases.

Predeclared illustrative candidate-prioritization policy (not a community standard): prioritize RF (or logistic symmetrically) only if its equal-study mean B-arm advantage is at least 0.02 AUROC, it has the same strict advantage in all three observed studies, and the D-arm mean advantage has the same sign. Otherwise retain both candidates for another independent study. Prioritization is an action for further validation, not proof of population-wide superiority. Report bootstrap intervals alongside the action without retroactively changing this rule.

For selection-scope response, the legacy fixed numerical threshold 0.02 is an **illustrative prespecified baseline**, with no claimed community-standard or FPR meaning. Exceedance means that the defined replay changed AUROC by more than 0.02; nonexceedance is not a clean verdict. A calibration will be joined only if feature construction, original-label law, model, selector, metric and LOSO protocol are applicable. The historical 195-feature family-fallback geometry and synthetic null do not satisfy this compatibility rule. If no compatible independently assessed calibration exists, report `calibration_not_established` and request target-matched calibration or independent outcome validation; do not assign calibrated present/absent from an unrelated curve.

Good fold-internal preprocessing and group-aware validation are valid conventional evidence. The new report must acknowledge unchanged actions or no demonstrated additional value when appropriate; a ranking reversal or claimed benefit is not required.

## Pilot, execution and outputs

Pilot: perform loader integrity checks, inspect class counts and IDs, and fit one B-arm fold per model with predictions withheld from output. Pilot checks finite predictions, two classes, no train/test overlap, training-only selected features and the no-heldout-label dependence of B. No pilot score is used to tune choices. After checks pass, record the protocol and implementation SHA256 and run the unchanged main analysis. A failed integrity check may be repaired with a dated protocol amendment before main fitting.

Outputs: input/feature manifest, pilot and execution logs, frozen protocol hash, actual OOF predictions, fold membership and selected features, train/test sizes and runtime, per-study and equal-study results, paired conditional intervals, matched-size controls, a decision table with required fields, and a mechanistic scope/limitations report. Files live only in `research_code/revision/e2`, `protocol/e2`, `results_raw/e2`, and `results_summary/e2`.

## Wheat boundary

The official BGLR manual (https://cran.r-project.org/web/packages/BGLR/BGLR.pdf, wheat sections) describes average grain yield in four mega-environments, with Y column 1 corresponding to environment 1; it does not identify per-line year or trial site here. The distributed marker matrix has prior full-data filtering/imputation. Wheat is therefore not promoted to a second independent deployment case in E2. It remains suitable for clearly bounded reference-geometry experiments elsewhere.

## Integrity amendment before any model fit

The first loader pilot stopped before fitting because four of 92 Zackular OTU sample columns have no matching row among the source's 90 metadata rows (88 overlap). Zeller: 129 of 129 OTU columns match; Baxter: 490 of 490 match. The loader now explicitly excludes and counts unmatched OTU samples through the same original-label eligibility gate (H/CRC required); it also reports metadata rows without an OTU column. No label is inferred for an unmatched sample. No model score was inspected, and model, target, selection, uncertainty and action rules are unchanged. This amendment precedes the successful integrity pilot and main-analysis freeze.
