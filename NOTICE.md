# Code and data attribution

## Licenses

The analysis code is distributed with the inherited Apache License 2.0 text in `LICENSE`, within the BioBenchGuard research-code lineage. The same license covers every other file outside `results_summary/` except the license texts: the test fixtures, the protocols, the execution and reproduction records, the manifests and the documentation.

The derived numerical tables under `results_summary/` (the tables of Data S1) are distributed under the Creative Commons Attribution-NonCommercial 4.0 International license (CC BY-NC 4.0). Copyright 2026 DeepFarm Inc. Reuse requires credit to the article named in `README.md` and is permitted for noncommercial purposes only; the terms are those of `LICENSE-DATA`.

Third-party source datasets keep their own terms. Neither license above applies to the third-party source datasets or to third-party analysis software.

## Contents

This package contains analysis source, recorded scientific protocols, numerical dependency pins and derived figure/source tables. Raw biological inputs and the original publication authors' R implementation are not included. Dataset provenance, attribution and permitted access are described in `DATA_ACCESS.md` and `reproduction/README_REPRODUCTION.md`; numerical replication and portability are distinguished in `reproduction/SCIENTIFIC_REPRODUCTION_REPORT.md`.

## Source datasets

The derived tables were computed from the following public datasets. Cite the sources when the tables are reused.

| Dataset | Source and citation | Source terms |
|---|---|---|
| Parkinson voice recordings (E1) | Little, M. (2007). Parkinsons [Dataset]. UCI Machine Learning Repository. https://doi.org/10.24432/C59C74 | CC BY 4.0 |
| Wheat lines and markers (E3–E5) | CIMMYT data distributed as the object `wheat` of the BGLR R package: Pérez, P., and de los Campos, G. (2014). Genetics 198, 483–495. https://doi.org/10.1534/genetics.114.164442 | Terms of the BGLR distribution; no separate redistribution right is assumed |
| Rice diversity panel (E3) | Zhao, K., et al. (2011). Nature Communications 2, 467. https://doi.org/10.1038/ncomms1467; Rice Diversity Project | Terms of the host; no separate redistribution right is assumed |
| Colorectal-cancer microbiome studies (E2–E4) | microbiomeHD: Duvallet, C., et al. (2017). Nature Communications 8, 1784. https://doi.org/10.1038/s41467-017-01973-8; Zenodo record 840333, https://doi.org/10.5281/zenodo.840333. Cohorts: Zeller et al. (2014), https://doi.org/10.15252/msb.20145645; Baxter et al. (2016), https://doi.org/10.1186/s13073-016-0290-3; Zackular et al. (2014), https://doi.org/10.1158/1940-6207.CAPR-14-0129 | CC BY-NC 4.0 |

The E1 analysis design follows Chaibub Neto, E., et al. (2019). npj Digital Medicine 2, 99. https://doi.org/10.1038/s41746-019-0178-x.
