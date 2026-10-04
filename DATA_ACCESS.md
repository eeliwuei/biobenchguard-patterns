# Source-data access and identity

Raw datasets and reconstructed sample-level matrices are not bundled. Use the cited primary sources and retain their conditions. Compare acquired files to `RAW_INPUT_SHA256.csv` before fitting; an accessible newer dataset is not automatically the analyzed version. Directory layout and preprocessing commands are in [the reproduction guide](reproduction/README_REPRODUCTION.md).

| Source | Acquisition route | Identity and terms |
|---|---|---|
| Parkinsons | [UCI dataset 174](https://archive.ics.uci.edu/dataset/174/parkinsons) | `parkinsons.data`,195 rows, 32 parsed subjects; CC BY 4.0. |
| BGLR wheat | [BGLR source object](https://raw.githubusercontent.com/gdlc/BGLR-R/master/data/wheat.RData), [primary documentation](https://cran.r-project.org/web/packages/BGLR/BGLR.pdf) | Fresh retrieval on 8 September 2026 reproduced SHA256 `8710523389007dd8b7c88a5deb74ee16aea0bae158496bc11f0626fb43b5ea21`; 477,712 bytes. Source labels 1, 2, 4, 5 are all retained in E5. The package license does not independently establish a right to redistribute the dataset. |
| MicrobiomeHD | [Zenodo record 840333](https://zenodo.org/records/840333) and `research_code/revision/e2/fetch_public_crc.py` | Three original CRC archives; fetcher verifies published MD5 values. Dataset CC BY-NC 4.0. Cite Duvallet2017 and each original cohort below. |
| Rice diversity panel | [Zhao et al. primary article](https://doi.org/10.1038/ncomms1467) | Obtain the original 413-line PLINK genotype file and the matching germplasm/population metadata; exact hashes are listed. Fresh downloads from the pinned PIQUE mirror and the archived original germplasm resource reproduced both analyzed hashes. The original-domain direct HTTPS route is unavailable; use the pinned routes below. No standalone redistribution permission is assumed. |

Original cohort attribution: [Zeller2014](https://doi.org/10.15252/msb.20145645), [Baxter2016](https://doi.org/10.1186/s13073-016-0290-3), [Zackular2014](https://doi.org/10.1158/1940-6207.CAPR-14-0129), and the [microbiomeHD collection paper](https://doi.org/10.1038/s41467-017-01973-8).

The derived source tables contain diagnostic statistics, not a relicensing of the source measurements. Public reuse, particularly commercial use of microbiomeHD, remains governed by the source terms.

## Exact input acquisition commands

```sh
.venv/bin/python research_code/revision/fetch_reference_inputs.py --raw-root raw
.venv/bin/python research_code/revision/e2/fetch_public_crc.py --destination raw
```

The first command downloads wheat, rice genotype, rice population metadata and UCI inputs, verifies archive/member SHA256 values, refuses inconsistent existing inputs and writes an acquisition receipt. It extracts only the explicitly named member bytes without extracting arbitrary archive paths. The second obtains the three original CRC archives and verifies their source MD5 values. Together they populate the raw-file layout described in the reproduction guide; pass `raw/e1/parkinsons.data` as the E1 input.

The rice genotype mirror is pinned to [PIQUE commit 097a945](https://raw.githubusercontent.com/tony-travis/PIQUE/097a945f1a715078a9f877469ea60f7497f7375d/test/RiceDiversity_44K_Genotypes_PLINK.tar.xz). Its tar.xz SHA256 is `0bdd54e44a831b197fcf0bdf5d1cec7055c721cb5ab94d7c33a005e889ab62a1`; the extracted PED member matches the input manifest. The population metadata comes from [the archived original RiceDiversity CSV](https://web.archive.org/web/20180408144810id_/http://ricediversity.org:80/data/sets/44kgwas/RiceDiversity.44K.germplasm.csv), 39,229 bytes. Its SHA256 is `5cfc7bdcdc755581366eb7c66ce0f576a9013e04fa6884aed566cd18b9a30df6`. These are source-identity checks; neither a mirror nor successful retrieval creates a new data license.
