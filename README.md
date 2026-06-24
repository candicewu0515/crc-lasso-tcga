# crc-lasso-tcga

A reproducible LASSO gene-signature pipeline for classifying colorectal
adenocarcinoma (COAD) versus normal colon tissue from TCGA RNA-seq data.

Companion code for the ICBB 2026 paper *"LASSO-Based Gene Signature
Identification for Colorectal Cancer Classification Using TCGA Expression Data"*
(Xia Wu, University of Iowa).

## Overview

The pipeline combines:

1. **Differential expression** — `limma` empirical-Bayes moderated t-tests
   (tumor vs. normal; |log2FC| > 2, BH-adjusted *p* < 0.05).
2. **Unsupervised analysis** — PCA on all genes vs. the differentially
   expressed genes (DEGs).
3. **Feature selection** — L1-regularized (LASSO) logistic regression with
   10-fold cross-validation (`glmnet`), reporting `lambda.min` and
   `lambda.1se` signatures and out-of-fold (pre-validated) ROC/AUC.
4. **Functional enrichment** — GO biological-process and KEGG
   over-representation of the up- and down-regulated DEGs (`clusterProfiler`).
5. **Harder endpoint (MSI status)** — predicting microsatellite-instability
   status (MSI-H vs MSS) by **nested** 10-fold CV, with DEG selection done
   inside each fold to avoid leakage (`msi_prediction.R`).

Tumor / normal labels are taken **directly from the TCGA barcode sample-type
code** (`-01` = primary tumor, `-11` = solid tissue normal), giving 286 tumor
and 41 normal samples.

## Data

Expression data are **not** redistributed here. Download the TCGA COAD cohort
from [UCSC Xena](https://xenabrowser.net/):

- Dataset: `TCGA.COAD.sampleMap/HiSeqV2` (gene-level, log2(norm_count+1),
  20,530 genes × 329 samples)

```bash
mkdir -p data
curl -L "https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/HiSeqV2.gz" \
  | gunzip > data/HiSeqV2
```

MSI labels for `msi_prediction.R` come from the cBioPortal PanCancer Atlas
(MANTIS score; MSI-H = score ≥ 0.4):

```bash
curl -s "https://www.cbioportal.org/api/studies/coadread_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=SAMPLE&attributeId=MSI_SCORE_MANTIS&projection=DETAILED" \
  | python3 -c "import json,sys; d=json.load(sys.stdin); print('sampleId\tMSI_SCORE_MANTIS'); [print(x['sampleId'],x['value'],sep='\t') for x in d if x['value'] not in ('','NA')]" \
  > data/MSI_SCORE_MANTIS.tsv
```

## Run

```bash
Rscript crc_lasso_analysis.R   # tumor vs normal + DEG + enrichment
Rscript msi_prediction.R       # harder endpoint: MSI-H vs MSS (nested CV)
```

Requires R (≥4.0) with `limma`, `glmnet`, `ggplot2`, `pheatmap`, `gridExtra`,
`clusterProfiler` and `org.Hs.eg.db`.
A fixed seed (`set.seed(2026)`) makes the cross-validation reproducible.
Outputs are written to `figures/` and `output/`.

## Results (summary)

| Step | Result |
|------|--------|
| DEGs | 1,618 (358 up, 1,260 down) |
| PCA PC1 variance | 14.5% (all genes) → 40.6% (DEGs) |
| Signature | 16 genes (`lambda.min`) / 15 genes (`lambda.1se`) |
| Cross-validated AUC (tumor vs normal) | 1.00 |
| **MSI-H vs MSS, nested-CV AUC** | **0.928** (MLH1 + immune signature) |

Signature genes include up-regulated *CDH3, KRT80, ETV4, ESM1, FOXQ1* and
down-regulated candidate suppressors *OTOP2, CDH10*. Enrichment of the
up-regulated genes recovers canonical CRC pathways (Wnt, cadherin, Hippo
signaling); down-regulated genes reflect loss of normal colonic transport and
metabolic functions.

> **Note on the AUC.** Tumor-vs-normal separation in bulk RNA-seq is an
> intrinsically easy task, so the cross-validated AUC of 1.00 is expected and
> serves as a proof-of-concept sanity check rather than a claim of
> clinical-grade performance. The MSI-H vs MSS task is the harder, clinically
> meaningful benchmark: under nested CV (no feature-selection leakage) it
> reaches a realistic AUC of 0.928, with `MLH1` and interferon/antigen-
> presentation genes leading the signature — consistent with the hypermutated,
> immune-infiltrated MSI-H phenotype.

## License

MIT — see [LICENSE](LICENSE).
