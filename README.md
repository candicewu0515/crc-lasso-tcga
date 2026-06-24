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

## Run

```bash
Rscript crc_lasso_analysis.R
```

Requires R (≥4.0) with `limma`, `glmnet`, `ggplot2`, `pheatmap`, `gridExtra`.
A fixed seed (`set.seed(2026)`) makes the cross-validation reproducible.
Outputs are written to `figures/` and `output/`.

## Results (summary)

| Step | Result |
|------|--------|
| DEGs | 1,618 (358 up, 1,260 down) |
| PCA PC1 variance | 14.5% (all genes) → 40.6% (DEGs) |
| Signature | 16 genes (`lambda.min`) / 15 genes (`lambda.1se`) |
| Cross-validated AUC | 1.00 |

Signature genes include up-regulated *CDH3, KRT80, ETV4, ESM1, FOXQ1* and
down-regulated candidate suppressors *OTOP2, CDH10*. Enrichment of the
up-regulated genes recovers canonical CRC pathways (Wnt, cadherin, Hippo
signaling); down-regulated genes reflect loss of normal colonic transport and
metabolic functions.

> **Note on the AUC.** Tumor-vs-normal separation in bulk RNA-seq is an
> intrinsically easy task, so the cross-validated AUC of 1.00 is expected and
> serves as a proof-of-concept sanity check rather than a claim of
> clinical-grade performance. Harder endpoints (subtyping, staging, MSI) and
> external validation are the natural next steps.

## License

MIT — see [LICENSE](LICENSE).
