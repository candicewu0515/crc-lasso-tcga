# =============================================================================
# Harder endpoint: predicting microsatellite-instability (MSI) status in TCGA
# COAD from expression, using the same Random-Forest / Boruta pipeline.
# MSI labels: cBioPortal PanCancer Atlas MANTIS score (MSI-H = MANTIS >= 0.4).
# Honest performance via NESTED 10-fold CV: DEG selection + Boruta + RF are all
# fit inside each training fold only (no information leakage to the test fold).
# Author: Xia Wu
# =============================================================================
set.seed(2026)
suppressMessages({ library(limma); library(ggplot2); library(Boruta); library(ranger) })
dir.create("figures", showWarnings = FALSE); dir.create("output", showWarnings = FALSE)

# ---- Load expression + MSI labels ------------------------------------------
raw  <- read.delim("data/HiSeqV2", check.names = FALSE)
expr <- as.matrix(raw[, -1]); rownames(expr) <- raw[[1]]
expr <- expr[apply(expr, 1, function(x) var(x) > 0), ]
orig_rn <- rownames(expr)                       # original gene symbols
rownames(expr) <- make.names(orig_rn)           # syntactically valid (for Boruta)
gmap <- setNames(orig_rn, rownames(expr))       # clean -> original symbol

mantis <- read.delim("data/MSI_SCORE_MANTIS.tsv", check.names = FALSE)
lab <- setNames(mantis$MSI_SCORE_MANTIS, mantis$sampleId)
common <- intersect(colnames(expr), names(lab))
expr <- expr[, common]
msi  <- factor(ifelse(lab[common] >= 0.4, "MSI_H", "MSS"), levels = c("MSS", "MSI_H"))
y <- as.numeric(msi == "MSI_H")
cat(sprintf("Samples: %d (MSI-H=%d, MSS=%d)\n", length(y), sum(y), sum(1 - y)))

auc_fun <- function(s, l){ r <- rank(s); np <- sum(l==1); nn <- sum(l==0)
  (sum(r[l==1]) - np*(np+1)/2)/(np*nn) }
roc_pts <- function(s,l){ o <- order(s, decreasing=TRUE)
  data.frame(fpr=c(0,cumsum(l[o]==0)/sum(l==0)), tpr=c(0,cumsum(l[o]==1)/sum(l==1))) }
deg_select <- function(mat, grp){
  d <- model.matrix(~0 + grp); colnames(d) <- levels(grp)
  f <- eBayes(contrasts.fit(lmFit(mat, d), makeContrasts(MSI_H - MSS, levels = d)))
  tt <- topTable(f, coef = 1, number = Inf)
  g <- rownames(tt)[abs(tt$logFC) > 1 & tt$adj.P.Val < 0.05]
  if (length(g) < 20) g <- rownames(tt)[order(tt$P.Value)][1:50]
  g
}
boruta_sel <- function(Xtr, grp){
  b <- TentativeRoughFix(Boruta(x = Xtr, y = grp, maxRuns = 60, doTrace = 0))
  f <- getSelectedAttributes(b, withTentative = FALSE)
  if (length(f) < 5) f <- colnames(Xtr)              # fallback: all DEGs
  f
}

# ---- NESTED 10-fold CV: DEG + Boruta + RF inside each training fold ---------
K <- 10
folds <- integer(length(y))
for (cl in 0:1){ ix <- which(y == cl); folds[ix] <- sample(rep(1:K, length.out = length(ix))) }
oof <- rep(NA, length(y)); nsel <- c()
for (k in 1:K){
  tr <- folds != k; te <- folds == k
  deg  <- deg_select(expr[, tr], droplevels(msi[tr]))
  Xtr  <- t(expr[deg, tr]); Xte <- t(expr[deg, te])
  set.seed(2026 + k)
  feat <- boruta_sel(Xtr, droplevels(msi[tr]))
  nsel <- c(nsel, length(feat))
  rf <- ranger(x = Xtr[, feat, drop = FALSE], y = droplevels(msi[tr]),
               probability = TRUE, num.trees = 500, seed = 2026)
  oof[te] <- predict(rf, Xte[, feat, drop = FALSE])$predictions[, "MSI_H"]
}
auc_nested <- auc_fun(oof, y)
cat(sprintf("Nested-CV features/fold: median=%d  |  HONEST nested-CV AUC = %.3f\n",
            as.integer(median(nsel)), auc_nested))

# ---- Final signature on full data (for reporting) --------------------------
deg_all <- deg_select(expr, msi)
Xall <- t(expr[deg_all, ])
set.seed(2026)
bA <- TentativeRoughFix(Boruta(x = Xall, y = msi, maxRuns = 100, doTrace = 0))
impd <- attStats(bA); impd <- impd[impd$decision == "Confirmed", ]
impd <- impd[order(-impd$meanImp), ]
sig_clean <- head(rownames(impd), 22)
sig_genes <- unname(gmap[sig_clean])
cat(sprintf("Boruta confirmed %d genes; top MSI signature (%d): %s\n",
            nrow(impd), length(sig_genes), paste(sig_genes, collapse = ", ")))
write.csv(data.frame(gene = unname(gmap[rownames(impd)]), impd, row.names = NULL),
          "output/msi_boruta_confirmed.csv", row.names = FALSE)

# ---- Figure: ROC + predicted-probability boxplot (Fig 7) -------------------
rp <- roc_pts(oof, y)
p_roc <- ggplot(rp, aes(fpr, tpr)) + geom_line(colour="#d95f02", linewidth=.9) +
  geom_abline(linetype="dashed", colour="grey60") +
  annotate("text", x=.62, y=.12, label=sprintf("Nested-CV AUC = %.3f", auc_nested)) +
  labs(title="MSI-H vs MSS (random forest)", x="False positive rate", y="True positive rate") +
  theme_bw(base_size=12)
p_box <- ggplot(data.frame(oof, msi), aes(msi, oof, colour=msi)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(width=.15, alpha=.5, size=.9) +
  scale_colour_manual(values=c(MSS="#1b9e77", MSI_H="#d95f02")) +
  labs(x=NULL, y="Predicted P(MSI-H)",
       subtitle=sprintf("Wilcoxon p = %.1e", wilcox.test(oof ~ msi)$p.value)) +
  theme_bw(base_size=12) + theme(legend.position="none")
png("figures/fig7_msi.png", width=9, height=4, units="in", res=300)
gridExtra::grid.arrange(p_roc, p_box, ncol=2); dev.off()

sink("output/msi_metrics.txt")
cat("MSI-H vs MSS prediction (TCGA COAD) — Random Forest / Boruta\n")
cat("============================================================\n")
cat(sprintf("Samples: %d (MSI-H=%d, MSS=%d)\n", length(y), sum(y), sum(1-y)))
cat(sprintf("Honest nested-CV AUC: %.3f\n", auc_nested))
cat(sprintf("Wilcoxon p (P(MSI-H) ~ group): %.2e\n", wilcox.test(oof ~ msi)$p.value))
cat(sprintf("Signature (%d genes): %s\n", length(sig_genes), paste(sig_genes, collapse=", ")))
sink()
cat("DONE.\n")
