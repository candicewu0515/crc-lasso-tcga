# =============================================================================
# Harder endpoint: predicting microsatellite-instability (MSI) status in TCGA
# COAD from expression. MSI labels: cBioPortal PanCancer Atlas MANTIS score
# (MSI-H = MANTIS >= 0.4), cross-checked against MSIsensor (95% concordant).
# Honest performance via NESTED 10-fold CV (DEG selection inside each fold).
# Author: Xia Wu
# =============================================================================
set.seed(2026)
suppressMessages({ library(limma); library(glmnet); library(ggplot2) })
dir.create("figures", showWarnings = FALSE); dir.create("output", showWarnings = FALSE)

# ---- Load expression + MSI labels ------------------------------------------
raw  <- read.delim("data/HiSeqV2", check.names = FALSE)
expr <- as.matrix(raw[, -1]); rownames(expr) <- raw[[1]]
expr <- expr[apply(expr, 1, function(x) var(x) > 0), ]

mantis <- read.delim("data/MSI_SCORE_MANTIS.tsv", check.names = FALSE)
lab <- setNames(mantis$MSI_SCORE_MANTIS, mantis$sampleId)
common <- intersect(colnames(expr), names(lab))
expr <- expr[, common]
msi  <- factor(ifelse(lab[common] >= 0.4, "MSI_H", "MSS"), levels = c("MSS", "MSI_H"))
y <- as.numeric(msi == "MSI_H")
cat(sprintf("Samples: %d (MSI-H=%d, MSS=%d)\n", length(y), sum(y), sum(1 - y)))

auc_fun <- function(s, l){ r <- rank(s); np <- sum(l==1); nn <- sum(l==0)
  (sum(r[l==1]) - np*(np+1)/2)/(np*nn) }
deg_select <- function(mat, grp){
  d <- model.matrix(~0 + grp); colnames(d) <- levels(grp)
  f <- eBayes(contrasts.fit(lmFit(mat, d), makeContrasts(MSI_H - MSS, levels = d)))
  tt <- topTable(f, coef = 1, number = Inf)
  g <- rownames(tt)[abs(tt$logFC) > 1 & tt$adj.P.Val < 0.05]
  if (length(g) < 10) g <- rownames(tt)[order(tt$P.Value)][1:50]
  g
}

# ---- NESTED 10-fold CV: DEG + LASSO trained only on the inner training set --
K <- 10
folds <- integer(length(y))
for (cl in 0:1){ ix <- which(y == cl); folds[ix] <- sample(rep(1:K, length.out = length(ix))) }
oof <- rep(NA, length(y))               # out-of-fold predicted probabilities
nsel <- c()
for (k in 1:K){
  tr <- folds != k; te <- folds == k
  deg <- deg_select(expr[, tr], droplevels(msi[tr]))
  nsel <- c(nsel, length(deg))
  Xtr <- t(expr[deg, tr]); Xte <- t(expr[deg, te])
  cv  <- cv.glmnet(Xtr, y[tr], family = "binomial", alpha = 1, nfolds = 10)
  oof[te] <- as.numeric(predict(cv, newx = Xte, s = cv$lambda.min, type = "response"))
}
auc_nested <- auc_fun(oof, y)
cat(sprintf("Nested-CV genes/fold: median=%d  |  HONEST nested-CV AUC = %.3f\n",
            as.integer(median(nsel)), auc_nested))

# ---- Final signature on full data (for reporting) --------------------------
deg_all <- deg_select(expr, msi)
Xall <- t(expr[deg_all, ])
cvf  <- cv.glmnet(Xall, y, family = "binomial", alpha = 1, nfolds = 10)
sig_genes <- { b <- coef(cvf, s = cvf$lambda.1se); rownames(b)[as.numeric(b) != 0][-1] }
cat(sprintf("Final MSI signature (lambda.1se): %d genes: %s\n",
            length(sig_genes), paste(sig_genes, collapse = ", ")))

# ---- Figure: ROC + predicted-probability boxplot ---------------------------
roc_pts <- function(s,l){ o <- order(s, decreasing=TRUE)
  data.frame(fpr=c(0,cumsum(l[o]==0)/sum(l==0)), tpr=c(0,cumsum(l[o]==1)/sum(l==1))) }
rp <- roc_pts(oof, y)
p_roc <- ggplot(rp, aes(fpr, tpr)) + geom_line(colour="#d95f02", linewidth=.9) +
  geom_abline(linetype="dashed", colour="grey60") +
  annotate("text", x=.62, y=.12, label=sprintf("Nested-CV AUC = %.3f", auc_nested)) +
  labs(title="MSI-H vs MSS", x="False positive rate", y="True positive rate") +
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
cat("MSI-H vs MSS prediction (TCGA COAD)\n===================================\n")
cat(sprintf("Samples: %d (MSI-H=%d, MSS=%d)\n", length(y), sum(y), sum(1-y)))
cat(sprintf("Honest nested-CV AUC: %.3f\n", auc_nested))
cat(sprintf("Final signature (%d genes): %s\n", length(sig_genes), paste(sig_genes, collapse=", ")))
sink()
cat("DONE.\n")
