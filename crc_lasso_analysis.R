# =============================================================================
# LASSO-based gene signature for colorectal cancer (TCGA COAD) classification
# Corrected re-analysis: tumor vs normal defined by TCGA barcode sample-type
#   -01 = primary tumor, -11 = solid tissue normal
# Data: UCSC Xena, TCGA.COAD.sampleMap/HiSeqV2 (log2(norm_count+1)),
#       20,530 genes x 329 samples.  https://xenabrowser.net/
# Author: Xia Wu
# =============================================================================
set.seed(2026)
suppressMessages({
  library(limma); library(glmnet); library(ggplot2); library(pheatmap)
})

dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)

# ---- 1. Load expression matrix ---------------------------------------------
raw   <- read.delim("data/HiSeqV2", check.names = FALSE)
genes <- raw[[1]]
expr  <- as.matrix(raw[, -1]); rownames(expr) <- genes
cat(sprintf("Loaded %d genes x %d samples\n", nrow(expr), ncol(expr)))

# ---- 2. CORRECT tumor/normal labels from barcode suffix --------------------
suffix <- sub(".*-(\\d+)$", "\\1", colnames(expr))
keep   <- suffix %in% c("01", "11")            # primary tumor vs normal
expr   <- expr[, keep]
group  <- factor(ifelse(suffix[keep] == "01", "Tumor", "Normal"),
                 levels = c("Normal", "Tumor"))
cat("Sample groups (corrected):\n"); print(table(group))

# Drop genes with no variance
expr <- expr[apply(expr, 1, function(x) var(x, na.rm = TRUE) > 0), ]

# ---- 3. Differential expression (limma) ------------------------------------
design <- model.matrix(~0 + group); colnames(design) <- levels(group)
fit  <- lmFit(expr, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Tumor - Normal, levels = design)))
deg  <- topTable(fit2, coef = 1, number = Inf)
deg$gene <- rownames(deg)

logFC_t <- 2; padj_t <- 0.05
deg$change <- with(deg, ifelse(adj.P.Val < padj_t & logFC >  logFC_t, "Up",
                        ifelse(adj.P.Val < padj_t & logFC < -logFC_t, "Down", "Stable")))
sig <- subset(deg, change != "Stable")
cat(sprintf("DEGs (|log2FC|>%g & adj.P<%g): %d  (Up=%d, Down=%d)\n",
            logFC_t, padj_t, nrow(sig), sum(sig$change=="Up"), sum(sig$change=="Down")))
write.csv(deg, "output/limma_all_genes.csv", row.names = FALSE)
write.csv(sig, "output/DEGs.csv", row.names = FALSE)

# ---- 4. Volcano (Fig 1) -----------------------------------------------------
p_vol <- ggplot(deg, aes(logFC, -log10(adj.P.Val), colour = change)) +
  geom_point(alpha = .5, size = 1.6) +
  scale_colour_manual(values = c(Down="#2c7fb8", Stable="#bdbdbd", Up="#d7191c")) +
  geom_vline(xintercept = c(-logFC_t, logFC_t), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(padj_t), linetype = "dashed", colour = "grey50") +
  labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adj.P), colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave("figures/fig1_volcano.png", p_vol, width = 5, height = 4.2, dpi = 300)

# ---- 5. PCA: all genes vs DEGs (Fig 2a) ------------------------------------
pca_df <- function(mat) {
  pc <- prcomp(t(mat), scale. = TRUE)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 2)
  list(d = data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], group = group), ve = ve)
}
pa <- pca_df(expr); pd <- pca_df(expr[sig$gene, ])
cat(sprintf("PCA PC1 variance: all genes=%.2f%%, DEGs=%.2f%%\n", pa$ve[1], pd$ve[1]))
mk <- function(o, ttl) ggplot(o$d, aes(PC1, PC2, colour = group)) +
  geom_point(alpha = .7, size = 1.5) + stat_ellipse() +
  scale_colour_manual(values = c(Normal="#1b9e77", Tumor="#d95f02")) +
  labs(title = ttl, x = sprintf("PC1 (%.1f%%)", o$ve[1]),
       y = sprintf("PC2 (%.1f%%)", o$ve[2]), colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
png("figures/fig2a_pca.png", width = 9, height = 4, units = "in", res = 300)
gridExtra::grid.arrange(mk(pa, "All 20,530 genes"), mk(pd, sprintf("%d DEGs", nrow(sig))), ncol = 2)
dev.off()

# ---- 6. Heatmap of DEGs (Fig 2b) -------------------------------------------
ann <- data.frame(Group = group); rownames(ann) <- colnames(expr)
top_deg <- sig$gene[order(-abs(sig$logFC))]; top_deg <- head(top_deg, 50)
pheatmap(expr[top_deg, ], scale = "row", show_colnames = FALSE,
         annotation_col = ann, fontsize_row = 5,
         annotation_colors = list(Group = c(Normal="#1b9e77", Tumor="#d95f02")),
         color = colorRampPalette(c("#2c7fb8","white","#d7191c"))(100),
         filename = "figures/fig2b_heatmap.png", width = 6, height = 7)

# ---- 7. LASSO logistic regression (cross-validated) ------------------------
X <- t(expr[sig$gene, ]); y <- as.numeric(group == "Tumor")
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                   nfolds = 10, type.measure = "deviance", keep = TRUE)

idx_min <- cvfit$index["min", 1]; idx_1se <- cvfit$index["1se", 1]
sel <- function(s) { b <- coef(cvfit, s = s); rownames(b)[as.numeric(b) != 0][-1] }
g_min <- sel(cvfit$lambda.min); g_1se <- sel(cvfit$lambda.1se)
cat(sprintf("lambda.min=%.4f -> %d genes; lambda.1se=%.4f -> %d genes\n",
            cvfit$lambda.min, length(g_min), cvfit$lambda.1se, length(g_1se)))

# Honest, cross-validated (pre-validated) predictions for ROC/AUC
auc_fun <- function(score, lab){ r <- rank(score); np <- sum(lab==1); nn <- sum(lab==0)
  (sum(r[lab==1]) - np*(np+1)/2) / (np*nn) }
roc_pts <- function(score, lab){ o <- order(score, decreasing=TRUE)
  tp <- cumsum(lab[o]==1)/sum(lab==1); fp <- cumsum(lab[o]==0)/sum(lab==0)
  data.frame(fpr=c(0,fp), tpr=c(0,tp)) }
pre_min <- cvfit$fit.preval[, idx_min]; pre_1se <- cvfit$fit.preval[, idx_1se]
auc_min <- auc_fun(pre_min, y); auc_1se <- auc_fun(pre_1se, y)
cat(sprintf("Cross-validated AUC: lambda.min=%.3f, lambda.1se=%.3f\n", auc_min, auc_1se))

# Fig 3a: CV curve
png("figures/fig3a_lasso_cv.png", width = 5, height = 4.2, units = "in", res = 300)
plot(cvfit); dev.off()

# Fig 3b: ROC
r1 <- roc_pts(pre_min, y); r2 <- roc_pts(pre_1se, y)
p_roc <- ggplot() +
  geom_line(data=r1, aes(fpr,tpr,colour="lambda.min"), linewidth=.8) +
  geom_line(data=r2, aes(fpr,tpr,colour="lambda.1se"), linewidth=.8) +
  geom_abline(linetype="dashed", colour="grey60") +
  annotate("text", x=.6, y=.18, label=sprintf("AUC(min)=%.3f", auc_min), colour="#d95f02") +
  annotate("text", x=.6, y=.08, label=sprintf("AUC(1se)=%.3f", auc_1se), colour="#1f78b4") +
  scale_colour_manual(values=c("lambda.min"="#d95f02","lambda.1se"="#1f78b4")) +
  labs(x="False positive rate", y="True positive rate", colour=NULL) +
  theme_bw(base_size=12) + theme(legend.position=c(.7,.4))
ggsave("figures/fig3b_roc.png", p_roc, width=4.6, height=4.2, dpi=300)

# Predicted-probability boxplot + Wilcoxon
prob <- plogis(pre_min)
wt <- wilcox.test(prob ~ group)
p_box <- ggplot(data.frame(prob, group), aes(group, prob, colour=group)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(width=.15, alpha=.5, size=.8) +
  scale_colour_manual(values=c(Normal="#1b9e77", Tumor="#d95f02")) +
  labs(x=NULL, y="Predicted probability (tumor)",
       subtitle=sprintf("Wilcoxon p = %.2e", wt$p.value)) +
  theme_bw(base_size=12) + theme(legend.position="none")
ggsave("figures/fig_boxplot.png", p_box, width=4, height=4.2, dpi=300)

# ---- 8. Save metrics --------------------------------------------------------
sink("output/metrics.txt")
cat("Corrected TCGA COAD analysis — metrics\n=====================================\n")
cat(sprintf("Samples: Tumor(-01)=%d, Normal(-11)=%d\n", sum(group=="Tumor"), sum(group=="Normal")))
cat(sprintf("DEGs: %d (Up=%d, Down=%d)\n", nrow(sig), sum(sig$change=="Up"), sum(sig$change=="Down")))
cat(sprintf("PCA PC1: all=%.2f%%, DEG=%.2f%%\n", pa$ve[1], pd$ve[1]))
cat(sprintf("lambda.min=%.4f (%d genes): %s\n", cvfit$lambda.min, length(g_min), paste(g_min, collapse=", ")))
cat(sprintf("lambda.1se=%.4f (%d genes): %s\n", cvfit$lambda.1se, length(g_1se), paste(g_1se, collapse=", ")))
cat(sprintf("CV AUC: min=%.3f, 1se=%.3f\n", auc_min, auc_1se))
cat(sprintf("Wilcoxon p (prob ~ group): %.3e\n", wt$p.value))
sink()
cat("\nDONE. Figures in figures/, tables in output/\n")
