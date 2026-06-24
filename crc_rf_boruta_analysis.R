# =============================================================================
# Random-Forest / Boruta gene signature for colorectal cancer (TCGA COAD)
# Tumor vs normal defined by TCGA barcode sample-type (-01 tumor, -11 normal)
# Data: UCSC Xena, TCGA.COAD.sampleMap/HiSeqV2 (log2(norm_count+1)),
#       20,530 genes x 329 samples.  https://xenabrowser.net/
# Author: Xia Wu
# =============================================================================
set.seed(2026)
suppressMessages({
  library(limma); library(ggplot2); library(pheatmap)
  library(Boruta); library(ranger)
})
dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)

# ---- 1. Load + label --------------------------------------------------------
raw  <- read.delim("data/HiSeqV2", check.names = FALSE)
expr <- as.matrix(raw[, -1]); rownames(expr) <- raw[[1]]
suffix <- sub(".*-(\\d+)$", "\\1", colnames(expr))
keep   <- suffix %in% c("01", "11")
expr   <- expr[, keep]
group  <- factor(ifelse(suffix[keep] == "01", "Tumor", "Normal"), c("Normal", "Tumor"))
expr   <- expr[apply(expr, 1, function(x) var(x, na.rm = TRUE) > 0), ]
cat("Samples:\n"); print(table(group))

# ---- 2. Differential expression (limma) ------------------------------------
design <- model.matrix(~0 + group); colnames(design) <- levels(group)
fit2 <- eBayes(contrasts.fit(lmFit(expr, design),
               makeContrasts(Tumor - Normal, levels = design)))
deg  <- topTable(fit2, coef = 1, number = Inf); deg$gene <- rownames(deg)
logFC_t <- 2; padj_t <- 0.05
deg$change <- with(deg, ifelse(adj.P.Val < padj_t & logFC >  logFC_t, "Up",
                        ifelse(adj.P.Val < padj_t & logFC < -logFC_t, "Down", "Stable")))
sig <- subset(deg, change != "Stable")
cat(sprintf("DEGs: %d (Up=%d, Down=%d)\n", nrow(sig),
            sum(sig$change=="Up"), sum(sig$change=="Down")))
write.csv(deg, "output/limma_all_genes.csv", row.names = FALSE)
write.csv(sig, "output/DEGs.csv", row.names = FALSE)

# ---- 3. Volcano (Fig 1) -----------------------------------------------------
ggsave("figures/fig1_volcano.png",
  ggplot(deg, aes(logFC, -log10(adj.P.Val), colour = change)) +
    geom_point(alpha=.5, size=1.6) +
    scale_colour_manual(values=c(Down="#2c7fb8", Stable="#bdbdbd", Up="#d7191c")) +
    geom_vline(xintercept=c(-logFC_t,logFC_t), linetype="dashed", colour="grey50") +
    geom_hline(yintercept=-log10(padj_t), linetype="dashed", colour="grey50") +
    labs(x=expression(log[2]~fold~change), y=expression(-log[10]~adj.P), colour=NULL) +
    theme_bw(base_size=12) + theme(legend.position="top"),
  width=5, height=4.2, dpi=300)

# ---- 4. PCA all vs DEG (Fig 2a) --------------------------------------------
pca_df <- function(mat){ pc<-prcomp(t(mat), scale.=TRUE)
  ve<-round(100*pc$sdev^2/sum(pc$sdev^2),2)
  list(d=data.frame(PC1=pc$x[,1],PC2=pc$x[,2],group=group), ve=ve) }
pa<-pca_df(expr); pd<-pca_df(expr[sig$gene,])
cat(sprintf("PCA PC1: all=%.2f%%, DEG=%.2f%%\n", pa$ve[1], pd$ve[1]))
mk<-function(o,t) ggplot(o$d,aes(PC1,PC2,colour=group))+geom_point(alpha=.7,size=1.5)+
  stat_ellipse()+scale_colour_manual(values=c(Normal="#1b9e77",Tumor="#d95f02"))+
  labs(title=t,x=sprintf("PC1 (%.1f%%)",o$ve[1]),y=sprintf("PC2 (%.1f%%)",o$ve[2]),colour=NULL)+
  theme_bw(base_size=11)+theme(legend.position="top")
png("figures/fig2a_pca.png", width=9, height=4, units="in", res=300)
gridExtra::grid.arrange(mk(pa,"All 20,530 genes"), mk(pd,sprintf("%d DEGs",nrow(sig))), ncol=2)
dev.off()

# ---- 5. Heatmap (Fig 2b) ---------------------------------------------------
ann<-data.frame(Group=group); rownames(ann)<-colnames(expr)
top_deg<-head(sig$gene[order(-abs(sig$logFC))],50)
pheatmap(expr[top_deg,], scale="row", show_colnames=FALSE, annotation_col=ann,
  fontsize_row=5, annotation_colors=list(Group=c(Normal="#1b9e77",Tumor="#d95f02")),
  color=colorRampPalette(c("#2c7fb8","white","#d7191c"))(100),
  filename="figures/fig2b_heatmap.png", width=6, height=7)

# ---- 6. Boruta feature selection + Random Forest ---------------------------
Xall <- t(expr[sig$gene, ])                     # samples x DEGs
orig_names <- colnames(Xall)                    # original gene symbols
colnames(Xall) <- make.names(orig_names)        # syntactically valid for Boruta
name_map <- setNames(orig_names, colnames(Xall))# clean -> original symbol
set.seed(2026)
bor  <- Boruta(x = Xall, y = group, maxRuns = 100, doTrace = 0)
bor  <- TentativeRoughFix(bor)
conf <- getSelectedAttributes(bor, withTentative = FALSE)
impd <- attStats(bor); impd <- impd[impd$decision == "Confirmed", ]
impd <- impd[order(-impd$meanImp), ]
top  <- head(rownames(impd), 20)                # compact signature (clean names)
top_sym <- unname(name_map[top])                # original symbols for display
cat(sprintf("Boruta confirmed %d genes; top-20 signature: %s\n",
            length(conf), paste(top_sym, collapse=", ")))
write.csv(data.frame(gene=unname(name_map[rownames(impd)]), impd, row.names=NULL),
          "output/boruta_confirmed.csv", row.names=FALSE)

auc_fun <- function(s,l){ r<-rank(s); np<-sum(l==1); nn<-sum(l==0)
  (sum(r[l==1])-np*(np+1)/2)/(np*nn) }
roc_pts <- function(s,l){ o<-order(s,decreasing=TRUE)
  data.frame(fpr=c(0,cumsum(l[o]==0)/sum(l==0)), tpr=c(0,cumsum(l[o]==1)/sum(l==1))) }

# 10-fold stratified CV with random forest (ranger), out-of-fold probabilities
y <- as.numeric(group == "Tumor"); nfold <- 10
set.seed(2026); fold <- integer(length(y))
for (cl in c(0,1)){ ix<-which(y==cl); fold[ix]<-sample(rep(1:nfold, length.out=length(ix))) }
cv_prob <- function(feat){
  p <- numeric(length(y))
  for (k in 1:nfold){ tr<-fold!=k; te<-fold==k
    rf <- ranger(x=Xall[tr,feat,drop=FALSE], y=group[tr],
                 probability=TRUE, num.trees=500, seed=2026)
    p[te] <- predict(rf, Xall[te,feat,drop=FALSE])$predictions[,"Tumor"] }
  p }
prob_all <- cv_prob(conf); prob_top <- cv_prob(top)
auc_all <- auc_fun(prob_all, y); auc_top <- auc_fun(prob_top, y)
cat(sprintf("RF 10-fold CV AUC: all-confirmed(%d)=%.3f, top-20=%.3f\n",
            length(conf), auc_all, auc_top))

# Fig 4: Boruta top-20 importance, coloured by DE direction
dir <- deg$change[match(top_sym, deg$gene)]
bdf <- data.frame(gene=factor(top_sym, levels=rev(top_sym)),
                  imp=impd[top,"meanImp"], dir=dir)
ggsave("figures/fig4_boruta_importance.png",
  ggplot(bdf, aes(imp, gene, fill=dir)) + geom_col() +
    scale_fill_manual(values=c(Up="#d7191c", Down="#2c7fb8"), name=NULL) +
    labs(x="Boruta importance (mean Z-score)", y=NULL) +
    theme_bw(base_size=11) + theme(legend.position="top"),
  width=4.6, height=5, dpi=300)

# Fig 5: RF cross-validated ROC
r1<-roc_pts(prob_all,y); r2<-roc_pts(prob_top,y)
ggsave("figures/fig5_rf_roc.png",
  ggplot()+geom_line(data=r1,aes(fpr,tpr,colour="All confirmed"),linewidth=.8)+
    geom_line(data=r2,aes(fpr,tpr,colour="Top-20 signature"),linewidth=.8)+
    geom_abline(linetype="dashed",colour="grey60")+
    annotate("text",x=.62,y=.18,label=sprintf("AUC(all)=%.3f",auc_all),colour="#1b9e77")+
    annotate("text",x=.62,y=.08,label=sprintf("AUC(top-20)=%.3f",auc_top),colour="#d95f02")+
    scale_colour_manual(values=c("All confirmed"="#1b9e77","Top-20 signature"="#d95f02"))+
    labs(x="False positive rate",y="True positive rate",colour=NULL)+
    theme_bw(base_size=12)+theme(legend.position=c(.68,.4)),
  width=4.6, height=4.2, dpi=300)

# ---- 7. Functional enrichment (GO BP + KEGG) -------------------------------
suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
to_entrez<-function(s) suppressWarnings(bitr(s,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID)
gl<-list(`Up-regulated`=to_entrez(sig$gene[sig$change=="Up"]),
         `Down-regulated`=to_entrez(sig$gene[sig$change=="Down"]))
cc<-compareCluster(gl, fun="enrichGO", OrgDb=org.Hs.eg.db, ont="BP",
                   pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
cc<-clusterProfiler::simplify(cc, cutoff=0.6)
ggsave("figures/fig6_enrichment.png",
  dotplot(cc, showCategory=6, font.size=9)+
    ggplot2::theme(axis.text.y=ggplot2::element_text(size=8)),
  width=7.2, height=6, dpi=300)
write.csv(as.data.frame(cc),"output/GO_enrichment.csv", row.names=FALSE)
tryCatch({ k<-compareCluster(gl,fun="enrichKEGG",organism="hsa",pvalueCutoff=0.05)
  write.csv(as.data.frame(k),"output/KEGG_enrichment.csv", row.names=FALSE) },
  error=function(e) cat("KEGG skipped:", conditionMessage(e), "\n"))

# ---- 8. Metrics -------------------------------------------------------------
sink("output/metrics.txt")
cat("RF/Boruta TCGA COAD analysis — metrics\n======================================\n")
cat(sprintf("Samples: Tumor=%d, Normal=%d\n", sum(group=="Tumor"), sum(group=="Normal")))
cat(sprintf("DEGs: %d (Up=%d, Down=%d)\n", nrow(sig), sum(sig$change=="Up"), sum(sig$change=="Down")))
cat(sprintf("PCA PC1: all=%.2f%%, DEG=%.2f%%\n", pa$ve[1], pd$ve[1]))
cat(sprintf("Boruta confirmed: %d genes\n", length(conf)))
cat(sprintf("Top-20 signature: %s\n", paste(top_sym, collapse=", ")))
cat(sprintf("RF 10-fold CV AUC: all-confirmed=%.3f, top-20=%.3f\n", auc_all, auc_top))
sink()
cat("DONE\n")
