# =============================================================================
# Revision: Independent MGS reconstruction in each external cohort
# No cross-platform projection - each cohort analyzed independently
# =============================================================================

library(dplyr)
library(readxl)
library(princurve)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"

# ---- Load discovery TPM ----
tpm_raw <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/TPM.xlsx")
genes <- tpm_raw$gene_name
keep <- !duplicated(genes)
tpm_raw <- tpm_raw[keep, ]; genes <- genes[keep]

our_samples <- c('tpm.P1','tpm.P10','tpm.P11','tpm.P12','tpm.P13','tpm.P14','tpm.P15',
                 'tpm.P16','tpm.P17','tpm.P18','tpm.P19','tpm.P2','tpm.P20','tpm.P21',
                 'tpm.P22','tpm.P23','tpm.P24','tpm.P25','tpm.P26','tpm.P27','tpm.P28',
                 'tpm.P3','tpm.P30','tpm.P32','tpm.P34','tpm.P35','tpm.P36','tpm.P37',
                 'tpm.P38','tpm.P39','tpm.P4','tpm.P40','tpm.P41','tpm.P5','tpm.P6',
                 'tpm.P7','tpm.S5','tpm.S8')

tpm_disc <- as.matrix(tpm_raw[, our_samples])
rownames(tpm_disc) <- genes
log_tpm_disc <- log2(tpm_disc + 1)
cat(sprintf("Discovery: %d genes x %d samples\n", nrow(log_tpm_disc), ncol(log_tpm_disc)))

# ---- Load external cohorts (with gene symbols) ----
load_mat <- function(path) {
  dat <- read.csv(path, check.names = FALSE)
  g <- as.character(dat[, 1])
  kp <- g != "" & !is.na(g)
  dat <- dat[kp, ]; g <- g[kp]
  dup <- duplicated(g); dat <- dat[!dup, ]; g <- g[!dup]
  mat <- as.matrix(dat[, -1, drop = FALSE])
  rownames(mat) <- g
  return(mat)
}

gse141801 <- load_mat("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE141801_Gugel/GSE141801_Expression_Log2.csv")
gse39645 <- load_mat("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE39645_Torres/GSE39645_Expression_Log2.csv")
cat(sprintf("GSE141801: %d x %d\n", nrow(gse141801), ncol(gse141801)))
cat(sprintf("GSE39645: %d x %d\n", nrow(gse39645), ncol(gse39645)))

# ---- Function: Build MGS independently in a cohort ----
build_mgs_independent <- function(expr_mat, n_genes = 3000, tag = "") {
  # expr_mat: genes x samples, already in log2 space or linear
  mads <- apply(expr_mat, 1, mad, na.rm = TRUE)
  n_use <- min(n_genes, sum(mads > 0, na.rm = TRUE))
  top_genes <- names(sort(mads, decreasing = TRUE))[1:n_use]
  cat(sprintf("  [%s] Using %d genes\n", tag, n_use))

  pca <- prcomp(t(expr_mat[top_genes, ]), center = TRUE, scale. = TRUE)
  var5 <- 100 * sum(summary(pca)$importance[2, 1:5])
  cat(sprintf("  [%s] PC1=%.1f%%, PC2=%.1f%%, PC1-5=%.1f%%\n",
              tag, 100*summary(pca)$importance[2,1], 100*summary(pca)$importance[2,2], var5))

  fit <- principal_curve(as.matrix(pca$x[, 1:min(5, ncol(pca$x))]), stretch = 0)
  mgs <- fit$lambda
  mgs <- (mgs - min(mgs)) / (max(mgs) - min(mgs))

  return(list(
    mgs = setNames(mgs, colnames(expr_mat)),
    pca = pca,
    fit = fit,
    var5 = var5,
    top_genes = top_genes,
    pc1_var = 100 * summary(pca)$importance[2, 1]
  ))
}

# ---- 1. Discovery MGS ----
cat("\n========== 1. Discovery MGS ==========\n")
mgs_disc <- build_mgs_independent(log_tpm_disc, 3000, "DISCOVERY")

# Subtype annotation
cl <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]
sample_ids_clean <- gsub("^tpm\\.", "", our_samples)
subtypes_disc <- setNames(cl$Subtype, cl$SampleID)[sample_ids_clean]

mgs_disc_df <- data.frame(Sample = names(mgs_disc$mgs), MGS = mgs_disc$mgs, Subtype = subtypes_disc)
cat("Discovery MGS by subtype:\n")
print(mgs_disc_df %>% group_by(Subtype) %>% summarise(mean=mean(MGS), sd=sd(MGS)))
kw_disc <- kruskal.test(MGS ~ Subtype, data = mgs_disc_df)
cat(sprintf("Kruskal-Wallis P = %.2e\n", kw_disc$p.value))

# ---- 2. GSE141801 independent MGS ----
cat("\n========== 2. GSE141801 Independent MGS ==========\n")
# Data is log2-normalized expression from microarray
mgs_141801_ind <- build_mgs_independent(gse141801, 3000, "GSE141801")

# Load clinical data
clin_141801 <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE141801_Gugel/Clinical_Data.csv")
clin_141801$MGS <- mgs_141801_ind$mgs[clin_141801$Sample_geo_accession]

cat("GSE141801 MGS by Barrett Class:\n")
print(clin_141801 %>% filter(!is.na(MGS)) %>% group_by(Class) %>%
  summarise(mean_MGS=mean(MGS), sd_MGS=sd(MGS), n=n()))
kw_class <- kruskal.test(MGS ~ Class, data = clin_141801[!is.na(clin_141801$MGS), ])
cat(sprintf("Kruskal-Wallis (Class) P = %.4f\n", kw_class$p.value))

cat("GSE141801 MGS by NF2:\n")
nf2_data <- clin_141801[clin_141801$NF2 %in% c("Y","N") & !is.na(clin_141801$MGS), ]
print(nf2_data %>% group_by(NF2) %>% summarise(mean_MGS=mean(MGS), sd=sd(MGS), n=n()))
wt_mgs <- wilcox.test(MGS ~ NF2, data = nf2_data)
cat(sprintf("Wilcoxon P = %.4f\n", wt_mgs$p.value))

cat("GSE141801 MGS by Size:\n")
print(clin_141801 %>% filter(!is.na(MGS) & !is.na(Size)) %>% group_by(Size) %>%
  summarise(mean_MGS=mean(MGS), sd=sd(MGS), n=n()))

# ---- 3. GSE39645 independent MGS ----
cat("\n========== 3. GSE39645 Independent MGS ==========\n")
mgs_39645_ind <- build_mgs_independent(gse39645, 3000, "GSE39645")

clin_39645 <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE39645_Torres/Clinical_Data.csv")
clin_39645$MGS <- mgs_39645_ind$mgs[clin_39645$Sample_geo_accession]

cat("GSE39645 MGS by Barrett Class:\n")
print(clin_39645 %>% filter(!is.na(MGS)) %>% group_by(Class) %>%
  summarise(mean_MGS=mean(MGS), sd=sd(MGS), n=n()))
kw_class_39 <- kruskal.test(MGS ~ Class, data = clin_39645[!is.na(clin_39645$MGS), ])
cat(sprintf("Kruskal-Wallis (Class) P = %.4f\n", kw_class_39$p.value))

cat("GSE39645 MGS by NF2:\n")
nf2_39 <- clin_39645[clin_39645$NF2 %in% c("Y","N") & !is.na(clin_39645$MGS), ]
if(nrow(nf2_39) > 5) {
  print(nf2_39 %>% group_by(NF2) %>% summarise(mean_MGS=mean(MGS), sd=sd(MGS), n=n()))
  cat(sprintf("Wilcoxon P = %.4f\n", wilcox.test(MGS ~ NF2, data = nf2_39)$p.value))
}

# ---- 4. Compare MGS distributions across cohorts ----
cat("\n========== 4. Cross-cohort MGS Comparison ==========\n")
cat(sprintf("%-15s %8s %8s %8s\n", "Cohort", "Mean", "SD", "IQR"))
cat(sprintf("%-15s %8.3f %8.3f [%.3f,%.3f]\n", "Discovery", mean(mgs_disc$mgs), sd(mgs_disc$mgs),
            quantile(mgs_disc$mgs, 0.25), quantile(mgs_disc$mgs, 0.75)))
cat(sprintf("%-15s %8.3f %8.3f [%.3f,%.3f]\n", "GSE141801", mean(mgs_141801_ind$mgs, na.rm=TRUE),
            sd(mgs_141801_ind$mgs, na.rm=TRUE),
            quantile(mgs_141801_ind$mgs, 0.25, na.rm=TRUE), quantile(mgs_141801_ind$mgs, 0.75, na.rm=TRUE)))
cat(sprintf("%-15s %8.3f %8.3f [%.3f,%.3f]\n", "GSE39645", mean(mgs_39645_ind$mgs, na.rm=TRUE),
            sd(mgs_39645_ind$mgs, na.rm=TRUE),
            quantile(mgs_39645_ind$mgs, 0.25, na.rm=TRUE), quantile(mgs_39645_ind$mgs, 0.75, na.rm=TRUE)))

# ---- 5. Key question: Do external MGS gradients capture immune/stromal axes? ----
cat("\n========== 5. External MGS Biological Validation ==========\n")

# In GSE141801, compute immune/ECM signature scores and correlate with MGS
# Use the EstIMATE-like gene signatures
immune_genes <- c("CD2","CD3D","CD3E","CD4","CD8A","CD8B","GZMA","PRF1","CTLA4",
                   "HAVCR2","LAG3","PDCD1","CD68","CD163","CD14","CD74","HLA-DRA",
                   "HLA-DRB1","CCL5","CXCL9","CXCL10","CXCL13","CD79A","CD79B")
ecm_genes <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2",
               "COL6A3","FN1","LOX","LOXL1","LOXL2","SPARC","THBS2","TGFB1","CTGF")

im_avail_141801 <- intersect(immune_genes, rownames(gse141801))
ecm_avail_141801 <- intersect(ecm_genes, rownames(gse141801))

imm_score_141801 <- colMeans(gse141801[im_avail_141801, , drop = FALSE])
ecm_score_141801 <- colMeans(gse141801[ecm_avail_141801, , drop = FALSE])

clin_141801$Immune_Score <- imm_score_141801[clin_141801$Sample_geo_accession]
clin_141801$ECM_Score <- ecm_score_141801[clin_141801$Sample_geo_accession]

cat("GSE141801: Correlations with MGS (independent):\n")
mgs_cor_141801 <- clin_141801[!is.na(clin_141801$MGS), ]
cat(sprintf("  MGS vs Immune Score: rho=%.3f, P=%.4f\n",
            cor(mgs_cor_141801$MGS, mgs_cor_141801$Immune_Score, method="spearman"),
            cor.test(mgs_cor_141801$MGS, mgs_cor_141801$Immune_Score, method="spearman")$p.value))
cat(sprintf("  MGS vs ECM Score:    rho=%.3f, P=%.4f\n",
            cor(mgs_cor_141801$MGS, mgs_cor_141801$ECM_Score, method="spearman"),
            cor.test(mgs_cor_141801$MGS, mgs_cor_141801$ECM_Score, method="spearman")$p.value))

# Same for GSE39645
im_avail_39645 <- intersect(immune_genes, rownames(gse39645))
ecm_avail_39645 <- intersect(ecm_genes, rownames(gse39645))
imm_score_39645 <- colMeans(gse39645[im_avail_39645, , drop = FALSE])
ecm_score_39645 <- colMeans(gse39645[ecm_avail_39645, , drop = FALSE])
clin_39645$Immune_Score <- imm_score_39645[clin_39645$Sample_geo_accession]
clin_39645$ECM_Score <- ecm_score_39645[clin_39645$Sample_geo_accession]

cat("\nGSE39645: Correlations with MGS (independent):\n")
mgs_cor_39645 <- clin_39645[!is.na(clin_39645$MGS), ]
cat(sprintf("  MGS vs Immune Score: rho=%.3f, P=%.4f\n",
            cor(mgs_cor_39645$MGS, mgs_cor_39645$Immune_Score, method="spearman"),
            cor.test(mgs_cor_39645$MGS, mgs_cor_39645$Immune_Score, method="spearman")$p.value))
cat(sprintf("  MGS vs ECM Score:    rho=%.3f, P=%.4f\n",
            cor(mgs_cor_39645$MGS, mgs_cor_39645$ECM_Score, method="spearman"),
            cor.test(mgs_cor_39645$MGS, mgs_cor_39645$ECM_Score, method="spearman")$p.value))

# ---- 6. PC1 direction validation ----
cat("\n========== 6. PC1 Direction Validation ==========\n")
# In discovery: PC1 correlates with immune score, anti-correlates with purity
# Does the same hold in external cohorts?
disc_pc1 <- mgs_disc$pca$x[, 1]
cat(sprintf("Discovery PC1 vs MGS: rho=%.3f\n", cor(disc_pc1, mgs_disc$mgs, method="spearman")))

# MGS is arc-length along principal curve, which should follow PC1
# Check: in external cohorts, does MGS capture immune activation?
cat("\n--- PC1 loading analysis (top genes in discovery) ---\n")
pc1_loadings <- mgs_disc$pca$rotation[, 1]
top_pc1_pos <- names(sort(pc1_loadings, decreasing = TRUE))[1:20]
top_pc1_neg <- names(sort(pc1_loadings, decreasing = FALSE))[1:20]
cat("Top PC1 positive (PC1-high direction):\n")
cat(paste(top_pc1_pos, collapse=", "), "\n")
cat("Top PC1 negative (PC1-low direction):\n")
cat(paste(top_pc1_neg, collapse=", "), "\n")

# ---- 7. Summary table for paper ----
cat("\n========== 7. Summary: Independent MGS Across Cohorts ==========\n")
cat("Key results:\n")
cat(sprintf("  Discovery: MGS C1(%.2f±%.2f) < C2(%.2f±%.2f) < C3(%.2f±%.2f)\n",
            mean(mgs_disc_df$MGS[mgs_disc_df$Subtype=="C1"]),
            sd(mgs_disc_df$MGS[mgs_disc_df$Subtype=="C1"]),
            mean(mgs_disc_df$MGS[mgs_disc_df$Subtype=="C2"]),
            sd(mgs_disc_df$MGS[mgs_disc_df$Subtype=="C2"]),
            mean(mgs_disc_df$MGS[mgs_disc_df$Subtype=="C3"]),
            sd(mgs_disc_df$MGS[mgs_disc_df$Subtype=="C3"])))

mgs_141801_nona <- clin_141801[!is.na(clin_141801$MGS), ]
cat(sprintf("  GSE141801: Injury-like(%.2f±%.2f) vs nmSC Core(%.2f±%.2f)\n",
            mean(mgs_141801_nona$MGS[mgs_141801_nona$Class=="Injury-like"]),
            sd(mgs_141801_nona$MGS[mgs_141801_nona$Class=="Injury-like"]),
            mean(mgs_141801_nona$MGS[mgs_141801_nona$Class=="nmSC Core"]),
            sd(mgs_141801_nona$MGS[mgs_141801_nona$Class=="nmSC Core"])))

mgs_39645_nona <- clin_39645[!is.na(clin_39645$MGS), ]
cat(sprintf("  GSE39645: Injury-like(%.2f±%.2f) vs nmSC Core(%.2f±%.2f)\n",
            mean(mgs_39645_nona$MGS[mgs_39645_nona$Class=="Injury-like"]),
            sd(mgs_39645_nona$MGS[mgs_39645_nona$Class=="Injury-like"]),
            mean(mgs_39645_nona$MGS[mgs_39645_nona$Class=="nmSC Core"]),
            sd(mgs_39645_nona$MGS[mgs_39645_nona$Class=="nmSC Core"])))

# ---- Save ----
saveRDS(list(discovery = mgs_disc, gse141801 = mgs_141801_ind, gse39645 = mgs_39645_ind,
             clin_141801 = clin_141801, clin_39645 = clin_39645),
        file.path(work_dir, "mgs_independent_results.rds"))

cat("\n========== DONE ==========\n")
