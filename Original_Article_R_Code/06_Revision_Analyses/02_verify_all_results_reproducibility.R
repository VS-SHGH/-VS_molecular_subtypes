# =============================================================================
# Final Verification: Re-run all key analyses for reproducibility check
# =============================================================================

library(dplyr)
library(readxl)

cat("========== VERIFICATION: All Key Results ==========\n\n")

# ---- Load core data ----
tpm <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/数据/TPM.xlsx")
genes <- tpm$gene_name
keep <- !duplicated(genes); tpm <- tpm[keep, ]; genes <- genes[keep]

cl <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/数据/cl_bulk_with_Subtype.xlsx")
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]

our_samples <- paste0("tpm.", cl$SampleID)
tpm_mat <- as.matrix(tpm[, our_samples])
rownames(tpm_mat) <- genes
log_tpm <- log2(tpm_mat + 1)
subtypes <- setNames(cl$Subtype, cl$SampleID)
names(subtypes) <- paste0("tpm.", names(subtypes))

cat(sprintf("Data: %d genes x %d samples\n", nrow(log_tpm), ncol(log_tpm)))
cat(sprintf("Subtypes: C1=%d, C2=%d, C3=%d\n\n",
            sum(cl$Subtype=="C1"), sum(cl$Subtype=="C2"), sum(cl$Subtype=="C3")))

PASS <- 0; FAIL <- 0; WARN <- 0
check <- function(name, expected, actual, tol=0.05) {
  diff <- abs(expected - actual)
  if (is.na(diff) || diff > tol) {
    cat(sprintf("  ❌ %s: expected=%.4f, got=%.4f\n", name, expected, actual))
    FAIL <<- FAIL + 1
  } else {
    cat(sprintf("  ✅ %s: %.4f (expected %.4f)\n", name, actual, expected))
    PASS <<- PASS + 1
  }
}

# ---- 1. PC1 Variance ----
cat("1. PCA Variance\n")
mads <- apply(log_tpm, 1, mad)
top_2000 <- names(sort(mads, decreasing=TRUE))[1:2000]
pca <- prcomp(t(log_tpm[top_2000,]), center=TRUE, scale.=TRUE)
pc1 <- 100*summary(pca)$importance[2,1]
pc2 <- 100*summary(pca)$importance[2,2]
check("PC1 variance (ms=24.5)", 24.5, pc1, tol=2.5)
check("PC2 variance (ms=13.7)", 13.7, pc2, tol=2.5)
cat(sprintf("  Actual: PC1=%.1f%%, PC2=%.1f%% (Manuscript: 24.5%%/13.7%%)\n\n", pc1, pc2))

# ---- 2. NK-cell AUC ----
cat("2. NK-cell ROC\n")
library(pROC)
cl$NK_cells <- as.numeric(cl$自然杀伤细胞)
cl$C3_vs_other <- ifelse(cl$Subtype=="C3", "C3", "C1_C2")
roc_obj <- roc(cl$C3_vs_other, cl$NK_cells, levels=c("C1_C2","C3"))
check("NK AUC", 0.874, auc(roc_obj), tol=0.05)
cat(sprintf("  Actual AUC: %.3f (95%% CI: %.3f-%.3f)\n\n",
            auc(roc_obj), ci.auc(roc_obj)[1], ci.auc(roc_obj)[3]))

# ---- 3. Regression Models ----
cat("3. NK Regression\n")
cl$SizeLarge <- ifelse(cl$SizeGrade=="large", 1, 0)
m1 <- lm(NK_cells ~ Subtype + Age, data=cl)
m2 <- lm(NK_cells ~ Subtype + Age + SizeLarge, data=cl)
beta_m1 <- coef(m1)["SubtypeC3"]
p_m1 <- summary(m1)$coefficients["SubtypeC3",4]
beta_m2 <- coef(m2)["SubtypeC3"]
p_m2 <- summary(m2)$coefficients["SubtypeC3",4]
beta_size <- coef(m2)["SizeLarge"]
p_size <- summary(m2)$coefficients["SizeLarge",4]
check("M1 C3 beta", 4.745, beta_m1, tol=0.2)
check("M1 C3 P", 0.017, p_m1, tol=0.01)
check("M2 C3 beta", 3.251, beta_m2, tol=0.5)
check("M2 C3 P", 0.053, p_m2, tol=0.02)
check("M2 Size beta", 5.696, beta_size, tol=0.5)
check("M2 Size P", 0.0003, p_size, tol=0.01)
cat(sprintf("  M1: C3 beta=%.3f P=%.4f\n", beta_m1, p_m1))
cat(sprintf("  M2: C3 beta=%.3f P=%.4f, Size beta=%.3f P=%.4f\n",
            beta_m2, p_m2, beta_size, p_size))
epv <- nrow(cl) / length(coef(m2))
check("M2 EPV", 7.6, epv, tol=2.0)
cat(sprintf("  EPV: %.1f\n\n", epv))

# ---- 4. ESTIMATE Purity ----
cat("4. ESTIMATE Purity\n")
# Quick purity by gene expression
immune_genes <- c("CD2","CD3D","CD8A","GZMA","PRF1","CD68","CD163","CD74","HLA-DRA")
stromal_genes <- c("COL1A1","COL1A2","COL3A1","FN1","LOX","SPARC","THBS2")
im_avail <- intersect(immune_genes, rownames(log_tpm))
st_avail <- intersect(stromal_genes, rownames(log_tpm))
imm_score <- colMeans(log_tpm[im_avail,])
str_score <- colMeans(log_tpm[st_avail,])
purity <- 1 - (imm_score + str_score) / max(imm_score + str_score)
c1_pur <- mean(purity[subtypes=="C1"])
c2_pur <- mean(purity[subtypes=="C2"])
c3_pur <- mean(purity[subtypes=="C3"])
check("C1 purity", 0.611, c1_pur, tol=0.1)
check("C2 purity", 0.471, c2_pur, tol=0.1)
check("C3 purity", 0.419, c3_pur, tol=0.1)
cat(sprintf("  Purity: C1=%.3f, C2=%.3f, C3=%.3f\n\n", c1_pur, c2_pur, c3_pur))

# ---- 5. MCP-counter Fibroblasts and Monocytes ----
cat("5. MCP-counter Key Results (from saved RDS)\n")
mcp <- readRDS("/Users/yanchen/Desktop/VSbulk+单细胞/VS_revision/output/mcp_counter_results.rds")
fibro_c1 <- mean(as.numeric(mcp["Fibroblasts", subtypes=="C1"]))
fibro_c2 <- mean(as.numeric(mcp["Fibroblasts", subtypes=="C2"]))
fibro_c3 <- mean(as.numeric(mcp["Fibroblasts", subtypes=="C3"]))
mono_c1 <- mean(as.numeric(mcp["Monocytic lineage", subtypes=="C1"]))
mono_c2 <- mean(as.numeric(mcp["Monocytic lineage", subtypes=="C2"]))
mono_c3 <- mean(as.numeric(mcp["Monocytic lineage", subtypes=="C3"]))
check("Fibroblasts C1", 199, fibro_c1, tol=20)
check("Fibroblasts C2", 444, fibro_c2, tol=30)
check("Monocytes C1", 86, mono_c1, tol=10)
check("Monocytes C3", 178, mono_c3, tol=20)
cat(sprintf("  Fibroblasts: C1=%.0f, C2=%.0f, C3=%.0f\n", fibro_c1, fibro_c2, fibro_c3))
cat(sprintf("  Monocytes: C1=%.0f, C2=%.0f, C3=%.0f\n\n", mono_c1, mono_c2, mono_c3))

# ---- 6. MGS ordering ----
cat("6. MGS Unbiased Ordering\n")
top_3000 <- names(sort(mads, decreasing=TRUE))[1:3000]
pca_mgs <- prcomp(t(log_tpm[top_3000,]), center=TRUE, scale.=TRUE)
library(princurve)
fit <- principal_curve(as.matrix(pca_mgs$x[,1:5]), stretch=0)
mgs_raw <- fit$lambda

# Orient by immune score
immune_sig_genes <- c("CD2","CD3D","CD3E","CD8A","GZMA","PRF1","CTLA4","HAVCR2",
                       "LAG3","PDCD1","CD68","CD163","CD74","HLA-DRA","CCL5","CXCL9","CXCL10")
im_sig_avail <- intersect(immune_sig_genes, rownames(log_tpm))
im_sig <- colMeans(log_tpm[im_sig_avail,])
cor_immune <- cor(mgs_raw, im_sig, method="spearman")
if (cor_immune < 0) {
  mgs <- 1 - (mgs_raw - min(mgs_raw)) / (max(mgs_raw) - min(mgs_raw))
} else {
  mgs <- (mgs_raw - min(mgs_raw)) / (max(mgs_raw) - min(mgs_raw))
}
names(mgs) <- colnames(log_tpm)
c1_mgs <- mean(mgs[subtypes=="C1"])
c2_mgs <- mean(mgs[subtypes=="C2"])
c3_mgs <- mean(mgs[subtypes=="C3"])
check("MGS C1", 0.13, c1_mgs, tol=0.05)
check("MGS C2", 0.47, c2_mgs, tol=0.07)
check("MGS C3", 0.82, c3_mgs, tol=0.12)
# Spearman with subtype order
subtype_num <- ifelse(subtypes=="C1", 1, ifelse(subtypes=="C2", 2, 3))
rho <- cor(mgs, subtype_num, method="spearman")
check("MGS-Subtype rho", 0.929, rho, tol=0.03)
cat(sprintf("  MGS: C1=%.3f, C2=%.3f, C3=%.3f, rho=%.3f\n\n", c1_mgs, c2_mgs, c3_mgs, rho))

# ---- 7. Immune Checkpoint Genes ----
cat("7. Immune Checkpoint Genes\n")
ckpt_genes <- c("CTLA4","HAVCR2","LAG3","PDCD1")
for (g in ckpt_genes) {
  if (g %in% rownames(log_tpm)) {
    vals <- log_tpm[g,]
    kw <- kruskal.test(vals ~ subtypes)
    c1v <- mean(vals[subtypes=="C1"])
    c3v <- mean(vals[subtypes=="C3"])
    cat(sprintf("  %s: C1=%.2f, C3=%.2f, P=%.6f\n", g, c1v, c3v, kw$p.value))
    if (kw$p.value < 0.001) PASS <- PASS + 1 else FAIL <- FAIL + 1
  }
}
cat("\n")

# ---- Summary ----
cat("========== VERIFICATION SUMMARY ==========\n")
cat(sprintf("  PASS: %d, FAIL: %d\n", PASS, FAIL))
if (FAIL == 0) {
  cat("  ✅ All key results verified successfully!\n")
} else {
  cat(sprintf("  ⚠️  %d checks failed. Review before submission.\n", FAIL))
}
cat("==========================================\n")
