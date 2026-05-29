# =============================================================================
# Revision: Unbiased MGS - direction defined by biology, NOT by subtype labels
#
# Key change: Instead of anchoring MGS=0 to C1 and MGS=1 to C3,
# we anchor based on (1) Schwann cell genes → "differentiated" end
#                    (2) Immune genes → "immune-infiltrated" end
# Then check: do C1/C2/C3 NATURALLY align along this unbiased axis?
# =============================================================================

library(dplyr)
library(readxl)
library(princurve)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"

# ---- Load data ----
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

tpm_mat <- as.matrix(tpm_raw[, our_samples])
rownames(tpm_mat) <- genes
log_tpm <- log2(tpm_mat + 1)

cl <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]
sample_ids_clean <- gsub("^tpm\\.", "", our_samples)
subtypes <- setNames(cl$Subtype, cl$SampleID)[sample_ids_clean]

# =============================================================================
# APPROACH:
# 1. PCA on top variable genes → get PC1 (24% variance, captures main biology)
# 2. Principal curve → get arc-length λ (NO subtype-based anchoring)
# 3. Check PC1 loading: which genes are at PC1+ end? Which at PC1-?
# 4. Orient MGS based on EXTERNAL biology:
#    - Schwann cell differentiation score → differentiated end (MGS=0)
#    - Immune score → immune-infiltrated end (MGS=1)
# 5. THEN check if C1/C2/C3 align along this unbiased axis
# =============================================================================

cat("========== Unbiased MGS Construction ==========\n")

# Step 1: Select top variable genes (subtype-agnostic)
mads <- apply(log_tpm, 1, mad)
top_3000 <- names(sort(mads, decreasing = TRUE))[1:3000]

# Step 2: PCA
pca <- prcomp(t(log_tpm[top_3000, ]), center = TRUE, scale. = TRUE)
cat(sprintf("PC1: %.1f%%, PC2: %.1f%%\n",
            100*summary(pca)$importance[2,1], 100*summary(pca)$importance[2,2]))

# Step 3: Analyze PC1 loadings - what IS this axis biologically?
pc1_load <- pca$rotation[, 1]
# Top positive and negative genes
pc1_pos <- names(sort(pc1_load, decreasing = TRUE))[1:50]
pc1_neg <- names(sort(pc1_load, decreasing = FALSE))[1:50]

cat("\n--- PC1 Biological Interpretation (gene loadings) ---\n")
cat("PC1+ end (top 50 genes):\n")
cat(paste(pc1_pos[1:25], collapse=", "), "\n")
cat("PC1- end (top 50 genes):\n")
cat(paste(pc1_neg[1:25], collapse=", "), "\n")

# Step 4: Compute EXTERNAL biological scores for anchoring
# Schwann cell differentiation score
sc_genes <- c("SOX10", "PLP1", "EGR2", "MPZ", "MBP", "PMP22", "MAL", "PRX", "CDH19")
sc_avail <- intersect(sc_genes, rownames(log_tpm))
sc_score <- colMeans(log_tpm[sc_avail, ])  # higher = more differentiated

# Immune score
immune_genes <- c("CD2","CD3D","CD3E","CD8A","GZMA","PRF1","CTLA4","HAVCR2",
                   "LAG3","PDCD1","CD68","CD163","CD74","HLA-DRA","CCL5","CXCL9","CXCL10")
im_avail <- intersect(immune_genes, rownames(log_tpm))
immune_score <- colMeans(log_tpm[im_avail, ])  # higher = more immune

# ECM score (for C2 direction check)
ecm_genes <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","COL6A3",
               "FN1","LOX","LOXL1","SPARC","THBS2","TGFB1","CTGF","POSTN")
ecm_avail <- intersect(ecm_genes, rownames(log_tpm))
ecm_score <- colMeans(log_tpm[ecm_avail, ])

# Step 5: Principal curve (NO anchoring to subtypes)
fit <- principal_curve(as.matrix(pca$x[, 1:5]), stretch = 0)
mgs_raw <- fit$lambda  # arc-length ordering, NOT anchored to C1/C3

# The raw MGS is just an ordering. Which direction should be 0 vs 1?
# Check: which end has higher SC score? That should be MGS=0 (differentiated)
# Which end has higher immune score? That should be MGS=1

cat("\n--- Anchoring MGS by Biology ---\n")

# Method 1: Orient so MGS correlates positively with immune score
# (i.e., "de-differentiated" = "immune-infiltrated" end is MGS=1)
cor_raw_immune <- cor(mgs_raw, immune_score, method = "spearman")
cat(sprintf("Raw MGS vs Immune Score: rho=%.3f\n", cor_raw_immune))

# If correlation is negative, flip MGS direction
if (cor_raw_immune < 0) {
  mgs_immune_anchored <- 1 - (mgs_raw - min(mgs_raw)) / (max(mgs_raw) - min(mgs_raw))
  cat("Flipping MGS to align with immune score...\n")
} else {
  mgs_immune_anchored <- (mgs_raw - min(mgs_raw)) / (max(mgs_raw) - min(mgs_raw))
  cat("MGS already aligned with immune score.\n")
}

# Verify
cat(sprintf("Anchored MGS vs Immune: rho=%.3f\n",
            cor(mgs_immune_anchored, immune_score, method = "spearman")))
cat(sprintf("Anchored MGS vs SC Score: rho=%.3f (should be negative)\n",
            cor(mgs_immune_anchored, sc_score, method = "spearman")))
cat(sprintf("Anchored MGS vs ECM Score: rho=%.3f\n",
            cor(mgs_immune_anchored, ecm_score, method = "spearman")))

# Step 6: THE KEY TEST - do subtypes NATURALLY align along this unbiased axis?
cat("\n========== THE KEY TEST: Do subtypes align along unbiased MGS? ==========\n")

mgs_unbiased <- setNames(mgs_immune_anchored, our_samples)
mgs_df <- data.frame(
  Sample = our_samples,
  MGS_unbiased = mgs_unbiased,
  Subtype = subtypes,
  SC_Score = sc_score,
  Immune_Score = immune_score,
  ECM_Score = ecm_score,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

cat("\nMGS (immune-anchored) by Subtype:\n")
mgs_summary <- mgs_df %>% group_by(Subtype) %>%
  summarise(
    MGS = mean(MGS_unbiased),
    SD = sd(MGS_unbiased),
    SC = mean(SC_Score),
    Immune = mean(Immune_Score),
    ECM = mean(ECM_Score)
  )
print(mgs_summary)

kw_mgs <- kruskal.test(MGS_unbiased ~ Subtype, data = mgs_df)
cat(sprintf("\nKruskal-Wallis (MGS by Subtype): P = %.2e\n", kw_mgs$p.value))

# Pairwise comparisons
cat("\nPairwise Dunn's test:\n")
# Manual pairwise
for (s1 in c("C1","C2")) {
  for (s2 in c("C2","C3")) {
    if (s1 < s2) {
      v1 <- mgs_df$MGS_unbiased[mgs_df$Subtype == s1]
      v2 <- mgs_df$MGS_unbiased[mgs_df$Subtype == s2]
      wt <- wilcox.test(v1, v2)
      cat(sprintf("  %s vs %s: W=%.1f, P=%.4f, medians: %.3f vs %.3f\n",
                  s1, s2, wt$statistic, wt$p.value, median(v1), median(v2)))
    }
  }
}

# Step 7: Compare with ORIGINAL approach (subtype-anchored)
cat("\n========== Comparison: Unbiased vs Original (subtype-anchored) MGS ==========\n")

# Original approach: anchor C1=0, C3=1
c1_samples <- our_samples[subtypes == "C1"]
c3_samples <- our_samples[subtypes == "C3"]
# This was done by ordering the principal curve lambda to match C1→C3
# Simulate: find which end of PC1 the C1 and C3 centroids are at
c1_centroid_pc1 <- mean(pca$x[subtypes == "C1", 1])
c3_centroid_pc1 <- mean(pca$x[subtypes == "C3", 1])
cat(sprintf("C1 centroid on PC1: %.3f\n", c1_centroid_pc1))
cat(sprintf("C3 centroid on PC1: %.3f\n", c3_centroid_pc1))
cat(sprintf("C2 centroid on PC1: %.3f\n", mean(pca$x[subtypes == "C2", 1])))

# The arc-length ordering from principal curve
# Does it NATURALLY order C1 → C2 → C3?
order_df <- data.frame(
  Sample = our_samples,
  Lambda_raw = fit$lambda,
  PC1 = pca$x[, 1],
  Subtype = subtypes
)
cat("\nNatural arc-length (lambda) by subtype (BEFORE any anchoring):\n")
print(order_df %>% group_by(Subtype) %>%
  summarise(Lambda = mean(Lambda_raw), PC1 = mean(PC1)))

# =============================================================================
# Answer: Does the gradient EXIST independently of subtype labeling?
# YES, if:
#   1. PC1 captures a coherent biological axis (immune vs Schwann cell)
#   2. The principal curve ordering is consistent with PC1
#   3. C1/C2/C3 naturally order along this axis (NOT manually forced)
# =============================================================================

cat("\n========== Verdict ==========\n")

# Rank correlation between unbiased MGS and subtype ordering (C1=1, C2=2, C3=3)
subtype_num <- ifelse(subtypes == "C1", 1, ifelse(subtypes == "C2", 2, 3))
rank_cor <- cor(mgs_unbiased, subtype_num, method = "spearman")
cat(sprintf("Spearman rho (unbiased MGS vs subtype order 1-2-3): %.3f\n", rank_cor))

# Is the ordering C1 < C2 < C3 preserved without anchoring?
# Check mean MGS values
mgs_order <- mgs_df %>% group_by(Subtype) %>% summarise(m = mean(MGS_unbiased)) %>% arrange(m)
cat(sprintf("Natural MGS order: %s\n", paste(mgs_order$Subtype, collapse=" < ")))

if (all(mgs_order$Subtype == c("C1", "C2", "C3"))) {
  cat("\n✅ C1-C2-C3 ordering is PRESERVED in unbiased MGS.\n")
  cat("   The gradient exists independently of subtype anchoring.\n")
  cat("   Anchoring by immune score naturally recapitulates the subtype order.\n")
} else {
  cat(sprintf("\n⚠️  Natural order is %s, NOT C1-C2-C3.\n",
              paste(mgs_order$Subtype, collapse=" < ")))
  cat("   The C1-C2-C3 ordering depends on how MGS is anchored.\n")
}

# Also check: if we anchor by SC score instead (differentiated=0)
sc_cor_raw <- cor(mgs_raw, sc_score, method = "spearman")
if (sc_cor_raw < 0) {
  mgs_sc_anchored <- 1 - (mgs_raw - min(mgs_raw)) / (max(mgs_raw) - min(mgs_raw))
} else {
  mgs_sc_anchored <- (mgs_raw - min(mgs_raw)) / (max(mgs_raw) - min(mgs_raw))
}
mgs_sc_df <- data.frame(Subtype = subtypes, MGS_SC = mgs_sc_anchored)
cat("\nSC-anchored MGS order:\n")
sc_order <- mgs_sc_df %>% group_by(Subtype) %>% summarise(m = mean(MGS_SC)) %>% arrange(m)
cat(sprintf("  %s\n", paste(sc_order$Subtype, collapse=" < ")))

# ---- Save ----
saveRDS(list(
  mgs_unbiased = mgs_df,
  pc1_loadings = pc1_load,
  sc_score = sc_score,
  immune_score = immune_score,
  ecm_score = ecm_score
), file.path(work_dir, "mgs_unbiased_results.rds"))

cat("\n========== Complete ==========\n")
