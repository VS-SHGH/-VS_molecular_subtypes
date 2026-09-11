# =============================================================================
# Revision Analysis 2: MGS External Validation + Clustering Sensitivity
# =============================================================================

library(dplyr)
library(princurve)
library(ConsensusClusterPlus)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"
set.seed(42)

# ---- Load data helper ----
load_external <- function(path) {
  dat <- read.csv(path, check.names = FALSE)
  genes <- as.character(dat[, 1])
  keep <- genes != "" & !is.na(genes)
  dat <- dat[keep, ]; genes <- genes[keep]
  dup <- duplicated(genes)
  dat <- dat[!dup, ]; genes <- genes[!dup]
  mat <- as.matrix(dat[, -1, drop = FALSE])
  rownames(mat) <- genes
  return(mat)
}

# ---- Load discovery data ----
tpm <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/clean_tpm.csv", row.names = 1, check.names = FALSE)
log_tpm <- log2(tpm + 1)
cat(sprintf("Discovery: %d genes x %d samples\n", nrow(log_tpm), ncol(log_tpm)))

# ---- Load external ----
gse141801 <- load_external("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE141801_Gugel/GSE141801_Expression_Log2.csv")
gse39645 <- load_external("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE39645_Torres/GSE39645_Expression_Log2.csv")
cat(sprintf("GSE141801: %d x %d\n", nrow(gse141801), ncol(gse141801)))
cat(sprintf("GSE39645: %d x %d\n", nrow(gse39645), ncol(gse39645)))

# ---- 1. Build MGS on discovery (3,000 MAD genes) ----
mads <- apply(log_tpm, 1, mad)
top_3000 <- names(sort(mads, decreasing = TRUE))[1:3000]

pca_d <- prcomp(t(log_tpm[top_3000, ]), center = TRUE, scale. = TRUE)
pc_d <- pca_d$x[, 1:5]
fit_d <- principal_curve(as.matrix(pc_d), stretch = 0)
mgs_d <- fit_d$lambda
mgs_d <- (mgs_d - min(mgs_d)) / (max(mgs_d) - min(mgs_d))
cat(sprintf("\nDiscovery MGS: [%.3f, %.3f]\n", min(mgs_d), max(mgs_d)))

# Load subtypes
cl <- readxl::read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]
s2s <- setNames(cl$Subtype, cl$SampleID)
mgs_by_subtype <- data.frame(Sample = colnames(log_tpm), MGS = mgs_d,
                              Subtype = s2s[colnames(log_tpm)])
cat("\nMGS by subtype (discovery):\n")
print(mgs_by_subtype %>% group_by(Subtype) %>%
  summarise(mean=mean(MGS), sd=sd(MGS), n=n()))

# ---- 2. External projection function ----
project_mgs <- function(ext_mat, log_tpm, top_genes, fit_obj, mgs_vec) {
  common <- intersect(rownames(log_tpm), rownames(ext_mat))
  top_avail <- intersect(top_genes, common)
  cat(sprintf("  Common genes: %d, top avail: %d\n", length(common), length(top_avail)))

  disc_sub <- log_tpm[top_avail, ]
  ext_sub <- ext_mat[top_avail, ]

  # Gene-wise means and SDs from discovery
  gmeans <- rowMeans(disc_sub)
  gsds <- apply(disc_sub, 1, sd)
  gsds[gsds == 0] <- 1e-10

  # PCA on discovery
  disc_scaled <- (disc_sub - gmeans) / gsds
  pca_m <- prcomp(t(disc_scaled), center = FALSE, scale. = FALSE)

  # Project external
  ext_scaled <- (ext_sub - gmeans) / gsds
  ext_pc <- as.matrix(t(ext_scaled)) %*% pca_m$rotation[, 1:min(5, ncol(pca_m$rotation))]

  # Nearest point on discovery curve
  curve_points <- fit_obj$s[fit_obj$ord, 1:min(5, ncol(ext_pc))]
  ext_mgs <- numeric(nrow(ext_pc))
  for (i in seq_len(nrow(ext_pc))) {
    dists <- sqrt(colSums((t(curve_points) - ext_pc[i, ])^2))
    ext_mgs[i] <- mgs_vec[fit_obj$ord[which.min(dists)]]
  }
  return(ext_mgs)
}

# ---- 3. Project GSE141801 ----
cat("\n--- GSE141801 projection ---\n")
mgs_141801 <- project_mgs(gse141801, log_tpm, top_3000, fit_d, mgs_d)
cat(sprintf("GSE141801 MGS: [%.3f, %.3f], mean=%.3f\n",
            min(mgs_141801), max(mgs_141801), mean(mgs_141801)))

# ---- 4. Project GSE39645 ----
cat("\n--- GSE39645 projection ---\n")
mgs_39645 <- project_mgs(gse39645, log_tpm, top_3000, fit_d, mgs_d)
cat(sprintf("GSE39645 MGS: [%.3f, %.3f], mean=%.3f\n",
            min(mgs_39645), max(mgs_39645), mean(mgs_39645)))

# ---- 5. Clustering sensitivity ----
cat("\n========== Clustering Sensitivity ==========\n")
gene_counts <- c(1000, 2000, 3000, 5000, 10000)
gene_counts <- gene_counts[gene_counts <= nrow(log_tpm)]

for (ng in gene_counts) {
  top_ng <- names(sort(mads, decreasing = TRUE))[1:ng]
  clust_input <- log_tpm[top_ng, ]
  clust_input <- sweep(clust_input, 1, apply(clust_input, 1, median, na.rm = TRUE))

  tryCatch({
    ccp <- ConsensusClusterPlus(as.matrix(clust_input), maxK = 6, reps = 80,
                                 pItem = 0.8, pFeature = 1,
                                 clusterAlg = "hc", distance = "pearson",
                                 seed = 42, plot = NULL)
    pac_vals <- sapply(2:6, function(k) {
      cm <- ccp[[k]]$consensusMatrix
      vals <- cm[lower.tri(cm)]
      mean(vals > 0.1 & vals < 0.9)
    })
    best_k <- which.min(pac_vals) + 1
    cat(sprintf("n_genes=%d: PAC(k=2-6) = %s | best k=%d (PAC=%.4f)\n",
                ng, paste(sprintf("%.4f", pac_vals), collapse=", "), best_k, min(pac_vals)))
  }, error = function(e) cat(sprintf("n_genes=%d: ERROR - %s\n", ng, e$message)))
}

# ---- 6. Gene overlap ----
cat("\n========== Gene Set Overlap ==========\n")
for (ng in c(1000, 2000, 3000, 5000)) {
  g_clust <- names(sort(mads, decreasing = TRUE))[1:ng]
  g_mgs <- names(sort(mads, decreasing = TRUE))[1:ng]
  ol <- length(intersect(g_clust, g_mgs))
  cat(sprintf("Both at %d genes: overlap = %d (%.1f%%)\n", ng, ol, 100*ol/ng))
}
# More relevant: 2000 for clustering vs 3000 for MGS
g2k <- names(sort(mads, decreasing = TRUE))[1:2000]
g3k <- names(sort(mads, decreasing = TRUE))[1:3000]
ol <- length(intersect(g2k, g3k))
cat(sprintf("Clustering(2k) vs MGS(3k): overlap = %d (%.1f%% of 2k)\n", ol, 100*ol/2000))

# ---- Save ----
saveRDS(list(mgs_discovery=mgs_d, mgs_gse141801=mgs_141801, mgs_gse39645=mgs_39645),
        file.path(work_dir, "mgs_results.rds"))
cat("\n========== Complete ==========\n")
