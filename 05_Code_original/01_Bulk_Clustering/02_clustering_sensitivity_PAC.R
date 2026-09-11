# =============================================================================
# Revision Analysis 4: Clustering Sensitivity + Purity Confounding
# =============================================================================

library(dplyr)
library(ConsensusClusterPlus)
library(princurve)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"
set.seed(42)

# ---- Load data ----
tpm_raw <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/clean_tpm.csv", check.names = FALSE)
# This file has sample IDs as column headers but NO gene symbol column
# Rows represent genes (indices only), cols represent samples
# Extract sample IDs from headers, data is the matrix
sample_ids <- colnames(tpm_raw)
tpm <- as.matrix(tpm_raw)
colnames(tpm) <- sample_ids
rownames(tpm) <- paste0("G", 1:nrow(tpm))  # Assign gene indices
log_tpm <- log2(tpm + 1)
cat(sprintf("Data: %d genes x %d samples\n", nrow(log_tpm), ncol(log_tpm)))

mads <- apply(log_tpm, 1, mad)

# ---- 1. Clustering sensitivity with different gene counts ----
cat("\n========== Clustering Sensitivity (gene counts) ==========\n")

gene_counts <- c(1000, 2000, 3000, 5000, 10000)
gene_counts <- gene_counts[gene_counts <= nrow(log_tpm)]
max_k <- 6

pac_matrix <- matrix(NA, nrow = length(gene_counts), ncol = max_k - 1,
                      dimnames = list(paste0("n", gene_counts), paste0("k", 2:max_k)))

for (i in seq_along(gene_counts)) {
  ng <- gene_counts[i]
  top_ng <- names(sort(mads, decreasing = TRUE))[1:ng]
  clust_input <- log_tpm[top_ng, ]
  clust_input <- sweep(clust_input, 1, apply(clust_input, 1, median, na.rm = TRUE))

  for (kv in 2:max_k) {
    tryCatch({
      ccp <- ConsensusClusterPlus(as.matrix(clust_input), maxK = kv, reps = 100,
                                   pItem = 0.8, pFeature = 1,
                                   clusterAlg = "hc", distance = "pearson",
                                   seed = 42, plot = NULL)
      cm <- ccp[[kv]]$consensusMatrix
      vals <- cm[lower.tri(cm)]
      pac_matrix[i, kv-1] <- mean(vals > 0.1 & vals < 0.9)
    }, error = function(e) {
      cat(sprintf("  n=%d, k=%d: ERROR\n", ng, kv))
    })
  }
  best_k <- which.min(pac_matrix[i, ]) + 1
  cat(sprintf("n_genes=%d: best_k=%d, PAC values: %s\n",
              ng, best_k, paste(sprintf("%.4f", pac_matrix[i, ]), collapse = ", ")))
}

# ---- 2. Subsampling stability ----
cat("\n========== Subsampling Stability (80% samples, 50 iterations) ==========\n")

top_2000 <- names(sort(mads, decreasing = TRUE))[1:2000]
n_iter <- 50
stability <- matrix(0, nrow = ncol(log_tpm), ncol = ncol(log_tpm))
n_subsample <- round(0.8 * ncol(log_tpm))

for (iter in 1:n_iter) {
  idx <- sample(1:ncol(log_tpm), n_subsample)
  sub_data <- log_tpm[top_2000, idx]
  sub_data <- sweep(sub_data, 1, apply(sub_data, 1, median, na.rm = TRUE))

  tryCatch({
    ccp <- ConsensusClusterPlus(as.matrix(sub_data), maxK = 3, reps = 50,
                                 pItem = 0.8, pFeature = 1,
                                 clusterAlg = "hc", distance = "pearson",
                                 seed = iter, plot = NULL)
    labels <- ccp[[3]]$consensusClass
    for (i in 1:(length(labels)-1)) {
      for (j in (i+1):length(labels)) {
        if (labels[i] == labels[j]) {
          stability[idx[i], idx[j]] <- stability[idx[i], idx[j]] + 1
          stability[idx[j], idx[i]] <- stability[idx[j], idx[i]] + 1
        }
      }
    }
  }, error = function(e) {})
}
cat(sprintf("Co-clustering frequency: mean=%.3f, median=%.3f\n",
            mean(stability[upper.tri(stability)])/n_iter,
            median(stability[upper.tri(stability)])/n_iter))

# ---- 3. Gradient-clustering overlap quantification ----
cat("\n========== Gene Set Overlap Analysis ==========\n")

for (clust_n in c(1000, 2000, 3000)) {
  for (mgs_n in c(1000, 2000, 3000, 5000)) {
    g_clust <- names(sort(mads, decreasing = TRUE))[1:clust_n]
    g_mgs <- names(sort(mads, decreasing = TRUE))[1:mgs_n]
    ol <- length(intersect(g_clust, g_mgs))
    cat(sprintf("Clust(%d) vs MGS(%d): overlap=%d (%.1f%% of clust)\n",
                clust_n, mgs_n, ol, 100*ol/clust_n))
  }
}

# ---- 4. PCA variance explained ----
cat("\n========== PCA Variance Explained ==========\n")
pca <- prcomp(t(log_tpm[top_2000, ]), center = TRUE, scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:10]
cat("Top 10 PCs:\n")
for (i in 1:10) {
  cat(sprintf("  PC%d: %.1f%%\n", i, 100*var_explained[i]))
}
cat(sprintf("Cumulative PC1-5: %.1f%%\n", 100*sum(var_explained[1:5])))

# ---- Save ----
saveRDS(list(pac_matrix = pac_matrix, var_explained = var_explained),
        file.path(work_dir, "clustering_sensitivity.rds"))
cat("\n========== Complete ==========\n")
