# Second-pass numerical rerun of the primary bulk clustering/PCA analysis.
# This deliberately produces no figures and never writes to the server RData.

suppressPackageStartupMessages({
  library(readxl)
  library(ConsensusClusterPlus)
})

source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
PASS_ROOT <- AUDIT_ROOT
OUT_DIR <- file.path(PASS_ROOT, "primary_clustering")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

step1_env <- new.env(parent = emptyenv())
load(file.path(SERVER_RDATA_DIR, "Step1_Clean_Data_VS3.RData"), envir = step1_env)
step2_env <- new.env(parent = emptyenv())
load(file.path(SERVER_RDATA_DIR, "Step2_Clustering.RData"), envir = step2_env)

reference_meta <- step2_env$meta_step2
reference <- setNames(as.character(reference_meta$Subtype), reference_meta$SampleID)

prepare_from_matrix <- function(tpm, source_name) {
  stopifnot(!is.null(rownames(tpm)), !is.null(colnames(tpm)))
  colnames(tpm) <- sub("^tpm\\.", "", colnames(tpm))
  tpm <- tpm[, intersect(colnames(tpm), names(reference)), drop = FALSE]
  sample_ids <- colnames(tpm)
  keep_samples <- !is.na(reference[sample_ids])
  tpm <- tpm[, keep_samples, drop = FALSE]
  sample_ids <- sub("^tpm\\.", "", colnames(tpm))
  if (!identical(sample_ids, names(reference)[match(sample_ids, names(reference))])) {
    stop("Sample alignment failed for ", source_name)
  }
  log_tpm <- log2(tpm + 1)
  mads <- apply(log_tpm, 1, mad, na.rm = TRUE)
  top_n <- min(5000L, nrow(log_tpm))
  top_genes <- names(sort(mads, decreasing = TRUE))[seq_len(top_n)]
  clust_input <- log_tpm[top_genes, , drop = FALSE]
  clust_input <- sweep(clust_input, 1, apply(clust_input, 1, median, na.rm = TRUE))

  set.seed(123456)
  ccp <- ConsensusClusterPlus(
    as.matrix(clust_input), maxK = 6, reps = 1000,
    pItem = 0.8, pFeature = 1, clusterAlg = "hc",
    distance = "pearson", seed = 123456,
    title = file.path(OUT_DIR, paste0(source_name, "_ccp")), plot = NULL
  )
  pac <- sapply(2:6, function(k) {
    vals <- ccp[[k]]$consensusMatrix[lower.tri(ccp[[k]]$consensusMatrix)]
    mean(vals > 0.1 & vals < 0.9)
  })
  names(pac) <- paste0("k=", 2:6)

  labels <- as.integer(ccp[[3]]$consensusClass)
  observed <- reference[sample_ids]
  perms <- rbind(
    c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
    c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
  )
  agreements <- apply(perms, 1, function(p) mean(p[labels] == as.integer(sub("C", "", observed))))
  best <- which.max(agreements)
  mapped <- paste0("C", perms[best, labels])

  pca <- prcomp(t(log_tpm[top_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
  pca_var <- 100 * summary(pca)$importance[2, 1:5]

  summary_row <- data.frame(
    source = source_name,
    genes = nrow(tpm), samples = ncol(tpm),
    top_genes = top_n,
    pc1 = pca_var[1], pc2 = pca_var[2], pc3 = pca_var[3],
    pc4 = pca_var[4], pc5 = pca_var[5],
    pac_k2 = pac["k=2"], pac_k3 = pac["k=3"], pac_k4 = pac["k=4"],
    pac_k5 = pac["k=5"], pac_k6 = pac["k=6"],
    pac_min_k = as.integer(sub("k=", "", names(which.min(pac)))),
    rerun_k3_sizes = paste(as.integer(table(labels)), collapse = "/"),
    locked_k3_sizes = paste(as.integer(table(observed)), collapse = "/"),
    best_label_agreement = max(agreements),
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(SampleID = sample_ids, locked = observed,
               rerun_cluster = labels, rerun_mapped = mapped),
    file.path(OUT_DIR, paste0(source_name, "_sample_labels.csv")), row.names = FALSE
  )
  saveRDS(list(source = source_name, pac = pac, pca_var = pca_var,
               ccp = ccp, sample_labels = summary_row),
          file.path(OUT_DIR, paste0(source_name, "_rerun.rds")))
  summary_row
}

step1_tpm <- step1_env$clean_tpm
step1_result <- prepare_from_matrix(step1_tpm, "Step1_clean_tpm")

tpm_xlsx <- as.data.frame(read_excel(BULK_TPM_XLSX), check.names = FALSE)
genes <- as.character(tpm_xlsx$gene_name)
keep <- !is.na(genes) & nzchar(genes) & !duplicated(genes)
tpm_xlsx_mat <- as.matrix(tpm_xlsx[keep, grep("^tpm\\.", names(tpm_xlsx)), drop = FALSE])
storage.mode(tpm_xlsx_mat) <- "numeric"
rownames(tpm_xlsx_mat) <- genes[keep]
xlsx_result <- prepare_from_matrix(tpm_xlsx_mat, "TPM_xlsx_gene_named")

summary <- rbind(step1_result, xlsx_result)
write.csv(summary, file.path(OUT_DIR, "primary_clustering_second_pass_summary.csv"), row.names = FALSE)
writeLines(capture.output({
  cat("Second-pass primary bulk clustering/PCA rerun\n")
  cat("R: ", R.version.string, "\n\n", sep = "")
  print(summary, row.names = FALSE)
}), file.path(OUT_DIR, "primary_clustering_second_pass_summary.txt"))
print(summary, row.names = FALSE)
