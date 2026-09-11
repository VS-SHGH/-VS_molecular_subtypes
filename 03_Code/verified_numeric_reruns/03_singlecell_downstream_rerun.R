# Discover Oncology VS: reproduce single-cell downstream results from the
# existing preprocessed Seurat object. This deliberately does not reconstruct
# FASTQ/10X preprocessing and does not create figures.

options(stringsAsFactors = FALSE)

source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
ROOT <- package_root
OUT <- file.path(AUDIT_ROOT, "singlecell_cellchat_rerun")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SC_QS <- SC_QS_INPUT
RDATA <- SERVER_RDATA_DIR
SC_CACHE <- file.path(RDATA, "Step06_sc_new_scored.rds")
MYELOID_CACHE <- file.path(RDATA, "Step07_all_myeloid.rds")
C3_CACHE <- file.path(RDATA, "Step07_C3high_myeloid.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(mclust)
})

message("Loading preprocessed Seurat object: ", SC_QS)
stopifnot(file.exists(SC_QS))
sc <- qread(SC_QS)
message("Loaded cells=", ncol(sc), "; genes=", nrow(sc),
        "; object_version=", as.character(sc@version))

# Match the manuscript code's annotation and module-score logic.
stopifnot("final_label" %in% colnames(sc@meta.data))
Idents(sc) <- "final_label"

deg_env <- new.env(parent = emptyenv())
load(file.path(RDATA, "Step3_DEA_Results.RData"), envir = deg_env)
stopifnot(all(c("deg_c1", "deg_c2", "deg_c3") %in% ls(deg_env)))

get_top_sig <- function(deg_df, n = 100) {
  x <- as.data.frame(deg_df)
  if (!"Gene" %in% colnames(x)) {
    gene_col <- grep("^gene$|^GeneName$|^SYMBOL$|^gene_name$",
                     colnames(x), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(gene_col)) x$Gene <- x[[gene_col]]
    else x$Gene <- rownames(x)
  }
  if (!"logFC" %in% colnames(x)) {
    lfc_col <- grep("log.*FC|avg_log2FC", colnames(x),
                    value = TRUE, ignore.case = TRUE)[1]
    if (is.na(lfc_col)) stop("No logFC column in DEG table")
    x$logFC <- x[[lfc_col]]
  }
  if (!"adj.P.Val" %in% colnames(x)) {
    p_col <- grep("adj\\.P\\.Val|padj|FDR|p_val_adj",
                  colnames(x), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(p_col)) stop("No adjusted p-value column in DEG table")
    x$adj.P.Val <- x[[p_col]]
  }
  x %>%
    filter(!is.na(adj.P.Val), !is.na(logFC),
           adj.P.Val < 0.05, logFC > 1) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = n) %>%
    pull(Gene) %>%
    unique()
}

gene_sets <- list(
  C1_Proliferative = get_top_sig(deg_env$deg_c1, 100),
  C2_Mesenchymal = get_top_sig(deg_env$deg_c2, 100),
  C3_Immune = get_top_sig(deg_env$deg_c3, 100)
)
sc_genes <- rownames(sc)
gene_sets_matched <- lapply(gene_sets, function(g) {
  idx <- match(toupper(g), toupper(sc_genes))
  unique(sc_genes[idx[!is.na(idx)]])
})
if (length(gene_sets_matched$C1_Proliferative) < 5) {
  rescue <- c("MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "MCM6", "CCNA2")
  idx <- match(toupper(rescue), toupper(sc_genes))
  gene_sets_matched$C1_Proliferative <- unique(c(
    gene_sets_matched$C1_Proliferative, sc_genes[idx[!is.na(idx)]]
  ))
}

write.csv(
  data.frame(
    signature = names(gene_sets),
    raw_n = vapply(gene_sets, length, integer(1)),
    matched_n = vapply(gene_sets_matched, length, integer(1)),
    stringsAsFactors = FALSE
  ),
  file.path(OUT, "module_signature_gene_counts.csv"), row.names = FALSE
)

message("Running Seurat AddModuleScore with seed=42...")
old_cols <- grep("Score_C|Bulk_Score", colnames(sc@meta.data), value = TRUE)
if (length(old_cols)) {
  sc@meta.data <- sc@meta.data[, setdiff(colnames(sc@meta.data), old_cols), drop = FALSE]
}
sc <- AddModuleScore(sc, features = gene_sets_matched,
                     name = "Bulk_Score", seed = 42)
sc$Score_C1 <- sc$Bulk_Score1
sc$Score_C2 <- sc$Bulk_Score2
sc$Score_C3 <- sc$Bulk_Score3
score_mat <- sc@meta.data[, c("Score_C1", "Score_C2", "Score_C3"), drop = FALSE]
sc$Bulk_Class <- factor(
  c("C1", "C2", "C3")[max.col(score_mat, ties.method = "first")],
  levels = c("C1", "C2", "C3")
)

# UMAP is stored as metadata in the downloaded object; reconstruct the same
# reduction used by the manuscript script for identity/coordinate checks.
if (!"umap" %in% Reductions(sc) && all(c("UMAP_1", "UMAP_2") %in% colnames(sc@meta.data))) {
  sc[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(sc@meta.data[, c("UMAP_1", "UMAP_2")]),
    key = "UMAP_", assay = DefaultAssay(sc)
  )
}

score_summary <- sc@meta.data %>%
  mutate(cell_type = as.character(Idents(sc))) %>%
  group_by(cell_type) %>%
  summarise(
    n_cells = n(),
    Score_C1_mean = mean(Score_C1),
    Score_C2_mean = mean(Score_C2),
    Score_C3_mean = mean(Score_C3),
    .groups = "drop"
  )
write.csv(score_summary, file.path(OUT, "singlecell_score_by_celltype.csv"), row.names = FALSE)
write.csv(as.data.frame(table(sc$Bulk_Class)),
          file.path(OUT, "singlecell_bulk_class_counts.csv"), row.names = FALSE)

# Compare the fresh score calculation with the saved Step06 object.
if (file.exists(SC_CACHE)) {
  message("Comparing module scores with saved Step06 cache...")
  cached <- readRDS(SC_CACHE)
  common_cells <- intersect(colnames(sc), colnames(cached))
  score_cols <- c("Score_C1", "Score_C2", "Score_C3")
  common_cols <- intersect(score_cols, intersect(colnames(sc@meta.data), colnames(cached@meta.data)))
  fresh <- as.matrix(sc@meta.data[common_cells, common_cols, drop = FALSE])
  old <- as.matrix(cached@meta.data[common_cells, common_cols, drop = FALSE])
  max_diff <- max(abs(fresh - old), na.rm = TRUE)
  class_agree <- mean(as.character(sc$Bulk_Class[common_cells]) ==
                        as.character(cached$Bulk_Class[common_cells]))
  write.csv(data.frame(
    metric = c("common_cells", "max_abs_score_difference", "Bulk_Class_agreement"),
    value = c(length(common_cells), max_diff, class_agree)
  ), file.path(OUT, "singlecell_score_cache_comparison.csv"), row.names = FALSE)
  rm(cached, fresh, old)
  gc(verbose = FALSE)
}

# Reproduce the C3-high myeloid GMM split used by the manuscript code.
myeloid_levels <- grep("myeloid|macrophage|monocyte|dc|dendritic",
                       levels(Idents(sc)), value = TRUE, ignore.case = TRUE)
if (!length(myeloid_levels)) stop("No myeloid identity found")
sc_myeloid <- subset(sc, idents = myeloid_levels)
set.seed(123)  # fixed in the original manuscript-era GMM code
gmm_fit <- Mclust(sc_myeloid$Score_C3, G = 2, verbose = FALSE)
means <- sort(as.numeric(gmm_fit$parameters$mean))
vars <- as.numeric(gmm_fit$parameters$variance$sigmasq)
props <- as.numeric(gmm_fit$parameters$pro)
if (length(vars) == 1) vars <- rep(vars, 2)
if (length(props) == 1) props <- rep(props, 2)
ord <- order(as.numeric(gmm_fit$parameters$mean))
vars <- vars[ord]
props <- props[ord]
cutoff <- mean(means)
sc_myeloid$C3_status <- ifelse(sc_myeloid$Score_C3 > cutoff, "C3-high", "C3-low")
sc_c3high <- subset(sc_myeloid, subset = C3_status == "C3-high")

write.csv(data.frame(
  metric = c("myeloid_identity_levels", "myeloid_cells", "C3high_cells", "C3low_cells",
             "GMM_low_mean", "GMM_high_mean", "GMM_cutoff", "GMM_low_variance",
             "GMM_high_variance", "GMM_low_proportion", "GMM_high_proportion"),
  value = c(paste(myeloid_levels, collapse = "|"), ncol(sc_myeloid), ncol(sc_c3high),
            ncol(sc_myeloid) - ncol(sc_c3high), means[1], means[2], cutoff,
            vars[1], vars[2], props[1], props[2])
), file.path(OUT, "c3high_myeloid_gmm_summary.csv"), row.names = FALSE)

gmm_cell_labels <- data.frame(
  cell = colnames(sc_myeloid),
  Score_C3 = sc_myeloid$Score_C3,
  C3_status = as.character(sc_myeloid$C3_status),
  stringsAsFactors = FALSE
)
write.csv(gmm_cell_labels, file.path(OUT, "c3high_myeloid_cell_labels.csv"), row.names = FALSE)

if (file.exists(MYELOID_CACHE)) {
  cached_m <- readRDS(MYELOID_CACHE)
  common <- intersect(colnames(sc_myeloid), colnames(cached_m))
  cache_status <- if ("C3_status" %in% colnames(cached_m@meta.data)) {
    as.character(cached_m@meta.data[common, "C3_status"])
  } else rep(NA_character_, length(common))
  fresh_status <- as.character(sc_myeloid@meta.data[common, "C3_status"])
  write.csv(data.frame(
    metric = c("common_myeloid_cells", "C3_status_agreement"),
    value = c(length(common), mean(fresh_status == cache_status, na.rm = TRUE))
  ), file.path(OUT, "c3high_myeloid_cache_comparison.csv"), row.names = FALSE)
  rm(cached_m)
  gc(verbose = FALSE)
}

# Save only metadata and numeric objects for audit; the original and server
# cache objects remain untouched.
saveRDS(sc@meta.data, file.path(OUT, "singlecell_scored_metadata.rds"))
saveRDS(sc_myeloid@meta.data, file.path(OUT, "myeloid_scored_metadata.rds"))
saveRDS(sc_c3high@meta.data, file.path(OUT, "c3high_myeloid_scored_metadata.rds"))

write.csv(data.frame(
  metric = c("cells", "genes", "cell_types", "umap_reduction_present"),
  value = c(ncol(sc), nrow(sc), length(unique(Idents(sc))), "umap" %in% Reductions(sc))
), file.path(OUT, "singlecell_object_summary.csv"), row.names = FALSE)

message("Single-cell downstream rerun finished. Output: ", OUT)
