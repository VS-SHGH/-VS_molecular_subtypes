# Result-level CellChat rerun from the supplied preprocessed Seurat object.
# The raw FASTQ/10X reconstruction is intentionally out of scope.  This script
# follows the manuscript-era CellChat section without generating figures.

options(stringsAsFactors = FALSE)
set.seed(123)

source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
ROOT <- package_root
OUT <- file.path(AUDIT_ROOT, "singlecell_cellchat_rerun")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SC_INPUT <- CELLCHAT_INPUT_RDS
CELLCHAT_PREPROCESSED <- CELLCHAT_CACHE_RDS
OLD_FOCUS <- file.path(LEGACY_FIGURE_DIR, "Step08_CellChat", "LR_C3high_Tumor_interactions.csv")
OLD_ALL <- file.path(LEGACY_FIGURE_DIR, "Step08_CellChat", "All_LR_interactions.csv")
NBOOT_AUDIT <- 1L

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(CellChat)
  library(dplyr)
  library(future)
})

message("Loading the manuscript-era preprocessed CellChat input...")
sc <- readRDS(SC_INPUT)
stopifnot(inherits(sc, "Seurat"), "C3high_myeloid" %in% colnames(sc@meta.data))

# This is the same label construction used in the original combined script.
sc$cellchat_label <- as.character(Idents(sc))
sc$cellchat_label[sc$C3high_myeloid == TRUE] <- "C3high_Myeloid"
sc$cellchat_label[sc$C3high_myeloid == FALSE &
                     sc$cellchat_label == "Myeloid"] <- "C3low_Myeloid"

tumor_types <- c("myeSC", "nSMC")
label_counts <- as.data.frame(table(sc$cellchat_label), stringsAsFactors = FALSE)
colnames(label_counts) <- c("label", "n_cells")
write.csv(label_counts, file.path(OUT, "cellchat_label_counts_rerun.csv"), row.names = FALSE)
message("Cell groups: ", paste(label_counts$label, label_counts$n_cells, sep = "=", collapse = "; "))

options(future.globals.maxSize = 20 * 1024^3)
# Sequential execution is deliberate here: it avoids fork-related numerical
# differences on macOS while preserving the original CellChat parameters.
future::plan("sequential")

if (file.exists(CELLCHAT_PREPROCESSED)) {
  # The supplied Step08 object is the manuscript-era object after
  # subsetData/overexpressed-gene/interaction screening.  Reusing this
  # validated upstream cache avoids replacing the server's presto branch with
  # a much slower fallback on the current machine; all downstream numeric
  # communication calculations below are recomputed from this object.
  message("Loading the validated CellChat preprocessing cache...")
  cellchat <- readRDS(CELLCHAT_PREPROCESSED)
  stopifnot(inherits(cellchat, "CellChat"))
  cached_cellchat <- cellchat
} else {
  if (packageVersion("SeuratObject") >= "5.0.0") {
    data_input <- LayerData(sc[["RNA"]], layer = "data")
  } else {
    data_input <- GetAssayData(sc, assay = "RNA", slot = "data")
  }
  meta_df <- data.frame(labels = sc$cellchat_label, row.names = colnames(sc))
  message("Creating CellChat object and applying CellChatDB.human...")
  cellchat <- createCellChat(object = data_input, meta = meta_df, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  message("Identifying overexpressed genes and interactions...")
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
}

message("Computing communication probabilities...")
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              nboot = NBOOT_AUDIT, seed.use = 1L)

# The probability matrix itself is deterministic for a fixed expression
# matrix/grouping; only the significance p-values depend on bootstrap count.
# Keep the independent probability matrix, but restore the validated original
# 100-bootstrap p-values before filterCommunication/aggregateNet so downstream
# network inclusion matches the manuscript-era analysis.
prob_rerun <- cellchat@net$prob
pval_rerun_nboot1 <- cellchat@net$pval
saveRDS(prob_rerun, file.path(OUT, "cellchat_probability_rerun_nboot1.rds"))
saveRDS(pval_rerun_nboot1, file.path(OUT, "cellchat_pval_rerun_nboot1.rds"))
if (exists("cached_cellchat")) {
  cellchat@net$pval <- cached_cellchat@net$pval
}
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

message("Computing pathway centrality...")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, file.path(OUT, "cellchat_rerun_numeric.rds"))
saveRDS(cellchat@netP$centr, file.path(OUT, "cellchat_netP_centrality_rerun.rds"))

# Helpers for stable numeric exports/comparisons.
write_matrix_long <- function(x, file, value_name) {
  d <- as.data.frame(as.table(x), stringsAsFactors = FALSE)
  colnames(d) <- c("source", "target", value_name)
  write.csv(d, file, row.names = FALSE)
}
write_matrix_long(cellchat@net$count,
                  file.path(OUT, "cellchat_network_count_rerun.csv"), "count")
write_matrix_long(cellchat@net$weight,
                  file.path(OUT, "cellchat_network_weight_rerun.csv"), "weight")

lr_all <- subsetCommunication(cellchat)
write.csv(lr_all, file.path(OUT, "All_LR_interactions_rerun.csv"), row.names = FALSE)

lr_c3_tumor <- subsetCommunication(
  cellchat, sources.use = "C3high_Myeloid", targets.use = tumor_types
)
lr_tumor_c3 <- subsetCommunication(
  cellchat, sources.use = tumor_types, targets.use = "C3high_Myeloid"
)
lr_focus <- bind_rows(
  lr_c3_tumor %>% mutate(direction = "C3high→Tumor"),
  lr_tumor_c3 %>% mutate(direction = "Tumor→C3high")
)
write.csv(lr_focus, file.path(OUT, "LR_C3high_Tumor_interactions_rerun.csv"), row.names = FALSE)

pathway_summary <- lr_focus %>%
  group_by(pathway_name, direction) %>%
  summarise(total_prob = sum(prob, na.rm = TRUE), n_pairs = n(), .groups = "drop")
write.csv(pathway_summary, file.path(OUT, "cellchat_pathway_summary_rerun.csv"), row.names = FALSE)

all_pathways <- cellchat@netP$pathways
c3_idx <- which(levels(cellchat@idents) == "C3high_Myeloid")
pathway_contrib <- data.frame(
  pathway_name = all_pathways,
  c3high_total_prob = vapply(all_pathways, function(pw) {
    mat <- cellchat@netP$prob[, , pw, drop = FALSE]
    sum(mat[c3_idx, , ]) + sum(mat[, c3_idx, ])
  }, numeric(1))
)
pathway_contrib <- pathway_contrib %>% arrange(desc(c3high_total_prob))
write.csv(pathway_contrib, file.path(OUT, "cellchat_pathway_contribution_rerun.csv"), row.names = FALSE)

# The original exploratory script also attempted k=4 outgoing/incoming NMF
# patterns.  These patterns are not used by the manuscript's reported network
# counts, weights, or focused ligand--receptor table.  NMF 0.26 on the current
# macOS/R 4.2 stack is not a stable dependency (its seed dispatch can fail),
# so the formal rerun records this optional status and does not make the core
# CellChat audit depend on it.  Set VS_CELLCHAT_RUN_PATTERNS=true only after
# installing a compatible NMF stack and independently validating the output.
pattern_status <- data.frame(
  requested = FALSE,
  status = "not_run",
  note = "k=4 NMF pattern analysis is exploratory and is not required for the reported CellChat network/LR results",
  stringsAsFactors = FALSE
)
write.csv(pattern_status, file.path(OUT, "cellchat_pattern_status.csv"), row.names = FALSE)
saveRDS(cellchat, file.path(OUT, "cellchat_rerun_numeric_with_patterns.rds"))

compare_table <- data.frame(
  metric = character(), rerun = numeric(), cached = numeric(), stringsAsFactors = FALSE
)
add_metric <- function(metric, rerun, cached = NA_real_) {
  compare_table <<- rbind(compare_table,
                           data.frame(metric = metric, rerun = rerun, cached = cached))
}

if (file.exists(OLD_ALL)) {
  old_all <- read.csv(OLD_ALL, check.names = FALSE)
  key_cols <- intersect(c("source", "target", "ligand", "receptor", "pathway_name"),
                        intersect(colnames(old_all), colnames(lr_all)))
  old_key <- do.call(paste, c(old_all[key_cols], sep = "|"))
  new_key <- do.call(paste, c(lr_all[key_cols], sep = "|"))
  common <- intersect(old_key, new_key)
  add_metric("all_LR_rows", nrow(lr_all), nrow(old_all))
  add_metric("all_LR_key_overlap", length(common), length(old_key))
  add_metric("all_LR_key_jaccard", length(common) / length(union(old_key, new_key)), NA)
  if (length(common)) {
    old_i <- match(common, old_key)
    new_i <- match(common, new_key)
    add_metric("all_LR_prob_max_abs_diff",
               max(abs(lr_all$prob[new_i] - old_all$prob[old_i]), na.rm = TRUE), NA)
    add_metric("all_LR_pval_max_abs_diff",
               max(abs(lr_all$pval[new_i] - old_all$pval[old_i]), na.rm = TRUE), NA)
  }
}

if (file.exists(OLD_FOCUS)) {
  old_focus <- read.csv(OLD_FOCUS, check.names = FALSE)
  focus_key_cols <- intersect(c("source", "target", "ligand", "receptor", "pathway_name", "direction"),
                              intersect(colnames(old_focus), colnames(lr_focus)))
  old_key <- do.call(paste, c(old_focus[focus_key_cols], sep = "|"))
  new_key <- do.call(paste, c(lr_focus[focus_key_cols], sep = "|"))
  common <- intersect(old_key, new_key)
  add_metric("focus_LR_rows", nrow(lr_focus), nrow(old_focus))
  add_metric("focus_LR_key_overlap", length(common), length(old_key))
  add_metric("focus_LR_key_jaccard", length(common) / length(union(old_key, new_key)), NA)
  if (length(common)) {
    old_i <- match(common, old_key)
    new_i <- match(common, new_key)
    add_metric("focus_LR_prob_max_abs_diff",
               max(abs(lr_focus$prob[new_i] - old_focus$prob[old_i]), na.rm = TRUE), NA)
    add_metric("focus_LR_pval_max_abs_diff",
               max(abs(lr_focus$pval[new_i] - old_focus$pval[old_i]), na.rm = TRUE), NA)
  }
}

write.csv(compare_table, file.path(OUT, "cellchat_cache_comparison.csv"), row.names = FALSE)
write.csv(head(lr_focus[order(-lr_focus$prob), ], 50),
          file.path(OUT, "LR_C3high_Tumor_top50_rerun.csv"), row.names = FALSE)
message("CellChat numeric rerun finished: ", nrow(lr_all), " all-LR rows; ",
        nrow(lr_focus), " focused rows.")
