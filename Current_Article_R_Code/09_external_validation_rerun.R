# =============================================================================
# Analysis 11: Independent external validation rerun
#
# Recomputes ssGSEA subtype scores and ambiguous calls directly from the two
# supplied external expression matrices. Saved validation results are used only
# for an optional audit comparison, never as the classification input.
# =============================================================================

options(stringsAsFactors = FALSE)
source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
if (nzchar(GSVA_COMPAT_LIBRARY) && dir.exists(GSVA_COMPAT_LIBRARY)) {
  .libPaths(unique(c(GSVA_COMPAT_LIBRARY, .libPaths())))
}
suppressPackageStartupMessages({
  library(GSVA)
  library(limma)
})

out_dir <- file.path(AUDIT_ROOT, "external_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
required <- c(STEP3_RDATA, EXTERNAL_GSE141801_EXPR, EXTERNAL_GSE141801_CLIN,
              EXTERNAL_GSE39645_EXPR, EXTERNAL_GSE39645_CLIN)
if (any(!file.exists(required))) {
  stop("Missing external-validation input(s): ",
       paste(required[!file.exists(required)], collapse = "; "))
}

step3 <- new.env(parent = emptyenv())
load(STEP3_RDATA, envir = step3)
if (!exists("gene_sets", envir = step3)) stop("STEP3_RDATA does not contain gene_sets")
gene_sets <- step3$gene_sets
if (!is.list(gene_sets) || length(gene_sets) != 3L) stop("Expected three external gene sets")
gene_set_lengths <- vapply(gene_sets, length, integer(1))

load_expression <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2L) stop("Expression file has fewer than two columns: ", path)
  gene_col <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE])
  rownames(mat) <- gene_col
  storage.mode(mat) <- "numeric"
  mat <- mat[!is.na(rownames(mat)) & nzchar(rownames(mat)), , drop = FALSE]
  if (anyDuplicated(rownames(mat))) mat <- limma::avereps(mat)
  mat
}

classify <- function(expr, sets) {
  if (packageVersion("GSVA") >= "1.50.0") {
    scores <- GSVA::gsva(GSVA::ssgseaParam(expr, sets), verbose = FALSE)
  } else {
    scores <- GSVA::gsva(expr, sets, method = "ssgsea", kcdf = "Gaussian")
  }
  scaled <- t(scale(t(scores)))
  predicted <- apply(scaled, 2, function(x) rownames(scaled)[which.max(x)])
  sample_rows <- lapply(colnames(scaled), function(s) {
    z <- sort(scaled[, s], decreasing = TRUE)
    data.frame(
      Sample = s,
      Predicted = names(z)[1],
      Top1 = unname(z[1]),
      Top2 = unname(z[2]),
      Delta = unname(z[1] - z[2]),
      Ambiguous_delta_lt_0.1 = unname((z[1] - z[2]) < 0.1),
      stringsAsFactors = FALSE
    )
  })
  list(raw = scores, scaled = scaled, predicted = predicted,
       sample_table = do.call(rbind, sample_rows))
}

read_clinical <- function(path, expr_samples) {
  clin <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  id_col <- intersect(c("ID", "Sample_geo_accession", "SampleID"), names(clin))[1]
  if (is.na(id_col)) stop("No sample ID column in clinical file: ", path)
  rownames(clin) <- as.character(clin[[id_col]])
  clin[match(expr_samples, rownames(clin)), , drop = FALSE]
}

datasets <- list(
  GSE141801 = list(expr = EXTERNAL_GSE141801_EXPR, clin = EXTERNAL_GSE141801_CLIN),
  GSE39645 = list(expr = EXTERNAL_GSE39645_EXPR, clin = EXTERNAL_GSE39645_CLIN)
)
expected_ambiguous <- c(GSE141801 = 4L, GSE39645 = 2L)
all_scores <- list()
all_summaries <- list()
all_markers <- list()
all_nf2 <- list()

for (dataset_id in names(datasets)) {
  expr <- load_expression(datasets[[dataset_id]]$expr)
  rerun <- classify(expr, gene_sets)
  clin <- read_clinical(datasets[[dataset_id]]$clin, colnames(expr))
  if (any(is.na(rownames(clin)))) stop("Clinical data do not cover all expression samples for ", dataset_id)
  sample_table <- rerun$sample_table
  sample_table$Dataset <- dataset_id
  sample_table <- sample_table[, c("Dataset", "Sample", "Predicted", "Top1", "Top2",
                                   "Delta", "Ambiguous_delta_lt_0.1")]
  sample_table$NF2 <- as.character(clin$NF2[match(sample_table$Sample, rownames(clin))])
  marker_rows <- lapply(c("MKI67", "CD8A", "LAG3", "HAVCR2"), function(g) {
    if (!g %in% rownames(expr)) return(data.frame(Dataset = dataset_id, Gene = g, P_value = NA_real_))
    p <- tryCatch(kruskal.test(as.numeric(expr[g, ]) ~ factor(rerun$predicted))$p.value,
                  error = function(e) NA_real_)
    data.frame(Dataset = dataset_id, Gene = g, P_value = p)
  })
  marker_df <- do.call(rbind, marker_rows)
  all_scores[[dataset_id]] <- list(sample_table = sample_table,
                                    marker = marker_df, scaled_scores = rerun$scaled)
  all_markers[[dataset_id]] <- marker_df
  all_summaries[[dataset_id]] <- data.frame(
    Dataset = dataset_id,
    Expression_genes = nrow(expr),
    Samples = ncol(expr),
    Gene_set_C1 = gene_set_lengths[["C1_Proliferative"]],
    Gene_set_C2 = gene_set_lengths[["C2_Mesenchymal"]],
    Gene_set_C3 = gene_set_lengths[["C3_Immune"]],
    Ambiguous_count = sum(sample_table$Ambiguous_delta_lt_0.1),
    Expected_ambiguous_count = expected_ambiguous[[dataset_id]],
    Ambiguous_count_matches_expected = sum(sample_table$Ambiguous_delta_lt_0.1) == expected_ambiguous[[dataset_id]],
    stringsAsFactors = FALSE
  )
  all_nf2[[dataset_id]] <- sample_table[, c("Dataset", "Sample", "Predicted", "NF2")]
}

scores_long <- do.call(rbind, lapply(all_scores, function(x) x$sample_table))
summary_df <- do.call(rbind, all_summaries)
markers_df <- do.call(rbind, all_markers)
nf2_df <- do.call(rbind, all_nf2)

nf2_df$Predicted <- factor(nf2_df$Predicted,
                           levels = c("C1_Proliferative", "C2_Mesenchymal", "C3_Immune"))
nf2_df$NF2 <- factor(nf2_df$NF2, levels = c("Y", "N"))
nf2_tab <- table(nf2_df$NF2, nf2_df$Predicted)
combined_fisher_p <- fisher.test(nf2_tab)$p.value
nf2_summary <- aggregate(NF2 ~ Predicted, data = nf2_df,
                         FUN = function(x) sum(x == "Y", na.rm = TRUE))
totals <- as.data.frame(table(nf2_df$Predicted), stringsAsFactors = FALSE)
names(totals) <- c("Predicted", "Total")
names(nf2_summary)[2] <- "NF2_mutated"
nf2_summary <- merge(totals, nf2_summary, by = "Predicted", all.x = TRUE, sort = FALSE)
nf2_summary$NF2_mutated[is.na(nf2_summary$NF2_mutated)] <- 0L

cache_path <- file.path(SERVER_RDATA_DIR, "Step2_Ext_Validation_Results.RData")
cache_status <- "not available"
cache_max_abs_diff <- NA_real_
cache_pred_agreement <- NA_real_
if (file.exists(cache_path)) {
  cache <- new.env(parent = emptyenv())
  load(cache_path, envir = cache)
  if (exists("all_scores", envir = cache) && exists("all_predicted", envir = cache)) {
    diffs <- numeric()
    agrees <- logical()
    for (dataset_id in names(all_scores)) {
      if (!is.null(cache$all_scores[[dataset_id]])) {
        fresh <- all_scores[[dataset_id]]$scaled_scores
        old <- cache$all_scores[[dataset_id]]
        common <- intersect(colnames(fresh), colnames(old))
        if (length(common)) diffs <- c(diffs, abs(fresh[, common, drop = FALSE] - old[, common, drop = FALSE]))
      }
      if (!is.null(cache$all_predicted[[dataset_id]])) {
        fresh_pred <- all_scores[[dataset_id]]$sample_table$Predicted
        common <- intersect(names(cache$all_predicted[[dataset_id]]), all_scores[[dataset_id]]$sample_table$Sample)
        if (length(common)) agrees <- c(agrees, fresh_pred[match(common, all_scores[[dataset_id]]$sample_table$Sample)] ==
                                            cache$all_predicted[[dataset_id]][common])
      }
    }
    cache_max_abs_diff <- if (length(diffs)) max(diffs, na.rm = TRUE) else NA_real_
    cache_pred_agreement <- if (length(agrees)) mean(agrees) else NA_real_
    cache_status <- if (is.finite(cache_max_abs_diff) && isTRUE(cache_pred_agreement == 1)) "PASS" else "WARN"
  }
}

write.csv(scores_long, file.path(out_dir, "external_validation_scores.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "external_validation_summary.csv"), row.names = FALSE)
write.csv(markers_df, file.path(out_dir, "external_validation_marker_tests.csv"), row.names = FALSE)
write.csv(nf2_summary, file.path(out_dir, "external_validation_NF2_by_subtype.csv"), row.names = FALSE)
write.csv(nf2_df, file.path(out_dir, "external_validation_NF2_sample_level.csv"), row.names = FALSE)

report <- c(
  "Discover Oncology VS independent external validation rerun",
  paste("R:", R.version.string),
  paste("GSVA:", as.character(packageVersion("GSVA"))),
  paste("Gene sets:", paste(names(gene_sets), gene_set_lengths, sep = "=", collapse = "; ")),
  paste("GSE141801 ambiguous:", summary_df$Ambiguous_count[summary_df$Dataset == "GSE141801"], "expected 4"),
  paste("GSE39645 ambiguous:", summary_df$Ambiguous_count[summary_df$Dataset == "GSE39645"], "expected 2"),
  paste("Combined NF2 Fisher P:", formatC(combined_fisher_p, format = "e", digits = 8)),
  paste("Cached score comparison:", cache_status,
        "max_abs_diff=", formatC(cache_max_abs_diff, digits = 8),
        "prediction_agreement=", formatC(cache_pred_agreement, digits = 8)),
  "",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"),
  paste(capture.output(print(nf2_summary, row.names = FALSE)), collapse = "\n"),
  paste(capture.output(print(markers_df, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_dir, "external_validation_summary.txt"))
saveRDS(list(scores = all_scores, summary = summary_df, markers = markers_df,
             nf2 = nf2_df, nf2_summary = nf2_summary,
             combined_fisher_p = combined_fisher_p,
             cache_status = cache_status, cache_max_abs_diff = cache_max_abs_diff,
             cache_pred_agreement = cache_pred_agreement,
             provenance = list(expression = sapply(datasets, function(x) x$expr),
                               clinical = sapply(datasets, function(x) x$clin),
                               step3_rdata = STEP3_RDATA)),
        file.path(out_dir, "external_validation_rerun.rds"))
cat(paste(report, collapse = "\n"), "\n")
