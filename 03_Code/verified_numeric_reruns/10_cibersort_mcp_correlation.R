# =============================================================================
# Analysis 12: Formal CIBERSORT--MCP-counter cross-validation
#
# MCP-counter is recalculated from the supplied 38-sample TPM matrix and
# compared with the supplied CIBERSORT result table. Supplementary Table 3 is only a
# rounded reporting target for an audit comparison.
# =============================================================================

options(stringsAsFactors = FALSE)
source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
suppressPackageStartupMessages({
  library(readxl)
  library(MCPcounter)
})

out_dir <- file.path(AUDIT_ROOT, "cibersort_mcp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
required <- c(BULK_TPM_XLSX, CLINICAL_XLSX, CIBERSORT_CSV)
if (any(!file.exists(required))) {
  stop("Missing CIBERSORT--MCP input(s): ",
       paste(required[!file.exists(required)], collapse = "; "))
}

# ---- Exact 38-sample discovery cohort ----
tpm_raw <- read_excel(BULK_TPM_XLSX)
genes <- as.character(tpm_raw$gene_name)
keep_gene <- !is.na(genes) & nzchar(genes) & !duplicated(genes)
tpm_raw <- tpm_raw[keep_gene, , drop = FALSE]
genes <- genes[keep_gene]
cl <- read_excel(CLINICAL_XLSX)
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", , drop = FALSE]
clinical_ids <- unique(as.character(cl$SampleID))
available_tpm <- grep("^tpm\\.", names(tpm_raw), value = TRUE)
available_ids <- sub("^tpm\\.", "", available_tpm)
our_samples <- available_tpm[available_ids %in% clinical_ids]
if (length(our_samples) != 38L) stop("Expected 38 clinical TPM samples")
tpm_mat <- as.matrix(tpm_raw[, our_samples, drop = FALSE])
storage.mode(tpm_mat) <- "numeric"
rownames(tpm_mat) <- genes
sample_ids <- sub("^tpm\\.", "", our_samples)

# ---- Fresh MCP-counter scores ----
mcp <- MCPcounter::MCPcounter.estimate(tpm_mat, featuresType = "HUGO_symbols")
mcp_rows <- c("CD8 T cells", "B lineage", "NK cells", "Monocytic lineage",
              "T cells", "Neutrophils", "Myeloid dendritic cells")
if (!all(mcp_rows %in% rownames(mcp))) stop("Required MCP-counter row missing")

# ---- CIBERSORT table aligned by sample ID ----
cib <- read.csv(CIBERSORT_CSV, check.names = FALSE, row.names = 1)
cib <- as.data.frame(cib, check.names = FALSE)
cib_ids <- rownames(cib)
if (!all(our_samples %in% cib_ids)) {
  stop("CIBERSORT table does not contain all clinical TPM sample IDs")
}
cib <- cib[our_samples, , drop = FALSE]

submission <- NULL
if (file.exists(SUBMISSION_TABLE3_CSV)) {
  submission <- read.csv(SUBMISSION_TABLE3_CSV, check.names = FALSE, stringsAsFactors = FALSE)
  required_table_cols <- c("Cell_Population", "MCP_CellType", "CIBERSORT_CellType",
                           "N", "Spearman_Rho", "P_value", "Mean_MCP", "Mean_CIBERSORT")
  if (!all(required_table_cols %in% names(submission))) {
    stop("Supplementary Table 3 is missing required columns")
  }
}

cib_vector <- function(label) {
  if (label == "B cells naive+memory") {
    return(as.numeric(cib[["B cells naive"]] + cib[["B cells memory"]]))
  }
  if (label == "NK resting+activated") {
    return(as.numeric(cib[["NK cells resting"]] + cib[["NK cells activated"]]))
  }
  if (label == "Macrophages M0+M1+M2") {
    return(as.numeric(cib[["Macrophages M0"]] + cib[["Macrophages M1"]] + cib[["Macrophages M2"]]))
  }
  if (!label %in% names(cib)) stop("Missing CIBERSORT column: ", label)
  as.numeric(cib[[label]])
}

mcp_vector <- function(label) {
  if (!label %in% rownames(mcp)) stop("Missing MCP-counter row: ", label)
  as.numeric(mcp[label, our_samples])
}

rows <- list(
  data.frame(Cell_Population = "CD8 T cells", MCP_CellType = "CD8 T cells",
             CIBERSORT_CellType = "T cells CD8"),
  data.frame(Cell_Population = "B cells naive", MCP_CellType = "B lineage",
             CIBERSORT_CellType = "B cells naive"),
  data.frame(Cell_Population = "NK cells resting", MCP_CellType = "NK cells",
             CIBERSORT_CellType = "NK cells resting"),
  data.frame(Cell_Population = "Macrophages M2", MCP_CellType = "Monocytic lineage",
             CIBERSORT_CellType = "Macrophages M2"),
  data.frame(Cell_Population = "Monocytes", MCP_CellType = "Monocytic lineage",
             CIBERSORT_CellType = "Monocytes"),
  data.frame(Cell_Population = "T cells CD4 memory resting", MCP_CellType = "T cells",
             CIBERSORT_CellType = "T cells CD4 memory resting"),
  data.frame(Cell_Population = "Neutrophils", MCP_CellType = "Neutrophils",
             CIBERSORT_CellType = "Neutrophils"),
  data.frame(Cell_Population = "Dendritic cells resting", MCP_CellType = "Myeloid dendritic cells",
             CIBERSORT_CellType = "Dendritic cells resting"),
  data.frame(Cell_Population = "NK cells activated", MCP_CellType = "NK cells",
             CIBERSORT_CellType = "NK cells activated"),
  data.frame(Cell_Population = "B cells memory", MCP_CellType = "B lineage",
             CIBERSORT_CellType = "B cells memory"),
  data.frame(Cell_Population = "B cells (naive+memory COMBINED)", MCP_CellType = "B lineage",
             CIBERSORT_CellType = "B cells naive+memory"),
  data.frame(Cell_Population = "NK cells (resting+activated COMBINED)", MCP_CellType = "NK cells",
             CIBERSORT_CellType = "NK resting+activated"),
  data.frame(Cell_Population = "Macrophages (M0+M1+M2 COMBINED)", MCP_CellType = "Monocytic lineage",
             CIBERSORT_CellType = "Macrophages M0+M1+M2")
)
rows <- do.call(rbind, rows)

result_rows <- lapply(seq_len(nrow(rows)), function(i) {
  r <- rows[i, ]
  x <- mcp_vector(r$MCP_CellType)
  y <- cib_vector(r$CIBERSORT_CellType)
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    Cell_Population = r$Cell_Population,
    MCP_CellType = r$MCP_CellType,
    CIBERSORT_CellType = r$CIBERSORT_CellType,
    N = length(x),
    Spearman_Rho = unname(test$estimate),
    P_value = test$p.value,
    Mean_MCP = mean(x),
    Mean_CIBERSORT = mean(y),
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, result_rows)

# Compare to the rounded Supplementary Table 3 target, if supplied.
result$Rho_abs_diff <- NA_real_
result$P_abs_diff <- NA_real_
result$Mean_MCP_abs_diff <- NA_real_
result$Mean_CIBERSORT_abs_diff <- NA_real_
result$Matches_rounded_submission <- NA
if (!is.null(submission)) {
  idx <- match(result$Cell_Population, submission$Cell_Population)
  if (anyNA(idx)) stop("Calculated row missing from Supplementary Table 3")
  result$Rho_abs_diff <- abs(result$Spearman_Rho - submission$Spearman_Rho[idx])
  result$P_abs_diff <- abs(result$P_value - submission$P_value[idx])
  result$Mean_MCP_abs_diff <- abs(result$Mean_MCP - submission$Mean_MCP[idx])
  result$Mean_CIBERSORT_abs_diff <- abs(result$Mean_CIBERSORT - submission$Mean_CIBERSORT[idx])
  p_tol <- pmax(5e-5, 0.05 * pmax(submission$P_value[idx], 1e-300))
  result$Matches_rounded_submission <- result$N == submission$N[idx] &
    result$Rho_abs_diff <= 0.001 &
    result$P_abs_diff <= p_tol &
    result$Mean_MCP_abs_diff <= 0.01 &
    result$Mean_CIBERSORT_abs_diff <= 0.001
}
overall_status <- if (all(result$Matches_rounded_submission %in% TRUE)) "PASS" else "WARN"

write.csv(result, file.path(out_dir, "cibersort_mcp_correlation.csv"), row.names = FALSE)
report <- c(
  "Discover Oncology VS CIBERSORT--MCP-counter formal rerun",
  paste("R:", R.version.string),
  paste("MCPcounter:", as.character(packageVersion("MCPcounter"))),
  paste("Samples:", length(our_samples)),
  paste("Submission comparison:", overall_status),
  paste(capture.output(print(result, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_dir, "cibersort_mcp_correlation_summary.txt"))
saveRDS(list(
  result = result,
  status = overall_status,
  input_tpm = normalizePath(BULK_TPM_XLSX, mustWork = TRUE),
  input_clinical = normalizePath(CLINICAL_XLSX, mustWork = TRUE),
  input_cibersort = normalizePath(CIBERSORT_CSV, mustWork = TRUE),
  samples = sample_ids
), file.path(out_dir, "cibersort_mcp_correlation.rds"))
cat(paste(report, collapse = "\n"), "\n")
