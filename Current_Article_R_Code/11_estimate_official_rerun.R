# =============================================================================
# Analysis 13: Official ESTIMATE score rerun
#
# This recalculates official ESTIMATE scores from the supplied Step1 cleaned
# TPM matrix and provides the score table used by the primary MGS composition
# analysis (08_mgs_variance_partitioning_estimate.R).
# =============================================================================

options(stringsAsFactors = FALSE)
source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
suppressPackageStartupMessages({
  library(readxl)
  library(estimate)
})

out_dir <- file.path(AUDIT_ROOT, "estimate_official")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
required <- c(STEP3_RDATA, CLINICAL_XLSX)
if (!file.exists(SERVER_RDATA_DIR)) stop("Missing server RData directory: ", SERVER_RDATA_DIR)
step1_path <- file.path(SERVER_RDATA_DIR, "Step1_Clean_Data_VS3.RData")
if (!file.exists(step1_path)) stop("Missing Step1 cleaned TPM RData: ", step1_path)
if (!file.exists(CLINICAL_XLSX)) stop("Missing clinical workbook: ", CLINICAL_XLSX)

step1 <- new.env(parent = emptyenv())
load(step1_path, envir = step1)
if (!is.matrix(step1$clean_tpm)) stop("Step1 RData does not contain matrix clean_tpm")
tpm <- step1$clean_tpm
if (ncol(tpm) != 38L) stop("Expected 38 Step1 samples")
if (anyDuplicated(rownames(tpm))) stop("Step1 gene symbols are not unique")
if (any(!is.finite(tpm))) stop("Non-finite TPM value in Step1 matrix")

tmp_dir <- file.path(out_dir, "estimate_tmp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
input_file <- file.path(tmp_dir, "tpm_for_estimate.txt")
filtered_file <- file.path(tmp_dir, "tpm_filtered_estimate.gct")
scores_file <- file.path(tmp_dir, "estimate_scores.gct")
input_df <- cbind(NAME = rownames(tpm), Description = rownames(tpm),
                  as.data.frame(tpm, check.names = FALSE))
write.table(input_df, input_file, sep = "\t", quote = FALSE, row.names = FALSE)
estimate::filterCommonGenes(input.f = input_file, output.f = filtered_file,
                            id = "GeneSymbol")
estimate::estimateScore(input.ds = filtered_file, output.ds = scores_file,
                        platform = "illumina")

raw <- read.table(scores_file, header = TRUE, skip = 2, sep = "\t",
                  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
score_rows <- raw[match(c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
                        rownames(raw)), , drop = FALSE]
if (nrow(score_rows) != 3L) stop("Official ESTIMATE score rows not found")
desc_cols <- grep("^Description", colnames(score_rows), value = TRUE)
score_mat <- as.matrix(score_rows[, !colnames(score_rows) %in% desc_cols, drop = FALSE])
storage.mode(score_mat) <- "numeric"
scores <- as.data.frame(t(score_mat), stringsAsFactors = FALSE)
colnames(scores) <- c("Stromal_Score", "Immune_Score", "ESTIMATE_Score")
scores$SampleID <- rownames(scores)
scores$Tumour_Purity <- pmin(pmax(
  cos(0.6049872018 + 0.0001467884 * scores$ESTIMATE_Score), 0), 1)
rownames(scores) <- NULL

cl <- read_excel(CLINICAL_XLSX)
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", , drop = FALSE]
cl$SampleID <- as.character(cl$SampleID)
scores$Subtype <- as.character(cl$Subtype[match(scores$SampleID, cl$SampleID)])
if (anyNA(scores$Subtype)) stop("Official ESTIMATE scores cannot be aligned to clinical subtypes")

summary_df <- aggregate(cbind(Stromal_Score, Immune_Score, ESTIMATE_Score,
                              Tumour_Purity) ~ Subtype, data = scores, FUN = mean)
p_df <- data.frame(
  metric = c("Stromal_Score_KW_P", "Immune_Score_KW_P", "ESTIMATE_Score_KW_P", "Tumour_Purity_KW_P"),
  P_value = c(
    kruskal.test(Stromal_Score ~ Subtype, scores)$p.value,
    kruskal.test(Immune_Score ~ Subtype, scores)$p.value,
    kruskal.test(ESTIMATE_Score ~ Subtype, scores)$p.value,
    kruskal.test(Tumour_Purity ~ Subtype, scores)$p.value
  )
)
write.csv(scores, file.path(out_dir, "estimate_official_scores.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "estimate_official_subtype_means.csv"), row.names = FALSE)
write.csv(p_df, file.path(out_dir, "estimate_official_subtype_kruskal.csv"), row.names = FALSE)
writeLines(capture.output({
  cat("Official ESTIMATE rerun from Step1 cleaned TPM\n")
  cat("Step1 genes:", nrow(tpm), "samples:", ncol(tpm), "\n")
  cat("Official score rows produced:", nrow(score_rows), "\n\n")
  print(summary_df, row.names = FALSE)
  print(p_df, row.names = FALSE)
}), file.path(out_dir, "estimate_official_summary.txt"))
saveRDS(list(scores = scores, subtype_means = summary_df, kruskal = p_df,
             provenance = list(step1_rdata = step1_path, clinical = CLINICAL_XLSX,
                               estimate_package = as.character(packageVersion("estimate")))),
        file.path(out_dir, "estimate_official_rerun.rds"))
cat(paste(readLines(file.path(out_dir, "estimate_official_summary.txt")), collapse = "\n"), "\n")
