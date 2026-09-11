# =============================================================================
# Analysis 10 (primary revision): MGS variance partitioning with official
# ESTIMATE composition scores.
#
# This script replaces the historical signature-proxy decomposition as the
# manuscript-facing analysis. Official ESTIMATE Immune and Stromal scores
# are entered sequentially. ESTIMATE Score and tumour purity are not entered
# into the same variance model because ESTIMATE Score is their sum and tumour
# purity is a deterministic transformation of ESTIMATE Score.
#
# The historical proxy implementation remains in
# 08_mgs_variance_partitioning.R for transparent sensitivity/provenance only;
# it is not called by the primary runner and its 71.8% value is not used by
# the revised manuscript.
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
work_dir <- Sys.getenv("VS_MGS_OUTPUT_DIR", unset = file.path(AUDIT_ROOT, "mgs"))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(BULK_TPM_XLSX, CLINICAL_XLSX, ESTIMATE_SCORES_CSV)
if (any(!file.exists(required_inputs))) {
  stop("Missing required input file(s): ",
       paste(required_inputs[!file.exists(required_inputs)], collapse = "; "))
}

mgs_path <- Sys.getenv("VS_MGS_RERUN_RDS", unset = file.path(work_dir, "mgs_unbiased_results.rds"))
if (!file.exists(mgs_path)) {
  stop("Run 07_mgs_unbiased_immune_anchoring.R first; missing: ", mgs_path)
}

# ---- Load the exact 38-sample discovery cohort ----
tpm_raw <- readxl::read_excel(BULK_TPM_XLSX)
genes <- as.character(tpm_raw$gene_name)
keep_gene <- !is.na(genes) & nzchar(genes) & !duplicated(genes)
tpm_raw <- tpm_raw[keep_gene, , drop = FALSE]
genes <- genes[keep_gene]

cl <- readxl::read_excel(CLINICAL_XLSX)
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", , drop = FALSE]
clinical_ids <- unique(as.character(cl$SampleID))
available_tpm <- grep("^tpm\\.", names(tpm_raw), value = TRUE)
available_ids <- sub("^tpm\\.", "", available_tpm)
our_samples <- available_tpm[available_ids %in% clinical_ids]
if (length(our_samples) != 38L) {
  stop("Expected 38 discovery samples, found ", length(our_samples))
}

tpm_mat <- as.matrix(tpm_raw[, our_samples, drop = FALSE])
storage.mode(tpm_mat) <- "numeric"
rownames(tpm_mat) <- genes
log_tpm <- log2(tpm_mat + 1)
sample_ids_clean <- sub("^tpm\\.", "", our_samples)
subtypes <- setNames(as.character(cl$Subtype), as.character(cl$SampleID))[sample_ids_clean]
subtype_counts <- table(subtypes)
expected_counts <- c(C1 = 13L, C2 = 12L, C3 = 13L)
if (anyNA(subtypes) || !all(unname(subtype_counts[names(expected_counts)]) ==
                             unname(expected_counts))) {
  stop("The selected clinical cohort is not the expected C1/C2/C3 = 13/12/13 cohort")
}

# ---- Formally rerun MGS input ----
mgs_obj <- readRDS(mgs_path)
if (!is.data.frame(mgs_obj$mgs_unbiased) ||
    !all(c("Sample", "MGS_unbiased") %in% names(mgs_obj$mgs_unbiased))) {
  stop("MGS RDS does not contain mgs_unbiased with Sample and MGS_unbiased columns")
}
mgs_tab <- mgs_obj$mgs_unbiased
mgs_keys <- sub("^tpm\\.", "", as.character(mgs_tab$Sample))
mgs_values <- setNames(as.numeric(mgs_tab$MGS_unbiased), mgs_keys)[sample_ids_clean]
if (anyNA(mgs_values)) stop("MGS results do not cover all 38 clinical samples")

# ---- Load official ESTIMATE output and enforce sample-level alignment ----
estimate_scores <- read.csv(ESTIMATE_SCORES_CSV, check.names = FALSE,
                            stringsAsFactors = FALSE)
required_scores <- c("SampleID", "Stromal_Score", "Immune_Score",
                     "ESTIMATE_Score", "Tumour_Purity")
if (!all(required_scores %in% names(estimate_scores))) {
  stop("Official ESTIMATE CSV is missing: ",
       paste(setdiff(required_scores, names(estimate_scores)), collapse = "; "))
}
estimate_scores$SampleID <- sub("^tpm\\.", "", as.character(estimate_scores$SampleID))
estimate_scores <- estimate_scores[match(sample_ids_clean, estimate_scores$SampleID), , drop = FALSE]
if (anyNA(estimate_scores$SampleID) || !identical(estimate_scores$SampleID, sample_ids_clean)) {
  stop("Official ESTIMATE scores cannot be aligned one-to-one to the 38 samples")
}
if (any(!vapply(estimate_scores[required_scores[-1]], function(x) {
  all(is.finite(as.numeric(x)))
}, logical(1)))) {
  stop("Non-finite official ESTIMATE score detected")
}

# ---- Schwann-cell score used in the manuscript's composition check ----
sc_genes <- c("SOX10", "PLP1", "EGR2", "MPZ", "MBP", "PMP22", "MAL", "PRX",
              "CDH19", "L1CAM", "GAP43", "NGFR", "GFAP", "S100B")
sc_avail <- intersect(sc_genes, rownames(log_tpm))
if (length(sc_avail) < 5L) stop("Too few Schwann-cell score genes available")
sc_score <- colMeans(log_tpm[sc_avail, , drop = FALSE])

df <- data.frame(
  Sample = sample_ids_clean,
  Subtype = unname(subtypes),
  MGS = unname(mgs_values),
  SC_Score = unname(sc_score[our_samples]),
  Tumour_Purity = as.numeric(estimate_scores$Tumour_Purity),
  Immune_Score = as.numeric(estimate_scores$Immune_Score),
  Stromal_Score = as.numeric(estimate_scores$Stromal_Score),
  ESTIMATE_Score = as.numeric(estimate_scores$ESTIMATE_Score),
  stringsAsFactors = FALSE
)
if (any(!vapply(df, function(x) {
  if (is.numeric(x)) all(is.finite(x)) else all(!is.na(x))
}, logical(1)))) {
  stop("Non-finite values detected in the official ESTIMATE analysis table")
}

# ---- Hierarchical variance partitioning ----
# Only the two primary ESTIMATE component scores enter this model. The order
# is fixed here and is also reported in the output; sequential Delta-R2 values
# are order-dependent, whereas the final R2 is the joint explanatory value.
m0 <- lm(MGS ~ 1, data = df)
m1 <- lm(MGS ~ Stromal_Score, data = df)
m2 <- lm(MGS ~ Stromal_Score + Immune_Score, data = df)
models <- list(
  "Null" = m0,
  "Stromal score only" = m1,
  "Official ESTIMATE stromal + immune scores" = m2
)

step_rows <- list()
prev_r2 <- 0
for (mname in names(models)) {
  r2 <- summary(models[[mname]])$r.squared
  delta <- r2 - prev_r2
  step_rows[[length(step_rows) + 1L]] <- data.frame(
    Model = mname, R2 = r2, Delta_R2 = delta,
    stringsAsFactors = FALSE
  )
  prev_r2 <- r2
}
step_df <- bind_rows(step_rows)
stromal_r2 <- summary(m1)$r.squared
immune_r2 <- summary(m2)$r.squared - stromal_r2
total_r2 <- summary(m2)$r.squared
residual_r2 <- 1 - total_r2
decomp <- c(stromal = stromal_r2, immune = immune_r2, residual = residual_r2)

# ---- Purity-adjusted Schwann-cell score ----
sc_purity_model <- lm(SC_Score ~ Tumour_Purity, data = df)
df$SC_purity_adjusted <- residuals(sc_purity_model)
df$MGS_purity_adjusted <- residuals(lm(MGS ~ Tumour_Purity, data = df))
sc_mgs_model <- lm(SC_Score ~ MGS + Tumour_Purity, data = df)
# The partial correlation residualizes both variables on the same covariate.
sc_partial_r <- unname(cor(df$SC_purity_adjusted, df$MGS_purity_adjusted,
                           method = "pearson"))
sc_partial_p <- summary(sc_mgs_model)$coefficients["MGS", "Pr(>|t|)"]
kw_raw <- kruskal.test(SC_Score ~ Subtype, data = df)
kw_adj <- kruskal.test(SC_purity_adjusted ~ Subtype, data = df)
sc_raw_mgs_test <- suppressWarnings(cor.test(df$SC_Score, df$MGS,
                                             method = "spearman", exact = FALSE))

adj_summary <- df %>%
  group_by(Subtype) %>%
  summarise(
    SC_raw = mean(SC_Score),
    SC_purity_adj = mean(SC_purity_adjusted),
    Tumour_Purity_mean = mean(Tumour_Purity),
    Immune_Score_mean = mean(Immune_Score),
    Stromal_Score_mean = mean(Stromal_Score),
    .groups = "drop"
  )

# ---- Correlations and auditable outputs ----
cor_vars <- c("MGS", "SC_Score", "Tumour_Purity", "Immune_Score", "Stromal_Score")
cor_m <- cor(df[, cor_vars, drop = FALSE], method = "spearman")
mgs_purity_test <- suppressWarnings(cor.test(df$MGS, df$Tumour_Purity,
                                              method = "spearman", exact = FALSE))

write.csv(df, file.path(work_dir, "variance_partitioning_input_official_estimate.csv"),
          row.names = FALSE)
write.csv(step_df, file.path(work_dir, "variance_partitioning_steps_official_estimate.csv"),
          row.names = FALSE)
write.csv(adj_summary, file.path(work_dir, "sc_score_purity_adjusted_official_estimate.csv"),
          row.names = FALSE)
write.csv(cor_m, file.path(work_dir, "variance_partitioning_spearman_official_estimate.csv"))
write.csv(data.frame(
  metric = c("MGS_vs_Tumour_Purity_Spearman_rho", "MGS_vs_Tumour_Purity_Spearman_P",
             "SC_vs_MGS_partial_Pearson_r_controlling_Tumour_Purity",
             "SC_vs_MGS_partial_P_controlling_Tumour_Purity"),
  value = c(unname(mgs_purity_test$estimate), mgs_purity_test$p.value,
            sc_partial_r, sc_partial_p)
), file.path(work_dir, "official_estimate_association_statistics.csv"), row.names = FALSE)

summary_lines <- c(
  "Official ESTIMATE MGS variance partitioning",
  paste0("n_samples=", nrow(df)),
  paste0("estimate_scores_csv=", normalizePath(ESTIMATE_SCORES_CSV, mustWork = TRUE)),
  paste0("stromal_component_percent=", sprintf("%.10f", 100 * stromal_r2)),
  paste0("immune_component_percent=", sprintf("%.10f", 100 * immune_r2)),
  paste0("total_estimate_components_percent=", sprintf("%.10f", 100 * total_r2)),
  paste0("residual_percent=", sprintf("%.10f", 100 * residual_r2)),
  paste0("MGS_vs_Tumour_Purity_rho=", sprintf("%.10f", unname(mgs_purity_test$estimate))),
  paste0("MGS_vs_Tumour_Purity_p=", format.pval(mgs_purity_test$p.value, digits = 10)),
  paste0("SC_vs_MGS_partial_r_controlling_Tumour_Purity=", sprintf("%.10f", sc_partial_r)),
  paste0("SC_vs_MGS_partial_p_controlling_Tumour_Purity=", format.pval(sc_partial_p, digits = 10)),
  paste0("SC_vs_MGS_raw_Spearman_rho=", sprintf("%.10f", unname(sc_raw_mgs_test$estimate))),
  paste0("SC_vs_MGS_raw_Spearman_p=", format.pval(sc_raw_mgs_test$p.value, digits = 10)),
  paste0("SC_raw_subtype_KW_p=", format.pval(kw_raw$p.value, digits = 10)),
  paste0("SC_purity_adjusted_subtype_KW_p=", format.pval(kw_adj$p.value, digits = 10)),
  "Note: ESTIMATE Score and tumour purity were not entered with Immune/ Stromal scores because they are derived quantities."
)
writeLines(summary_lines, file.path(work_dir, "official_estimate_variance_summary.txt"))

result <- list(
  provenance = list(
    input_tpm = normalizePath(BULK_TPM_XLSX, mustWork = TRUE),
    input_clinical = normalizePath(CLINICAL_XLSX, mustWork = TRUE),
    estimate_scores = normalizePath(ESTIMATE_SCORES_CSV, mustWork = TRUE),
    mgs_input = normalizePath(mgs_path, mustWork = TRUE),
    estimate_definition = "Official estimate 1.0.13 output; Tumour_Purity = cos(0.6049872018 + 0.0001467884 * ESTIMATE_Score)",
    variance_model_definition = "MGS ~ Stromal_Score, then MGS ~ Stromal_Score + Immune_Score; ESTIMATE_Score and Tumour_Purity excluded from the same model as derived quantities",
    n_samples = nrow(df), sample_ids = sample_ids_clean
  ),
  variance_partitioning = list(
    mgs_r2_decomposition = decomp,
    total_estimate_components_r2 = total_r2,
    sc_vs_mgs_partial_r_controlling_tumour_purity = sc_partial_r,
    sc_vs_mgs_partial_p_controlling_tumour_purity = sc_partial_p,
    sc_vs_mgs_raw_spearman_rho = unname(sc_raw_mgs_test$estimate),
    sc_vs_mgs_raw_spearman_p = sc_raw_mgs_test$p.value,
    sc_kw_raw = kw_raw$p.value,
    sc_kw_adj = kw_adj$p.value,
    mgs_vs_tumour_purity_rho = unname(mgs_purity_test$estimate),
    mgs_vs_tumour_purity_p = mgs_purity_test$p.value
  ),
  models = models,
  stepwise = step_df,
  adjusted_summary = adj_summary,
  spearman = cor_m,
  df = df
)
saveRDS(result, file.path(work_dir, "variance_partitioning_results_official_estimate.rds"))

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("========== Official ESTIMATE variance analysis complete ==========\n")
