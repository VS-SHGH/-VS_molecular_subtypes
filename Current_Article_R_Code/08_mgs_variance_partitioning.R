# =============================================================================
# Historical sensitivity analysis: variance partitioning of MGS using proxies
#
# This script is retained only for provenance of the earlier 71.8% sensitivity
# result. The manuscript-facing analysis is now
# 08_mgs_variance_partitioning_estimate.R, which uses official ESTIMATE
# Immune/Stromal scores. The composition covariates called *_proxy below are
# explicitly defined signature means; they are not official ESTIMATE output.
# =============================================================================

library(dplyr)
library(readxl)
library(MCPcounter)

source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
source(file.path(code_dir, "composition_proxy_functions.R"))
work_dir <- Sys.getenv("VS_MGS_OUTPUT_DIR", unset = file.path(AUDIT_ROOT, "mgs"))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(BULK_TPM_XLSX, CLINICAL_XLSX)
if (any(!file.exists(required_inputs))) {
  stop("Missing required input file(s): ",
       paste(required_inputs[!file.exists(required_inputs)], collapse = "; "))
}
mgs_path <- Sys.getenv("VS_MGS_RERUN_RDS", unset = file.path(work_dir, "mgs_unbiased_results.rds"))
if (!file.exists(mgs_path)) {
  stop("Run 07_mgs_unbiased_immune_anchoring.R first; missing: ", mgs_path)
}

# ---- Load the exact 38-sample discovery cohort ----
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
if (length(our_samples) != 38L) {
  stop("Expected 38 discovery samples, found ", length(our_samples))
}

tpm_mat <- as.matrix(tpm_raw[, our_samples, drop = FALSE])
storage.mode(tpm_mat) <- "numeric"
rownames(tpm_mat) <- genes
log_tpm <- log2(tpm_mat + 1)
sample_ids_clean <- sub("^tpm\\.", "", our_samples)
subtypes <- setNames(as.character(cl$Subtype), as.character(cl$SampleID))[sample_ids_clean]
if (anyNA(subtypes) || !identical(as.integer(table(subtypes)), c(13L, 12L, 13L))) {
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

# ---- Explicit ESTIMATE-like composition proxies ----
# These reproduce the historical proxy definition used for the manuscript's
# 71.8% decomposition while making the definition and provenance auditable.
stromal_proxy <- signature_mean(tpm_mat, ESTIMATE_LIKE_STROMAL_GENES)
immune_proxy <- signature_mean(tpm_mat, ESTIMATE_LIKE_IMMUNE_GENES)
proxy_scale <- max(stromal_proxy + immune_proxy, na.rm = TRUE)
purity_proxy <- 1 - (stromal_proxy + immune_proxy) / proxy_scale

# ---- Fresh MCP-counter calculation from the same linear TPM input ----
mcp <- MCPcounter::MCPcounter.estimate(tpm_mat, featuresType = "HUGO_symbols")
required_mcp <- c("Fibroblasts", "Monocytic lineage", "T cells", "CD8 T cells",
                  "B lineage", "NK cells", "Neutrophils", "Endothelial cells",
                  "Myeloid dendritic cells", "Cytotoxic lymphocytes")
if (!all(required_mcp %in% rownames(mcp))) {
  stop("Fresh MCP-counter result is missing: ",
       paste(setdiff(required_mcp, rownames(mcp)), collapse = "; "))
}

# ---- Schwann-cell score used in the manuscript's composition check ----
sc_genes <- c("SOX10", "PLP1", "EGR2", "MPZ", "MBP", "PMP22", "MAL", "PRX",
              "CDH19", "L1CAM", "GAP43", "NGFR", "GFAP", "S100B")
sc_avail <- intersect(sc_genes, rownames(log_tpm))
if (length(sc_avail) < 5L) stop("Too few Schwann-cell score genes available")
sc_score <- colMeans(log_tpm[sc_avail, , drop = FALSE])

# ---- Analysis summary table ----
df <- data.frame(
  Sample = sample_ids_clean,
  Subtype = unname(subtypes),
  MGS = unname(mgs_values),
  SC_Score = unname(sc_score[our_samples]),
  Purity_proxy = unname(purity_proxy),
  Immune_proxy = unname(immune_proxy),
  Stromal_proxy = unname(stromal_proxy),
  Fibroblasts = as.numeric(mcp["Fibroblasts", our_samples]),
  Monocytes = as.numeric(mcp["Monocytic lineage", our_samples]),
  Tcells = as.numeric(mcp["T cells", our_samples]),
  CD8_Tcells = as.numeric(mcp["CD8 T cells", our_samples]),
  Bcells = as.numeric(mcp["B lineage", our_samples]),
  NK_cells_mcp = as.numeric(mcp["NK cells", our_samples]),
  Neutrophils = as.numeric(mcp["Neutrophils", our_samples]),
  Endothelial = as.numeric(mcp["Endothelial cells", our_samples]),
  mDC = as.numeric(mcp["Myeloid dendritic cells", our_samples]),
  CytoLymph = as.numeric(mcp["Cytotoxic lymphocytes", our_samples]),
  stringsAsFactors = FALSE
)
if (any(!vapply(df, function(x) {
  if (is.numeric(x)) all(is.finite(x)) else all(!is.na(x))
}, logical(1)))) {
  stop("Non-finite values detected in the variance-partitioning input table")
}

# ---- Hierarchical variance partitioning ----
cat("========== Hierarchical Variance Partitioning of MGS ==========\n\n")
cat("Composition covariates: explicit ESTIMATE-like proxies + fresh MCP-counter\n")
m0 <- lm(MGS ~ 1, data = df)
m1 <- lm(MGS ~ Purity_proxy, data = df)
m2 <- lm(MGS ~ Purity_proxy + Immune_proxy, data = df)
m3 <- lm(MGS ~ Purity_proxy + Immune_proxy + Fibroblasts, data = df)
m4 <- lm(MGS ~ Purity_proxy + Immune_proxy + Fibroblasts + Monocytes + Tcells + Bcells, data = df)
models <- list(
  "Null" = m0,
  "Purity proxy only" = m1,
  "Purity proxy + immune proxy" = m2,
  "Purity proxy + immune proxy + fibroblasts" = m3,
  "Full (proxies + fibro + 3 cell types)" = m4
)
step_rows <- list()
prev_r2 <- 0
for (mname in names(models)) {
  r2 <- summary(models[[mname]])$r.squared
  delta <- r2 - prev_r2
  step_rows[[length(step_rows) + 1L]] <- data.frame(Model = mname, R2 = r2, Delta_R2 = delta)
  cat(sprintf("%-46s R2=%.3f  Delta_R2=%.3f\n", mname, r2, delta))
  prev_r2 <- r2
}
step_df <- bind_rows(step_rows)
residual_r2 <- 1 - summary(m4)$r.squared

# ---- Schwann-cell score association and partial R2 ----
cat("\n========== SC Score Variance Check ==========\n\n")
s1 <- lm(SC_Score ~ MGS, data = df)
s2 <- lm(SC_Score ~ MGS + Purity_proxy, data = df)
s3 <- lm(SC_Score ~ MGS + Purity_proxy + Immune_proxy + Fibroblasts + Monocytes, data = df)
ss_mgs_only <- sum(residuals(lm(SC_Score ~ Purity_proxy + Immune_proxy + Fibroblasts + Monocytes, data = df))^2)
ss_mgs_full <- sum(residuals(s3)^2)
partial_r2_mgs <- 1 - ss_mgs_full / ss_mgs_only
cat(sprintf("SC_Score ~ MGS: R2=%.3f\n", summary(s1)$r.squared))
cat(sprintf("SC_Score ~ MGS + composition: R2=%.3f\n", summary(s3)$r.squared))
cat(sprintf("Partial R2 of MGS after composition adjustment: %.3f\n", partial_r2_mgs))

# ---- Decomposition summary ----
total_r2 <- summary(m4)$r.squared
comp_r2 <- summary(m1)$r.squared
immune_r2 <- summary(m2)$r.squared - comp_r2
fibro_r2 <- summary(m3)$r.squared - summary(m2)$r.squared
other_r2 <- summary(m4)$r.squared - summary(m3)$r.squared
decomp <- c(purity = comp_r2, immune = immune_r2, fibro = fibro_r2,
            other = other_r2, residual = residual_r2)
cat("\nMGS variance decomposition (%):\n")
print(round(100 * decomp, 1))

# ---- Purity-adjusted SC score ----
df$SC_purity_adjusted <- residuals(lm(SC_Score ~ Purity_proxy, data = df))
adj_summary <- df %>%
  group_by(Subtype) %>%
  summarise(
    SC_raw = mean(SC_Score),
    SC_purity_adj = mean(SC_purity_adjusted),
    Purity_proxy_mean = mean(Purity_proxy),
    .groups = "drop"
  )
kw_raw <- kruskal.test(SC_Score ~ Subtype, data = df)
kw_adj <- kruskal.test(SC_purity_adjusted ~ Subtype, data = df)

# ---- Correlation matrix and auditable intermediate files ----
vars <- c("MGS", "SC_Score", "Purity_proxy", "Immune_proxy", "Fibroblasts",
          "Monocytes", "Tcells", "Bcells")
cor_m <- cor(df[, vars, drop = FALSE], method = "spearman")
write.csv(df, file.path(work_dir, "variance_partitioning_input.csv"), row.names = FALSE)
write.csv(step_df, file.path(work_dir, "variance_partitioning_steps.csv"), row.names = FALSE)
write.csv(adj_summary, file.path(work_dir, "sc_score_purity_adjusted_by_subtype.csv"), row.names = FALSE)
write.csv(cor_m, file.path(work_dir, "variance_partitioning_spearman_correlations.csv"))

result <- list(
  provenance = list(
    input_tpm = normalizePath(BULK_TPM_XLSX, mustWork = TRUE),
    input_clinical = normalizePath(CLINICAL_XLSX, mustWork = TRUE),
    mgs_input = normalizePath(mgs_path, mustWork = TRUE),
    mgs_definition = "07_mgs_unbiased_immune_anchoring.R; top 3000 MAD genes; PCA PC1-5; principal_curve stretch=0; immune anchor",
    composition_definition = "Explicit ESTIMATE-like stromal/immune signature means; purity proxy = 1 - (stromal + immune)/maximum combined score; these are proxies, not official ESTIMATE output",
    mcp_definition = "MCPcounter::MCPcounter.estimate on linear TPM, featuresType=HUGO_symbols",
    n_samples = nrow(df),
    sample_ids = sample_ids_clean
  ),
  variance_partitioning = list(
    mgs_r2_decomposition = decomp,
    sc_partial_r2_mgs = partial_r2_mgs,
    sc_kw_raw = kw_raw$p.value,
    sc_kw_adj = kw_adj$p.value
  ),
  models = models,
  stepwise = step_df,
  adjusted_summary = adj_summary,
  spearman = cor_m,
  df = df
)
saveRDS(result, file.path(work_dir, "variance_partitioning_results.rds"))
cat(sprintf("\nTOTAL explained by composition: %.1f%%\n", 100 * total_r2))
cat("========== Complete ==========\n")
