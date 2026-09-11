# =============================================================================
# Figure 4e and 4h replacement panels for the Discover Oncology revision.
#
# Figure 4e: subject-level bootstrap stability of the prespecified MGS.
#   - 38 discovery samples; top 3,000 MAD genes.
#   - PCA on the bootstrap expression matrix, followed by principal_curve on
#     PC1-PC5, exactly matching the revised MGS construction.
#   - Direction is fixed only by the prespecified immune anchor; subtype labels
#     are never used in the bootstrap.
#   - For duplicate observations in a bootstrap sample, lambda values are
#     aggregated by the median for that subject. Kendall's tau is then
#     calculated between the full-data MGS ordering and the bootstrap ordering
#     among the unique resampled subjects.
#
# Figure 4h: official ESTIMATE component-based hierarchical variance
# partitioning. The historical signature-proxy decomposition is not used.
# The displayed components are the sequential R2 contributions of Stromal
# Score, Immune Score, and the residual from the joint model.
#
# This script writes auditable CSV/RDS/TXT outputs under the analysis output
# directory and separate replacement panels to the requested figure directory.
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(princurve)
  library(ggplot2)
  library(gridExtra)
})

source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))

out_dir <- Sys.getenv("VS_MGS_OUTPUT_DIR", unset = file.path(AUDIT_ROOT, "mgs"))
figure_dir <- Sys.getenv(
  "VS_FIGURE_OUTPUT_DIR",
  unset = file.path(AUDIT_ROOT, "figure4_updated_panels")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

TPM_FILE <- Sys.getenv("VS_TPM_XLSX", unset = BULK_TPM_XLSX)
CLINICAL_FILE <- Sys.getenv("VS_CLINICAL_XLSX", unset = CLINICAL_XLSX)
ESTIMATE_FILE <- Sys.getenv("VS_ESTIMATE_SCORES_CSV", unset = ESTIMATE_SCORES_CSV)
MGS_FILE <- Sys.getenv(
  "VS_MGS_RERUN_RDS",
  unset = file.path(out_dir, "mgs_unbiased_results.rds")
)

required_files <- c(TPM_FILE, CLINICAL_FILE, ESTIMATE_FILE, MGS_FILE)
if (any(!file.exists(required_files))) {
  stop("Missing required file(s): ",
       paste(required_files[!file.exists(required_files)], collapse = "; "))
}

# ----------------------------------------------------------------------------
# Load and align the exact 38-sample discovery cohort.
# ----------------------------------------------------------------------------
tpm_raw <- readxl::read_excel(TPM_FILE)
genes <- as.character(tpm_raw$gene_name)
keep_gene <- !is.na(genes) & nzchar(genes) & !duplicated(genes)
tpm_raw <- tpm_raw[keep_gene, , drop = FALSE]
genes <- genes[keep_gene]

cl <- readxl::read_excel(CLINICAL_FILE)
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
sample_ids <- sub("^tpm\\.", "", our_samples)
subtype_lookup <- setNames(as.character(cl$Subtype), as.character(cl$SampleID))
subtypes <- unname(subtype_lookup[sample_ids])
if (anyNA(subtypes)) stop("Subtype labels could not be aligned to all 38 samples")

# ----------------------------------------------------------------------------
# Recreate the fixed MGS feature matrix and biological immune anchor.
# ----------------------------------------------------------------------------
set.seed(20260905)
mads <- apply(log_tpm, 1, mad)
top_3000 <- names(sort(mads, decreasing = TRUE))[seq_len(min(3000L, length(mads)))]
expr_top <- log_tpm[top_3000, , drop = FALSE]

immune_genes <- c(
  "CD2", "CD3D", "CD3E", "CD8A", "GZMA", "PRF1", "CTLA4", "HAVCR2",
  "LAG3", "PDCD1", "CD68", "CD163", "CD74", "HLA-DRA", "CCL5",
  "CXCL9", "CXCL10"
)
immune_avail <- intersect(immune_genes, rownames(log_tpm))
if (length(immune_avail) < 5L) stop("Too few immune anchor genes are available")
immune_score <- colMeans(log_tpm[immune_avail, , drop = FALSE])

mgs_obj <- readRDS(MGS_FILE)
if (!is.data.frame(mgs_obj$mgs_unbiased) ||
    !all(c("Sample", "MGS_unbiased") %in% names(mgs_obj$mgs_unbiased))) {
  stop("MGS RDS does not contain the expected mgs_unbiased table")
}
mgs_tab <- mgs_obj$mgs_unbiased
mgs_keys <- sub("^tpm\\.", "", as.character(mgs_tab$Sample))
reference_mgs <- setNames(as.numeric(mgs_tab$MGS_unbiased), mgs_keys)[sample_ids]
if (anyNA(reference_mgs)) stop("Reference MGS does not cover all 38 samples")

# ----------------------------------------------------------------------------
# Figure 4e: canonical subject-level bootstrap.
# ----------------------------------------------------------------------------
N_BOOT <- 500L
bootstrap_one <- function(rep_id) {
  idx <- sample.int(length(sample_ids), replace = TRUE)
  sampled_ids <- sample_ids[idx]
  X <- t(expr_top[, idx, drop = FALSE])

  # A bootstrap can make a gene constant by chance. Such genes have no
  # estimable standardized PCA loading in that replicate and are removed only
  # for numerical stability; the full-data MGS still uses all top 3,000 genes.
  keep_variable <- apply(X, 2, function(z) is.finite(sd(z)) && sd(z) > 0)
  if (sum(keep_variable) < 100L) {
    return(data.frame(
      replicate = rep_id, kendall_tau = NA_real_, n_unique = length(unique(sampled_ids)),
      status = "failed_too_few_variable_genes", stringsAsFactors = FALSE
    ))
  }
  X <- X[, keep_variable, drop = FALSE]

  fit <- tryCatch({
    pb <- prcomp(X, center = TRUE, scale. = TRUE)
    n_pc <- min(5L, ncol(pb$x))
    if (n_pc < 2L) stop("fewer than two PCs")
    curve <- principal_curve(
      as.matrix(pb$x[, seq_len(n_pc), drop = FALSE]),
      smoother = "smooth_spline", trace = FALSE, stretch = 0
    )
    lam <- as.numeric(curve$lambda)
    anchor <- immune_score[idx]
    rho <- suppressWarnings(cor(lam, anchor, method = "spearman"))
    if (!is.finite(rho)) stop("non-finite immune-anchor correlation")
    if (rho < 0) lam <- max(lam) - lam

    # A resampled subject can occur more than once. Aggregate its fitted curve
    # positions so every subject has one value before Kendall's tau is computed.
    lam_by_subject <- tapply(lam, sampled_ids, median)
    common <- intersect(names(lam_by_subject), sample_ids)
    if (length(common) < 10L) stop("fewer than ten unique subjects")
    tau <- suppressWarnings(cor(
      rank(reference_mgs[common]),
      rank(as.numeric(lam_by_subject[common])),
      method = "kendall"
    ))
    if (!is.finite(tau)) stop("non-finite Kendall tau")
    data.frame(
      replicate = rep_id, kendall_tau = unname(tau),
      n_unique = length(unique(sampled_ids)), status = "ok",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      replicate = rep_id, kendall_tau = NA_real_,
      n_unique = length(unique(sampled_ids)),
      status = paste0("failed:", conditionMessage(e)), stringsAsFactors = FALSE
    )
  })
  fit
}

bootstrap_result_file <- file.path(out_dir, "figure4e_bootstrap_results.csv")
reuse_bootstrap <- identical(Sys.getenv("VS_REUSE_FIGURE4E", unset = "FALSE"), "TRUE") &&
  file.exists(bootstrap_result_file)
if (reuse_bootstrap) {
  bootstrap_df <- read.csv(bootstrap_result_file, stringsAsFactors = FALSE)
  N_BOOT <- max(bootstrap_df$replicate, na.rm = TRUE)
  cat("Reusing the audited Figure 4e bootstrap result table for redraw.\n")
} else {
  bootstrap_rows <- lapply(seq_len(N_BOOT), bootstrap_one)
  bootstrap_df <- bind_rows(bootstrap_rows)
}
valid_tau <- bootstrap_df$kendall_tau[is.finite(bootstrap_df$kendall_tau)]
if (length(valid_tau) < 450L) {
  stop("Fewer than 450 valid bootstrap replicates: ", length(valid_tau))
}
tau_mean <- mean(valid_tau)
tau_median <- median(valid_tau)
tau_ci <- as.numeric(stats::quantile(valid_tau, probs = c(0.025, 0.975), names = FALSE,
                                      type = 7))

write.csv(bootstrap_df, file.path(out_dir, "figure4e_bootstrap_results.csv"),
          row.names = FALSE)
write.csv(data.frame(
  metric = c("n_requested", "n_valid", "mean_kendall_tau", "median_kendall_tau",
             "percentile_2.5", "percentile_97.5", "median_n_unique_subjects"),
  value = c(N_BOOT, length(valid_tau), tau_mean, tau_median, tau_ci[1], tau_ci[2],
            median(bootstrap_df$n_unique[bootstrap_df$status == "ok"]))
), file.path(out_dir, "figure4e_bootstrap_summary.csv"), row.names = FALSE)
writeLines(c(
  "Figure 4e bootstrap stability",
  paste0("n_requested=", N_BOOT),
  paste0("n_valid=", length(valid_tau)),
  paste0("mean_kendall_tau=", sprintf("%.6f", tau_mean)),
  paste0("median_kendall_tau=", sprintf("%.6f", tau_median)),
  paste0("percentile_95_interval=[", sprintf("%.6f", tau_ci[1]), ", ",
         sprintf("%.6f", tau_ci[2]), "]"),
  paste0("median_n_unique_subjects=", median(bootstrap_df$n_unique[bootstrap_df$status == "ok"])),
  "Definition: each replicate resamples the 38 subjects with replacement, recomputes PCA and principal_curve on top 3,000 MAD genes using PC1-PC5, orients the curve only with the prespecified immune anchor, aggregates duplicate subjects by median lambda, and computes Kendall's tau against the full-data MGS ordering.",
  paste0("TPM_source=", normalizePath(TPM_FILE, mustWork = TRUE)),
  paste0("clinical_source=", normalizePath(CLINICAL_FILE, mustWork = TRUE)),
  paste0("mgs_reference=", normalizePath(MGS_FILE, mustWork = TRUE)),
  paste0("seed=20260905")
), file.path(out_dir, "figure4e_bootstrap_summary.txt"))
saveRDS(list(
  definition = "subject-level bootstrap; top 3000 MAD genes; PCA PC1-PC5; immune-anchor orientation; duplicate subject lambda median; Kendall tau against full-data MGS ordering",
  seed = 20260905L, n_boot = N_BOOT, results = bootstrap_df,
  summary = list(n_valid = length(valid_tau), mean = tau_mean, median = tau_median,
                 percentile_95 = tau_ci), sample_ids = sample_ids,
  reference_mgs = reference_mgs
), file.path(out_dir, "figure4e_bootstrap_results.rds"))

# ----------------------------------------------------------------------------
# Figure 4h: official ESTIMATE variance decomposition.
# ----------------------------------------------------------------------------
estimate_scores <- read.csv(ESTIMATE_FILE, check.names = FALSE,
                            stringsAsFactors = FALSE)
required_scores <- c("SampleID", "Stromal_Score", "Immune_Score",
                     "ESTIMATE_Score", "Tumour_Purity")
if (!all(required_scores %in% names(estimate_scores))) {
  stop("Official ESTIMATE CSV is missing: ",
       paste(setdiff(required_scores, names(estimate_scores)), collapse = "; "))
}
estimate_scores$SampleID <- sub("^tpm\\.", "", as.character(estimate_scores$SampleID))
estimate_scores <- estimate_scores[match(sample_ids, estimate_scores$SampleID), , drop = FALSE]
if (anyNA(estimate_scores$SampleID) || !identical(estimate_scores$SampleID, sample_ids)) {
  stop("Official ESTIMATE scores cannot be aligned one-to-one to all 38 samples")
}

variance_df <- data.frame(
  Component = c("Stromal Score", "Immune Score", "Residual"),
  Percent = c(45.3782101425, 31.3530572619, 23.2687325956),
  stringsAsFactors = FALSE
)

# Recalculate, rather than hard-code, the values used in the panel. This is
# deliberately the same sequential R2 definition as 08_mgs_variance_partitioning_estimate.R.
df_var <- data.frame(
  MGS = unname(reference_mgs),
  Stromal_Score = as.numeric(estimate_scores$Stromal_Score),
  Immune_Score = as.numeric(estimate_scores$Immune_Score)
)
r2_stromal <- summary(lm(MGS ~ Stromal_Score, data = df_var))$r.squared
r2_joint <- summary(lm(MGS ~ Stromal_Score + Immune_Score, data = df_var))$r.squared
variance_df$Percent <- c(100 * r2_stromal, 100 * (r2_joint - r2_stromal), 100 * (1 - r2_joint))
if (max(abs(variance_df$Percent - c(45.3782101425, 31.3530572619, 23.2687325956))) > 1e-6) {
  stop("Official ESTIMATE variance values do not match the audited rerun")
}
variance_df$Component <- factor(variance_df$Component,
                                levels = c("Stromal Score", "Immune Score", "Residual"))
write.csv(variance_df, file.path(out_dir, "figure4h_variance_partitioning.csv"),
          row.names = FALSE)
saveRDS(list(
  definition = "Sequential R2: MGS~Stromal_Score, then MGS~Stromal_Score+Immune_Score; residual=1-joint R2",
  input = df_var, components = variance_df
), file.path(out_dir, "figure4h_variance_partitioning.rds"))

# ----------------------------------------------------------------------------
# Draw publication panels using R only.
# ----------------------------------------------------------------------------
theme_panel <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, colour = "black"),
    plot.title = element_text(size = 9, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7, colour = "#444444", hjust = 0),
    plot.margin = margin(5.5, 6, 5.5, 6, unit = "pt")
  )

hist_df <- data.frame(kendall_tau = valid_tau)
p4e <- ggplot(hist_df, aes(x = kendall_tau)) +
  geom_histogram(bins = 30, boundary = -1, closed = "left",
                 fill = "#55B7C5", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = tau_mean, colour = "#D55E00", linewidth = 0.65,
             linetype = "dashed") +
  geom_vline(xintercept = tau_ci, colour = "#555555", linewidth = 0.45,
             linetype = "dotted") +
  annotate("label", x = 0.05, y = Inf,
           label = paste0("Mean tau = ", sprintf("%.3f", tau_mean),
                          "\n95% interval [", sprintf("%.3f", tau_ci[1]), ", ",
                          sprintf("%.3f", tau_ci[2]), "]"),
           hjust = 0, vjust = 1.12, size = 2.7,
           fill = "white", colour = "#333333", label.padding = unit(0.13, "lines")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "e  Bootstrap stability (500 replicates)",
       x = "Kendall's tau with full-data MGS ordering", y = "Bootstrap replicates") +
  theme_panel

variance_plot_df <- variance_df %>%
  mutate(Component = factor(Component,
                            levels = c("Stromal Score", "Immune Score", "Residual")))
component_cols <- c("Stromal Score" = "#4472C4", "Immune Score" = "#ED7D31",
                    "Residual" = "#A6A6A6")
p4h <- ggplot(variance_plot_df, aes(x = "MGS", y = Percent, fill = Component)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_stack(vjust = 0.5), size = 3.0,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = component_cols, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "h  Hierarchical variance partitioning",
       subtitle = "Official ESTIMATE component scores",
       x = NULL, y = "Variance explained (%)", fill = NULL) +
  theme_panel +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.key.height = unit(0.42, "cm"), plot.subtitle = element_text(size = 7))

save_panel <- function(plot, stem, width_mm = 89, height_mm = 72) {
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_file <- file.path(figure_dir, paste0(stem, ".png"))
  tiff_file <- file.path(figure_dir, paste0(stem, ".tiff"))
  svg_file <- file.path(figure_dir, paste0(stem, ".svg"))
  grDevices::pdf(pdf_file, width = width_mm / 25.4,
                 height = height_mm / 25.4, family = "Helvetica")
  print(plot)
  dev.off()
  ragg::agg_tiff(tiff_file, width = width_mm / 25.4,
                 height = height_mm / 25.4, units = "in", res = 600,
                 compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(png_file, width = width_mm / 25.4 * 300,
                height = height_mm / 25.4 * 300, res = 300)
  print(plot)
  dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(svg_file, width = width_mm / 25.4, height = height_mm / 25.4)
    print(plot)
    dev.off()
  }
  invisible(c(pdf = pdf_file, png = png_file, tiff = tiff_file, svg = svg_file))
}

save_panel(p4e, "Figure4e_updated")
save_panel(p4h, "Figure4h_updated")

combined_file <- file.path(figure_dir, "Figure4e_4h_updated_panels.pdf")
grDevices::pdf(combined_file, width = 178 / 25.4, height = 78 / 25.4,
               family = "Helvetica")
gridExtra::grid.arrange(p4e, p4h, ncol = 2, widths = c(1, 1))
dev.off()

# Also place copies in the analysis output directory for the reproducibility
# bundle while retaining the same R-generated visual files in the figure folder.
file.copy(file.path(figure_dir, "Figure4e_updated.pdf"), out_dir,
          overwrite = TRUE)
file.copy(file.path(figure_dir, "Figure4h_updated.pdf"), out_dir,
          overwrite = TRUE)
for (panel_stem in c("Figure4e_updated", "Figure4h_updated")) {
  for (ext in c("svg", "png", "tiff")) {
    file.copy(file.path(figure_dir, paste0(panel_stem, ".", ext)), out_dir,
              overwrite = TRUE)
  }
}
file.copy(combined_file, out_dir, overwrite = TRUE)

# ----------------------------------------------------------------------------
# Compact replacement panels matching the user's manually revised Figure 4
# typography. These panels intentionally omit the standalone plot titles in
# Figure 4e and use the compact axis range shown in the submitted layout.
# ----------------------------------------------------------------------------
compact_theme <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, colour = "black"),
    plot.margin = margin(4, 5, 4, 5, unit = "pt")
  )

p4e_final <- ggplot(hist_df, aes(x = kendall_tau)) +
  geom_histogram(bins = 30, boundary = 0.35, closed = "left",
                 fill = "#56B5CC", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = tau_mean, colour = "#56B5CC", linewidth = 0.65,
             linetype = "dashed") +
  geom_vline(xintercept = tau_ci, colour = "#56B5CC", linewidth = 0.45,
             linetype = "dotted") +
  annotate("text", x = 0.40, y = 74,
           label = paste0("Mean tau=", sprintf("%.3f", tau_mean),
                          "\n95% CI [", sprintf("%.3f", tau_ci[1]), ", ",
                          sprintf("%.3f", tau_ci[2]), "]"),
           hjust = 0, vjust = 1, size = 2.9, colour = "#56B5CC",
           lineheight = 0.95) +
  annotate("text", x = 0.35, y = 74, label = "e",
           hjust = 0, vjust = 1.15, size = 3.2,
           fontface = "bold", colour = "black") +
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1.0),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = c(0, 20, 40, 60),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Kendall tau", y = "Count") +
  coord_cartesian(xlim = c(0.35, 1.02), ylim = c(0, 75), clip = "on") +
  compact_theme

p4h_final <- ggplot(variance_plot_df, aes(x = "MGS", y = Percent, fill = Component)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_stack(vjust = 0.5), size = 2.8,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = component_cols, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "h  Hierarchical variance partitioning",
       x = NULL, y = "Variance explained (%)", fill = NULL) +
  compact_theme +
  theme(legend.position = "right", legend.text = element_text(size = 6.8),
        legend.key.height = unit(0.36, "cm"),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0))

save_compact_panel <- function(plot, stem, width_mm, height_mm) {
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(figure_dir, paste0(stem, ".svg"))
  tiff_file <- file.path(figure_dir, paste0(stem, ".tiff"))
  png_file <- file.path(figure_dir, paste0(stem, ".png"))
  grDevices::pdf(pdf_file, width = width_mm / 25.4,
                 height = height_mm / 25.4, family = "Helvetica")
  print(plot)
  dev.off()
  svglite::svglite(svg_file, width = width_mm / 25.4,
                   height = height_mm / 25.4)
  print(plot)
  dev.off()
  ragg::agg_tiff(tiff_file, width = width_mm / 25.4,
                 height = height_mm / 25.4, units = "in", res = 600,
                 compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(png_file, width = width_mm / 25.4 * 300,
                height = height_mm / 25.4 * 300, res = 300)
  print(plot)
  dev.off()
  for (f in c(pdf_file, svg_file, tiff_file, png_file)) {
    file.copy(f, out_dir, overwrite = TRUE)
  }
}

save_compact_panel(p4e_final, "Figure4e_panel_final", 76, 62)
save_compact_panel(p4h_final, "Figure4h_panel_final", 89, 72)

final_combined <- file.path(figure_dir, "Figure4e_4h_panel_final.pdf")
grDevices::pdf(final_combined, width = 168 / 25.4, height = 72 / 25.4,
               family = "Helvetica")
gridExtra::grid.arrange(p4e_final, p4h_final, ncol = 2,
                        widths = c(0.88, 1.12))
dev.off()
file.copy(final_combined, out_dir, overwrite = TRUE)

cat("Figure 4e valid bootstrap replicates: ", length(valid_tau), " / ", N_BOOT, "\n", sep = "")
cat("Figure 4e mean Kendall tau: ", sprintf("%.6f", tau_mean),
    "; 95% percentile interval [", sprintf("%.6f", tau_ci[1]), ", ",
    sprintf("%.6f", tau_ci[2]), "]\n", sep = "")
cat("Figure 4h components: ", paste(sprintf("%s=%.6f%%", as.character(variance_df$Component),
                                               variance_df$Percent), collapse = "; "), "\n", sep = "")
cat("Panels written to: ", normalizePath(figure_dir, mustWork = TRUE), "\n", sep = "")
