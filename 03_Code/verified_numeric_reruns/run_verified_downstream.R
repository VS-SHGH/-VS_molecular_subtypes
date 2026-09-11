# Run the verified downstream analyses in dependency order.
# This runner intentionally excludes raw FASTQ/10X preprocessing, which was
# not independently reconstructable from the supplied files.

options(stringsAsFactors = FALSE)
runner_arg <- grep("^--file=", commandArgs(), value = TRUE)
runner_path <- if (length(runner_arg)) sub("^--file=", "", runner_arg[1]) else "run_verified_downstream.R"
runner_dir <- dirname(normalizePath(runner_path, mustWork = FALSE))
scripts <- c(
  "06_generate_supplementary_table4.R",
  "14_clinical_outcome_models.R",
  "01_primary_clustering_rerun.R",
  "11_estimate_official_rerun.R",
  "07_mgs_unbiased_immune_anchoring.R",
  "08_mgs_variance_partitioning_estimate.R",
  "09_external_validation_rerun.R",
  "10_cibersort_mcp_correlation.R",
  "03_singlecell_downstream_rerun.R",
  "04_myeloid_deg_gsea_rerun.R",
  "05_cellchat_numeric_rerun.R"
)
for (script in scripts) {
  message("Running ", script)
  source(file.path(runner_dir, script), local = new.env(parent = globalenv()))
}
message("Verified downstream rerun sequence completed.")
