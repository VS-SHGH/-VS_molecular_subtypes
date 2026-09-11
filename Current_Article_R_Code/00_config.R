# Portable configuration for the verified downstream reruns.
# Set the VS_* environment variables when running outside the prepared package.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else getwd()
code_dir <- if (dir.exists(script_path)) normalizePath(script_path) else dirname(normalizePath(script_path))
package_root <- Sys.getenv("VS_UPLOAD_ROOT", unset = normalizePath(file.path(code_dir, "..", ".."), mustWork = FALSE))

AUX_R_LIBRARY <- Sys.getenv("VS_AUX_R_LIBRARY", unset = "")
if (nzchar(AUX_R_LIBRARY) && dir.exists(AUX_R_LIBRARY)) {
  .libPaths(unique(c(AUX_R_LIBRARY, .libPaths())))
}

AUDIT_ROOT <- Sys.getenv("VS_OUTPUT_DIR", unset = file.path(package_root, "outputs"))
SERVER_RDATA_DIR <- Sys.getenv("VS_RDATA_DIR", unset = file.path(package_root, "input", "RData"))
BULK_TPM_XLSX <- Sys.getenv("VS_BULK_TPM_XLSX", unset = file.path(package_root, "input", "TPM.xlsx"))
CLINICAL_XLSX <- Sys.getenv("VS_CLINICAL_XLSX", unset = file.path(package_root, "input", "cl_bulk_with_Subtype.xlsx"))
SURGERY_XLSX <- Sys.getenv("VS_SURGERY_XLSX", unset = file.path(package_root, "input", "VS_surgery_analysis.xlsx"))
SC_QS_INPUT <- Sys.getenv("VS_SC_QS", unset = file.path(package_root, "input", "sc_nc_raw.qs"))
CELLCHAT_INPUT_RDS <- Sys.getenv("VS_CELLCHAT_INPUT_RDS", unset = file.path(SERVER_RDATA_DIR, "Step07_sc_for_cellchat.rds"))
CELLCHAT_CACHE_RDS <- Sys.getenv("VS_CELLCHAT_CACHE_RDS", unset = file.path(SERVER_RDATA_DIR, "Step08_cellchat_object.rds"))
LEGACY_FIGURE_DIR <- Sys.getenv("VS_LEGACY_FIGURE_DIR", unset = file.path(package_root, "input", "Figures"))
STEP3_RDATA <- Sys.getenv("VS_STEP3_RDATA", unset = file.path(SERVER_RDATA_DIR, "Step3_DEA_Results.RData"))
ESTIMATE_SCORES_CSV <- Sys.getenv("VS_ESTIMATE_SCORES_CSV", unset = file.path(AUDIT_ROOT, "estimate_official", "estimate_official_scores.csv"))
EXTERNAL_GSE141801_EXPR <- Sys.getenv("VS_EXTERNAL_GSE141801_EXPR", unset = file.path(package_root, "input", "GSE141801_Expression_Log2.csv"))
EXTERNAL_GSE141801_CLIN <- Sys.getenv("VS_EXTERNAL_GSE141801_CLIN", unset = file.path(package_root, "input", "GSE141801_Clinical_Data.csv"))
EXTERNAL_GSE39645_EXPR <- Sys.getenv("VS_EXTERNAL_GSE39645_EXPR", unset = file.path(package_root, "input", "GSE39645_Expression_Log2.csv"))
EXTERNAL_GSE39645_CLIN <- Sys.getenv("VS_EXTERNAL_GSE39645_CLIN", unset = file.path(package_root, "input", "GSE39645_Clinical_Data.csv"))
CIBERSORT_CSV <- Sys.getenv("VS_CIBERSORT_CSV", unset = file.path(package_root, "input", "CIBERSORT_Result.csv"))
SUBMISSION_TABLE3_CSV <- Sys.getenv("VS_SUBMISSION_TABLE3_CSV", unset = file.path(package_root, "02_Submission_Materials", "03_Supplementary_Tables", "Supplementary_Table3_CrossValidation.csv"))
GSVA_COMPAT_LIBRARY <- Sys.getenv("VS_GSVA_COMPAT_LIBRARY", unset = "")

dir.create(AUDIT_ROOT, recursive = TRUE, showWarnings = FALSE)
