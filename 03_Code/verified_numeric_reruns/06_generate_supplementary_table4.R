# Generate Supplementary Table 4 directly from the de-identified clinical
# workbook used for the discovery-cohort regression analysis.
#
# Required environment variable:
#   VS_CLINICAL_XLSX       path to cl_bulk_with_Subtype.xlsx
# Optional:
#   VS_SUPPLEMENT_DIR      output directory for the CSV table

options(stringsAsFactors = FALSE)

clinical_path <- Sys.getenv("VS_CLINICAL_XLSX", unset = "")
if (!nzchar(clinical_path) || !file.exists(clinical_path)) {
  stop("Set VS_CLINICAL_XLSX to the clinical workbook used for the discovery cohort.")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else getwd()
package_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
out_dir <- Sys.getenv("VS_SUPPLEMENT_DIR", unset = file.path(package_root, "02_Submission_Materials", "03_Supplementary_Tables"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(car)
  library(lmtest)
})

cl <- readxl::read_excel(clinical_path)
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", , drop = FALSE]
cl <- cl[!is.na(cl$Age), , drop = FALSE]
stopifnot(nrow(cl) == 38L)

# The clinical workbook records this laboratory variable as a percentage.
cl$NK_percentage <- as.numeric(cl$自然杀伤细胞)
cl$Age <- as.numeric(cl$Age)
cl$Subtype <- factor(cl$Subtype, levels = c("C1", "C2", "C3"))
cl$SizeLarge <- ifelse(cl$SizeGrade == "large", 1, 0)
stopifnot(sum(!is.na(cl$NK_percentage)) == 38L)

m1 <- lm(NK_percentage ~ Subtype + Age, data = cl)
m2 <- lm(NK_percentage ~ Subtype + Age + SizeLarge, data = cl)

adjusted_gvif_max <- function(model) {
  v <- car::vif(model)
  if (is.matrix(v)) max(v[, ncol(v)], na.rm = TRUE) else max(v, na.rm = TRUE)
}

model_rows <- function(model, model_label, predictors) {
  co <- summary(model)$coefficients
  ci <- confint(model, level = 0.95)
  bp <- lmtest::bptest(model)$p.value
  shapiro <- shapiro.test(residuals(model))$p.value
  diagnostics <- list(
    r2 = summary(model)$r.squared,
    shapiro = shapiro,
    bp = bp,
    gvif = adjusted_gvif_max(model),
    cook = max(cooks.distance(model), na.rm = TRUE)
  )
  out <- lapply(seq_along(predictors), function(i) {
    predictor <- names(predictors)[i]
    term <- unname(predictors[[i]])
    data.frame(
      Model = model_label,
      Predictor = predictor,
      Beta = unname(co[term, "Estimate"]),
      SE = unname(co[term, "Std. Error"]),
      `95% CI` = sprintf("%.3f to %.3f", ci[term, 1], ci[term, 2]),
      `P value` = unname(co[term, "Pr(>|t|)"]),
      `Model R2` = if (i == 1) diagnostics$r2 else NA_real_,
      `Shapiro-Wilk P` = if (i == 1) diagnostics$shapiro else NA_real_,
      `Breusch-Pagan P` = if (i == 1) diagnostics$bp else NA_real_,
      `Adjusted GVIF max` = if (i == 1) diagnostics$gvif else NA_real_,
      `Cook's D max` = if (i == 1) diagnostics$cook else NA_real_,
      check.names = FALSE
    )
  })
  do.call(rbind, out)
}

rows_m1 <- model_rows(
  m1,
  "Age-adjusted (NK percentage ~ Subtype + Age)",
  c("C3 (vs C1)" = "SubtypeC3", "Age" = "Age", "C2 (vs C1)" = "SubtypeC2")
)
rows_m2 <- model_rows(
  m2,
  "Fully-adjusted (NK percentage ~ Subtype + Age + Tumour Size)",
  c("C3 (vs C1)" = "SubtypeC3", "Age" = "Age", "Tumour Size (large vs small)" = "SizeLarge")
)

result <- rbind(rows_m1, rows_m2)
write.csv(result, file.path(out_dir, "Supplementary_Table4_Regression.csv"), row.names = FALSE, na = "")
writeLines(capture.output({
  cat("Supplementary Table 4 generated from clinical workbook\n")
  cat("n = ", nrow(cl), "\n", sep = "")
  cat("VIF field = maximum adjusted GVIF^(1/(2*df)); the factor df is accounted for.\n\n")
  print(result, row.names = FALSE)
}), file.path(out_dir, "Supplementary_Table4_Regression_generation.txt"))
print(result, row.names = FALSE)
