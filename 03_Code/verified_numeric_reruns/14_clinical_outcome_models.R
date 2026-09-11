# Outcome-specific clinical models used in the methodological audit.
#
# This script makes the model distinction in the manuscript executable:
#   * tumour texture: proportional-odds ordinal logistic regression;
#   * brainstem adhesion: Firth penalized logistic regression;
#   * recorded NK-cell percentage: ordinary linear models.
#
# It uses the de-identified clinical workbook and the separate surgery-
# assessment workbook. The raw FASTQ and single-cell preprocessing branches
# are unrelated to this script and are intentionally not called here.
#
# Required environment variables:
#   VS_CLINICAL_XLSX  path to cl_bulk_with_Subtype.xlsx
#   VS_SURGERY_XLSX   path to VS surgery analysis.xlsx
# Optional:
#   VS_OUTPUT_DIR     audit output directory (from 00_config.R)

options(stringsAsFactors = FALSE)

config_path <- Sys.getenv("VS_CONFIG_R", unset = "")
if (!nzchar(config_path) || !file.exists(config_path)) {
  candidates <- c(
    file.path(getwd(), "00_config.R"),
    file.path(dirname(normalizePath(commandArgs()[1], mustWork = FALSE)), "00_config.R")
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Cannot locate 00_config.R; set VS_CONFIG_R.")
  config_path <- candidates[1]
}
source(config_path)

clinical_path <- Sys.getenv("VS_CLINICAL_XLSX", unset = CLINICAL_XLSX)
surgery_path <- Sys.getenv("VS_SURGERY_XLSX", unset = SURGERY_XLSX)
if (!file.exists(clinical_path)) stop("Missing clinical workbook: ", clinical_path)
if (!file.exists(surgery_path)) stop("Missing surgery workbook: ", surgery_path)

out_dir <- file.path(AUDIT_ROOT, "clinical_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(MASS)
  library(logistf)
})

raw <- as.data.frame(readxl::read_excel(clinical_path))
surgery <- as.data.frame(readxl::read_excel(surgery_path))
required_clinical <- c("SampleID", "Subtype", "Age", "TumorSize_raw", "SizeGrade", "自然杀伤细胞")
required_surgery <- c("SampleID", "肿瘤质地", "脑干面粘连")
if (any(!required_clinical %in% names(raw))) {
  stop("Clinical workbook is missing: ", paste(setdiff(required_clinical, names(raw)), collapse = ", "))
}
if (any(!required_surgery %in% names(surgery))) {
  stop("Surgery workbook is missing: ", paste(setdiff(required_surgery, names(surgery)), collapse = ", "))
}

cl <- raw[raw$Subtype %in% c("C1", "C2", "C3"), , drop = FALSE]
su <- surgery
stopifnot(nrow(cl) == 38L, !anyDuplicated(cl$SampleID), !anyDuplicated(su$SampleID))
stopifnot(setequal(cl$SampleID, su$SampleID))
su <- su[match(cl$SampleID, su$SampleID), , drop = FALSE]
stopifnot(identical(as.character(su$SampleID), as.character(cl$SampleID)))

parse_size <- function(x) {
  values <- as.numeric(strsplit(gsub("×", "*", as.character(x), fixed = TRUE), "\\*")[[1]])
  stopifnot(length(values) == 3L, !anyNA(values))
  max(values)
}

size_mm <- vapply(cl$TumorSize_raw, parse_size, numeric(1))
# Author-confirmed correction: P32's first dimension in 2.2*13*14 is 22 mm.
size_mm[cl$SampleID == "P32"] <- 22

texture_map <- c("软" = "soft", "中等" = "medium", "硬" = "hard")
adhesion_map <- c("有" = 1, "无" = 0)
df <- data.frame(
  SampleID = cl$SampleID,
  Subtype = factor(cl$Subtype, levels = c("C1", "C2", "C3")),
  Age = as.numeric(cl$Age),
  size_raw = as.character(cl$TumorSize_raw),
  size_mm = size_mm,
  size_cm = size_mm / 10,
  size_provenance = ifelse(cl$SampleID == "P32", "author_confirmed_20260906", "parsed_original_mm"),
  texture = ordered(unname(texture_map[as.character(su$肿瘤质地)]), levels = c("soft", "medium", "hard")),
  adhesion = unname(adhesion_map[as.character(su$脑干面粘连)]),
  NK_percentage = as.numeric(cl$自然杀伤细胞),
  size_large = as.integer(cl$SizeGrade == "large"),
  check.names = FALSE
)
stopifnot(!anyNA(df), sum(df$adhesion) == 10L)

write.csv(
  data.frame(
    SampleID = "P32",
    raw_record = as.character(cl$TumorSize_raw[cl$SampleID == "P32"]),
    corrected_first_dimension_mm = 22,
    correction_basis = "author_confirmed_first_dimension",
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "clinical_size_correction.csv"), row.names = FALSE
)

coefficient_rows <- list()
joint_rows <- list()
model_cache <- list()

add_model <- function(model, label, outcome, kind) {
  if (kind == "firth") {
    beta <- coef(model)
    se <- sqrt(diag(model$var))
    low <- model$ci.lower
    high <- model$ci.upper
    p_value <- model$prob
    covariance <- model$var
    ci_method <- "penalized_profile_likelihood"
    estimate <- exp(beta)
    ci_low <- exp(low)
    ci_high <- exp(high)
  } else {
    beta <- coef(model)
    covariance <- vcov(model)[names(beta), names(beta), drop = FALSE]
    se <- sqrt(diag(covariance))[names(beta)]
    if (kind == "ordinal") {
      interval <- tryCatch(
        suppressMessages(confint(model)),
        error = function(e) NULL
      )
      if (is.null(interval) || any(!is.finite(interval))) {
        low <- beta - 1.96 * se
        high <- beta + 1.96 * se
        ci_method <- "Wald_fallback"
      } else {
        low <- interval[names(beta), 1]
        high <- interval[names(beta), 2]
        ci_method <- "profile_likelihood"
      }
      p_value <- 2 * pnorm(-abs(beta / se))
      estimate <- exp(beta)
      ci_low <- exp(low)
      ci_high <- exp(high)
    } else {
      interval <- confint(model)
      low <- interval[names(beta), 1]
      high <- interval[names(beta), 2]
      p_value <- summary(model)$coefficients[names(beta), 4]
      ci_method <- "t_interval"
      estimate <- beta
      ci_low <- low
      ci_high <- high
    }
  }

  subtype_terms <- grep("^Subtype", names(beta))
  subtype_cov <- covariance[subtype_terms, subtype_terms, drop = FALSE]
  subtype_beta <- beta[subtype_terms]
  joint_stat <- as.numeric(t(subtype_beta) %*% solve(subtype_cov, subtype_beta))
  joint_p <- pchisq(joint_stat, length(subtype_terms), lower.tail = FALSE)

  coefficient_rows[[label]] <<- data.frame(
    model = label,
    outcome = outcome,
    method = kind,
    n = nrow(df),
    term = names(beta),
    beta = unname(beta),
    SE = unname(se),
    estimate = unname(estimate),
    CI_low = unname(ci_low),
    CI_high = unname(ci_high),
    p = unname(p_value),
    CI_method = ci_method,
    stringsAsFactors = FALSE
  )
  joint_rows[[label]] <<- data.frame(
    model = label,
    outcome = outcome,
    method = kind,
    n = nrow(df),
    joint_subtype_Wald = joint_stat,
    df = length(subtype_terms),
    p = joint_p,
    stringsAsFactors = FALSE
  )
  model_cache[[label]] <<- model
}

for (label in c("unadjusted", "size_adjusted", "size_age_sensitivity")) {
  rhs <- switch(
    label,
    unadjusted = "Subtype",
    size_adjusted = "Subtype + size_cm",
    size_age_sensitivity = "Subtype + size_cm + Age"
  )
  firth_fit <- logistf::logistf(
    as.formula(paste("adhesion ~", rhs)), data = df, pl = TRUE,
    control = logistf::logistf.control(maxit = 1000),
    plcontrol = logistf::logistpl.control(maxit = 1000)
  )
  add_model(firth_fit, paste0("adhesion_", label), "brainstem_adhesion", "firth")

  ordinal_fit <- MASS::polr(
    as.formula(paste("texture ~", rhs)), data = df,
    Hess = TRUE, method = "logistic"
  )
  add_model(ordinal_fit, paste0("texture_", label), "texture_soft_medium_hard", "ordinal")
}

for (label in c("age_adjusted", "age_size_binary_legacy", "age_size_continuous")) {
  rhs <- switch(
    label,
    age_adjusted = "Subtype + Age",
    age_size_binary_legacy = "Subtype + Age + size_large",
    age_size_continuous = "Subtype + Age + size_cm"
  )
  linear_fit <- lm(as.formula(paste("NK_percentage ~", rhs)), data = df)
  add_model(linear_fit, paste0("NK_", label), "NK_recorded_laboratory_percentage", "linear")
}

coefficients <- do.call(rbind, coefficient_rows)
omnibus <- do.call(rbind, joint_rows)
omnibus$q_primary_two_surgical_endpoints <- NA_real_
primary <- omnibus$model %in% c("adhesion_size_adjusted", "texture_size_adjusted")
omnibus$q_primary_two_surgical_endpoints[primary] <- p.adjust(omnibus$p[primary], method = "BH")

write.csv(coefficients, file.path(out_dir, "clinical_model_coefficients.csv"), row.names = FALSE)
write.csv(omnibus, file.path(out_dir, "clinical_omnibus_tests.csv"), row.names = FALSE)
saveRDS(model_cache, file.path(out_dir, "clinical_model_fits.rds"))

# Proportional-odds sensitivity: compare subtype effects at the two cutpoints.
threshold_rows <- list()
for (threshold in c(1, 2)) {
  z <- df
  z$harder <- as.integer(as.integer(z$texture) > threshold)
  fit <- logistf::logistf(harder ~ Subtype + size_cm, data = z, pl = TRUE)
  threshold_rows[[as.character(threshold)]] <- data.frame(
    threshold = threshold,
    term = names(coef(fit)),
    OR = exp(coef(fit)),
    CI_low = exp(fit$ci.lower),
    CI_high = exp(fit$ci.upper),
    p = fit$prob,
    stringsAsFactors = FALSE
  )
}
write.csv(
  do.call(rbind, threshold_rows),
  file.path(out_dir, "texture_threshold_sensitivity.csv"), row.names = FALSE
)

writeLines(
  c(
    "Outcome-specific clinical models",
    paste0("n = ", nrow(df)),
    "Tumour texture: proportional-odds ordinal logistic regression (MASS::polr).",
    "Brainstem adhesion: Firth penalized logistic regression (logistf).",
    "NK-cell laboratory variable: recorded percentage; ordinary linear models.",
    "P32 first dimension: author-confirmed 22 mm.",
    "The script does not treat these models as a single generic linear-regression analysis."
  ),
  file.path(out_dir, "clinical_models_generation.txt")
)

print(coefficients[grepl("size_adjusted|age_size_continuous", coefficients$model) & grepl("Subtype", coefficients$term), ])
print(omnibus)
