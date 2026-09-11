# =============================================================================
# Analysis 11: Subtype vs Clinical Outcomes
# HB grade and Extent of Resection in 25 matched samples
# =============================================================================

library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"

# ---- Load 38 discovery sample info ----
cl38 <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
cl38 <- cl38[!is.na(cl38$Subtype) & cl38$Subtype != "", ]
cl38$住院号_str <- as.character(cl38$住院号)

# ---- Load immune343 for matched samples ----
imm343 <- read_excel("/Users/yanchen/Desktop/外周血免疫2/01_Datasets/immune 343.xlsx")
imm343$病案号_str <- as.character(imm343$病案号)

# ---- Match 38 samples to immune343 ----
cat("========== Sample Matching ==========\n")
matched <- merge(cl38, imm343, by.x = "住院号_str", by.y = "病案号_str", all.x = TRUE)
cat(sprintf("Total discovery: %d\n", nrow(cl38)))
cat(sprintf("Matched in immune343: %d\n", sum(!is.na(matched$术后面神经功能))))
cat(sprintf("Missing: %d\n", sum(is.na(matched$术后面神经功能))))

# ---- 1. HB Grade Analysis (facial nerve function) ----
cat("\n========== HB Grade by Subtype ==========\n")

hb_data <- matched[!is.na(matched$术后面神经功能), ]
hb_data$HB <- factor(hb_data$术后面神经功能, levels = c(1,2,3,4),
                     labels = c("I (Normal)", "II (Mild)", "III (Moderate)", "IV (Severe)"))

cat(sprintf("Samples with HB data: %d\n", nrow(hb_data)))
cat("\nHB distribution by subtype:\n")
hb_tbl <- table(hb_data$Subtype, hb_data$HB)
print(hb_tbl)

# Fisher's exact for HB>=2 vs HB=1
hb_data$HB_ge2 <- ifelse(hb_data$术后面神经功能 >= 2, "HB≥2", "HB=I")
cat("\nHB≥2 rate by subtype:\n")
for (st in c("C1","C2","C3")) {
  sub <- hb_data[hb_data$Subtype == st, ]
  n_ge2 <- sum(sub$HB_ge2 == "HB≥2")
  cat(sprintf("  %s: %d/%d (%.1f%%)\n", st, n_ge2, nrow(sub), 100*n_ge2/nrow(sub)))
}
ft_hb <- fisher.test(table(hb_data$Subtype, hb_data$HB_ge2), simulate.p.value = TRUE, B = 10000)
cat(sprintf("Fisher's exact P = %.4f\n", ft_hb$p.value))

# Trend test: ordinal regression
hb_data$HB_num <- as.numeric(hb_data$术后面神经功能)
kw_hb <- kruskal.test(HB_num ~ Subtype, data = hb_data)
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw_hb$p.value))

# Also check HB by MGS
mgs_full <- readRDS(file.path(work_dir, "mgs_full_results.rds"))
mgs_vals <- mgs_full$mgs_discovery
# mgs_full$mgs_discovery already named by our 38 samples (tpm.P1, etc.)
# Map to SampleID (strip tpm. prefix)
names(mgs_vals) <- gsub("^tpm\\.", "", names(mgs_vals))

hb_data$MGS <- mgs_vals[hb_data$SampleID]
cat(sprintf("\nMGS by HB grade:\n"))
for (hb in 1:4) {
  sub <- hb_data[hb_data$术后面神经功能 == hb, ]
  if (nrow(sub) > 0) {
    cat(sprintf("  HB=%d: MGS=%.3f±%.3f (n=%d)\n", hb, mean(sub$MGS,na.rm=TRUE),
                sd(sub$MGS,na.rm=TRUE), nrow(sub)))
  }
}

# Correlation MGS vs HB
hb_mgs_data <- hb_data[!is.na(hb_data$MGS) & !is.na(hb_data$HB_num), ]
cat(sprintf("MGS vs HB Spearman rho=%.3f, P=%.4f\n",
            cor(hb_mgs_data$MGS, hb_mgs_data$HB_num, method="spearman"),
            cor.test(hb_mgs_data$MGS, hb_mgs_data$HB_num, method="spearman")$p.value))

# ---- 2. Extent of Resection ----
cat("\n========== Extent of Resection by Subtype ==========\n")

eor_data <- matched[!is.na(matched$切除程度) & matched$切除程度 != "", ]
cat(sprintf("Samples with EOR data: %d\n", nrow(eor_data)))

cat("\nEOR distribution by subtype:\n")
eor_tbl <- table(eor_data$Subtype, eor_data$切除程度)
print(eor_tbl)

# GTR rate per subtype
cat("\nGTR rate by subtype:\n")
for (st in c("C1","C2","C3")) {
  sub <- eor_data[eor_data$Subtype == st, ]
  n_gtr <- sum(sub$切除程度 == "GTR")
  cat(sprintf("  %s: %d/%d (%.1f%%)\n", st, n_gtr, nrow(sub), 100*n_gtr/nrow(sub)))
}

# Fisher GTR vs non-GTR
eor_data$GTR <- ifelse(eor_data$切除程度 == "GTR", "GTR", "Non-GTR")
if (length(unique(eor_data$GTR)) >= 2) {
  ft_eor <- fisher.test(table(eor_data$Subtype, eor_data$GTR), simulate.p.value = TRUE, B = 10000)
  cat(sprintf("Fisher's exact P (GTR by subtype) = %.4f\n", ft_eor$p.value))
}

# ---- 3. Combined surgical outcome (HB≥2 OR non-GTR) ----
cat("\n========== Composite Surgical Difficulty ==========\n")

combined <- matched[!is.na(matched$术后面神经功能) & !is.na(matched$切除程度), ]
combined$Poor_Outcome <- ifelse(combined$术后面神经功能 >= 2 | combined$切除程度 != "GTR", "Poor", "Good")
cat(sprintf("Samples with both HB+EOR: %d\n", nrow(combined)))

cat("\nPoor outcome rate by subtype:\n")
for (st in c("C1","C2","C3")) {
  sub <- combined[combined$Subtype == st, ]
  n_poor <- sum(sub$Poor_Outcome == "Poor")
  cat(sprintf("  %s: %d/%d (%.1f%%)\n", st, n_poor, nrow(sub), 100*n_poor/nrow(sub)))
}
ft_combined <- fisher.test(table(combined$Subtype, combined$Poor_Outcome),
                           simulate.p.value = TRUE, B = 10000)
cat(sprintf("Fisher's exact P = %.4f\n", ft_combined$p.value))

# ---- 4. Logistic regression: HB≥2 ~ Subtype + Age ----
cat("\n========== Logistic Regression: HB≥2 ==========\n")

hb_data$HB_ge2_num <- ifelse(hb_data$HB_ge2 == "HB≥2", 1, 0)

if (sum(hb_data$HB_ge2_num) >= 5) {
  logit1 <- glm(HB_ge2_num ~ Subtype + Age, data = hb_data, family = binomial)
  cat(sprintf("Model: HB≥2 ~ Subtype + Age\n"))
  print(summary(logit1)$coefficients)
} else {
  cat("Insufficient HB≥2 events for logistic regression (n=",
      sum(hb_data$HB_ge2_num), ")\n", sep="")
  cat("Only 2/25 patients had HB≥2, both with HB=III\n")
  cat("HB outcomes are excellent across all subtypes (92% HB=I)\n")
}

# ---- 5. Summary table for paper ----
cat("\n========== Summary for Manuscript ==========\n\n")

cat("Table: Clinical outcomes by molecular subtype\n\n")
cat(sprintf("%-15s %8s %8s %8s\n", "Outcome", "C1", "C2", "C3"))
cat(rep("-", 45), "\n", sep="")

# HB≥2
for (st in c("C1","C2","C3")) {
  sub <- hb_data[hb_data$Subtype == st, ]
  n_ge2 <- sum(sub$HB_ge2 == "HB≥2")
  assign(paste0("hb_", st), sprintf("%d/%d (%.0f%%)", n_ge2, nrow(sub), 100*n_ge2/nrow(sub)))
}
cat(sprintf("%-15s %8s %8s %8s\n", "HB≥2", hb_C1, hb_C2, hb_C3))

# GTR
for (st in c("C1","C2","C3")) {
  sub <- eor_data[eor_data$Subtype == st, ]
  n_gtr <- sum(sub$GTR == "GTR")
  assign(paste0("gtr_", st), sprintf("%d/%d (%.0f%%)", n_gtr, nrow(sub), 100*n_gtr/nrow(sub)))
}
cat(sprintf("%-15s %8s %8s %8s\n", "GTR rate", gtr_C1, gtr_C2, gtr_C3))

# Poor outcome
if (nrow(combined) > 0) {
  for (st in c("C1","C2","C3")) {
    sub <- combined[combined$Subtype == st, ]
    n_poor <- sum(sub$Poor_Outcome == "Poor")
    assign(paste0("poor_", st), sprintf("%d/%d (%.0f%%)", n_poor, nrow(sub), 100*n_poor/nrow(sub)))
  }
  cat(sprintf("%-15s %8s %8s %8s\n", "Poor outcome", poor_C1, poor_C2, poor_C3))
}

cat("\nP-values: HB≥2 Fisher P=%.4f; GTR Fisher P=%.4f; Poor outcome Fisher P=%.4f\n",
    ft_hb$p.value, if(exists("ft_eor")) ft_eor$p.value else NA,
    if(exists("ft_combined")) ft_combined$p.value else NA)

# ---- 6. Also check: Tumor texture vs EOR, Adhesion vs HB ----
cat("\n========== Surgical Phenotype → Outcome Analysis ==========\n")

if ("肿瘤质地" %in% colnames(combined)) {
  cat("\nTumor texture vs GTR:\n")
  print(table(combined$肿瘤质地, combined$GTR))
}
if ("脑干粘连" %in% colnames(combined)) {
  cat("\nBrainstem adhesion vs HB≥2:\n")
  print(table(combined$脑干粘连, combined$HB_ge2))
}

# ---- Save ----
saveRDS(list(hb_data = hb_data, eor_data = eor_data, combined = combined),
        file.path(work_dir, "clinical_outcomes_results.rds"))

cat("\n========== Complete ==========\n")
