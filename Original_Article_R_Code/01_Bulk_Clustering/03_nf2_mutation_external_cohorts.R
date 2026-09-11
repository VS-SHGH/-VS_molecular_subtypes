# =============================================================================
# Revision Analysis 3: NF2 Mutation Analysis in External Cohorts
# =============================================================================

library(dplyr)

# ---- Load external cohort clinical data ----
gse141801_clin <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE141801_Gugel/Clinical_Data.csv")
gse39645_clin <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE39645_Torres/Clinical_Data.csv")

cat("=== GSE141801 (Gugel et al.) ===\n")
cat(sprintf("Total samples: %d\n", nrow(gse141801_clin)))
cat(sprintf("NF2 column available: %s\n", "NF2" %in% colnames(gse141801_clin)))
cat("NF2 distribution:\n")
print(table(gse141801_clin$NF2, useNA = "ifany"))

# Load external expression data and assign subtypes
gse141801_expr <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE141801_Gugel/GSE141801_Expression_Log2.csv",
                            row.names = 1, check.names = FALSE)
gse39645_expr <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE39645_Torres/GSE39645_Expression_Log2.csv",
                           row.names = 1, check.names = FALSE)

# The clinical data has Class column (Injury-like / nmSC Core) from Barrett et al.
# This is the predicted subtype from the original analysis
cat("\n=== GSE141801 Class (Barrett classification) ===\n")
print(table(gse141801_clin$Class, useNA = "ifany"))

cat("\n=== GSE39645 Class ===\n")
print(table(gse39645_clin$Class, useNA = "ifany"))

# ---- NF2 by Class (Barrett classification) ----
cat("\n========== NF2 Mutation by Barrett Class ==========\n")

# GSE141801
nf2_141801 <- gse141801_clin[!is.na(gse141801_clin$NF2) & gse141801_clin$NF2 != "", ]
cat(sprintf("\nGSE141801: %d samples with NF2 data\n", nrow(nf2_141801)))
tbl_141801 <- table(nf2_141801$Class, nf2_141801$NF2)
cat("NF2 by Class (GSE141801):\n")
print(tbl_141801)

# Fisher's exact test
if (nrow(tbl_141801) >= 2 && ncol(tbl_141801) >= 2) {
  ft_141801 <- fisher.test(tbl_141801, simulate.p.value = TRUE, B = 10000)
  cat(sprintf("Fisher's exact P = %.4f\n", ft_141801$p.value))
}

# Per-class NF2 mutation rate
cat("\nNF2 mutation rates by Class (GSE141801):\n")
for (cl in unique(nf2_141801$Class)) {
  sub <- nf2_141801[nf2_141801$Class == cl, ]
  n_mut <- sum(sub$NF2 == "Y")
  cat(sprintf("  %s: %d/%d (%.1f%%)\n", cl, n_mut, nrow(sub), 100*n_mut/nrow(sub)))
}

# GSE39645
nf2_39645 <- gse39645_clin[!is.na(gse39645_clin$NF2) & gse39645_clin$NF2 != "", ]
cat(sprintf("\nGSE39645: %d samples with NF2 data\n", nrow(nf2_39645)))
tbl_39645 <- table(nf2_39645$Class, nf2_39645$NF2)
cat("NF2 by Class (GSE39645):\n")
print(tbl_39645)

if (nrow(tbl_39645) >= 2 && ncol(tbl_39645) >= 2) {
  ft_39645 <- fisher.test(tbl_39645, simulate.p.value = TRUE, B = 10000)
  cat(sprintf("Fisher's exact P = %.4f\n", ft_39645$p.value))
}

cat("\nNF2 mutation rates by Class (GSE39645):\n")
for (cl in unique(nf2_39645$Class)) {
  sub <- nf2_39645[nf2_39645$Class == cl, ]
  n_mut <- sum(sub$NF2 == "Y")
  cat(sprintf("  %s: %d/%d (%.1f%%)\n", cl, n_mut, nrow(sub), 100*n_mut/nrow(sub)))
}

# ---- Combined analysis ----
cat("\n========== Combined External Cohorts ==========\n")
nf2_combined <- rbind(
  gse141801_clin[, c("NF2", "Class")],
  gse39645_clin[, c("NF2", "Class")]
)
nf2_combined <- nf2_combined[!is.na(nf2_combined$NF2) & nf2_combined$NF2 != "", ]
cat(sprintf("Combined: %d samples with NF2 data\n", nrow(nf2_combined)))

tbl_comb <- table(nf2_combined$Class, nf2_combined$NF2)
print(tbl_comb)
ft_comb <- fisher.test(tbl_comb, simulate.p.value = TRUE, B = 10000)
cat(sprintf("Fisher's exact P = %.4f\n", ft_comb$p.value))

cat("\nNF2 mutation rates (combined):\n")
for (cl in unique(nf2_combined$Class)) {
  sub <- nf2_combined[nf2_combined$Class == cl, ]
  n_mut <- sum(sub$NF2 == "Y")
  cat(sprintf("  %s: %d/%d (%.1f%%)\n", cl, n_mut, nrow(sub), 100*n_mut/nrow(sub)))
}

# ---- Sensitivity analysis: NF2 detection limitations ----
cat("\n========== Sensitivity Consideration ==========\n")
cat("NF2 inactivation in sporadic VS: 60-70% when including mutations + 22q LOH (Agnihotri 2016)\n")
cat("Microarray-based NF2 detection misses: copy-neutral LOH, promoter methylation, deep intronic mutations\n")
cat("Expected true NF2 inactivation rate: ~65%\n")
cat("Observed in external cohorts:\n")
total_nf2 <- sum(nf2_combined$NF2 == "Y")
total_samples <- nrow(nf2_combined)
cat(sprintf("  Observed: %d/%d = %.1f%%\n", total_nf2, total_samples, 100*total_nf2/total_samples))
cat(sprintf("  Estimated undetected: ~%.0f samples\n",
            0.65 * total_samples - total_nf2))

# ---- Size/NF2 relationship ----
cat("\n========== NF2 by Tumor Size (GSE141801) ==========\n")
gse141801_clin$NF2_bin <- ifelse(gse141801_clin$NF2 == "Y", "Mutant", "WT")
tbl_size <- table(gse141801_clin$Size, gse141801_clin$NF2_bin)
print(tbl_size)

cat("\n========== Analysis Complete ==========\n")
