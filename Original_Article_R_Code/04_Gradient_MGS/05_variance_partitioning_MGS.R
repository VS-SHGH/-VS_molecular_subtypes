# =============================================================================
# Analysis 10: Variance Partitioning of MGS
# Addresses R2: "What proportion of MGS variance is compositional vs cell-intrinsic?"
# =============================================================================

library(dplyr)
library(readxl)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"

# ---- Load TPM and all derived data ----
tpm_raw <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/TPM.xlsx")
genes <- tpm_raw$gene_name
keep <- !duplicated(genes)
tpm_raw <- tpm_raw[keep, ]; genes <- genes[keep]

our_samples <- c('tpm.P1','tpm.P10','tpm.P11','tpm.P12','tpm.P13','tpm.P14','tpm.P15',
                 'tpm.P16','tpm.P17','tpm.P18','tpm.P19','tpm.P2','tpm.P20','tpm.P21',
                 'tpm.P22','tpm.P23','tpm.P24','tpm.P25','tpm.P26','tpm.P27','tpm.P28',
                 'tpm.P3','tpm.P30','tpm.P32','tpm.P34','tpm.P35','tpm.P36','tpm.P37',
                 'tpm.P38','tpm.P39','tpm.P4','tpm.P40','tpm.P41','tpm.P5','tpm.P6',
                 'tpm.P7','tpm.S5','tpm.S8')

tpm_mat <- as.matrix(tpm_raw[, our_samples])
rownames(tpm_mat) <- genes
log_tpm <- log2(tpm_mat + 1)

cl <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]
sample_ids_clean <- gsub("^tpm\\.", "", our_samples)
subtypes <- setNames(cl$Subtype, cl$SampleID)[sample_ids_clean]

# Load previously computed data
mcp <- readRDS(file.path(work_dir, "mcp_counter_results.rds"))
est_df <- readRDS(file.path(work_dir, "estimate_results.rds"))
mgs_full <- readRDS(file.path(work_dir, "mgs_full_results.rds"))
mgs_d <- mgs_full$mgs_discovery
comp_df <- readRDS(file.path(work_dir, "composition_analysis.rds"))

# ---- 1. Schwann cell score (already computed) ----
sc_genes <- c("SOX10","PLP1","EGR2","MPZ","MBP","PMP22","MAL","PRX","CDH19","L1CAM",
              "GAP43","NGFR","GFAP","S100B")
sc_avail <- intersect(sc_genes, rownames(log_tpm))
sc_score <- colMeans(log_tpm[sc_avail, ])
names(sc_score) <- colnames(log_tpm)

# ---- 2. Build comprehensive multi-variable dataset ----
df <- data.frame(
  Sample = our_samples,
  Subtype = subtypes,
  MGS = mgs_d,
  SC_Score = sc_score[our_samples],
  Purity = est_df$Purity,
  Immune = est_df$Immune,
  Stromal = est_df$Stromal,
  Fibroblasts = as.numeric(mcp["Fibroblasts", our_samples]),
  Monocytes = as.numeric(mcp["Monocytic lineage", our_samples]),
  Tcells = as.numeric(mcp["T cells", our_samples]),
  CD8_Tcells = as.numeric(mcp["CD8 T cells", our_samples]),
  Bcells = as.numeric(mcp["B lineage", our_samples]),
  NK_cells_mcp = as.numeric(mcp["NK cells", our_samples]),
  Neutrophils = as.numeric(mcp["Neutrophils", our_samples]),
  Endothelial = as.numeric(mcp["Endothelial cells", our_samples]),
  mDC = as.numeric(mcp["Myeloid dendritic cells", our_samples]),
  CytoLymph = as.numeric(mcp["Cytotoxic lymphocytes", our_samples])
)

# ---- 3. Hierarchical variance partitioning ----
cat("========== Hierarchical Variance Partitioning of MGS ==========\n\n")

# Model progression
m0 <- lm(MGS ~ 1, data = df)  # null
m1 <- lm(MGS ~ Purity, data = df)
m2 <- lm(MGS ~ Purity + Immune, data = df)
m3 <- lm(MGS ~ Purity + Immune + Fibroblasts, data = df)
m4 <- lm(MGS ~ Purity + Immune + Fibroblasts + Monocytes + Tcells + Bcells, data = df)

models <- list(
  "Null" = m0,
  "Purity only" = m1,
  "Purity + Immune" = m2,
  "Purity + Immune + Fibroblasts" = m3,
  "Full (Purity+Immune+Fibro+3cell)" = m4
)

cat("Stepwise R² accumulation:\n")
cat(sprintf("%-40s %6s %6s %s\n", "Model", "R²", "ΔR²", "Interpretation"))
cat(rep("-", 75), "\n", sep="")

prev_r2 <- 0
for (mname in names(models)) {
  m <- models[[mname]]
  r2 <- summary(m)$r.squared
  delta <- r2 - prev_r2
  if (mname == "Null") {
    interp <- "Baseline"
  } else if (mname == "Purity only") {
    interp <- "Compositional (tumor cell fraction)"
  } else if (mname == "Purity + Immune") {
    interp <- "Immune infiltration (beyond purity)"
  } else if (mname == "Purity + Immune + Fibroblasts") {
    interp <- "Stromal/ECM remodeling"
  } else {
    interp <- "Residual cell-type effects"
  }
  cat(sprintf("%-40s %6.3f %6.3f %s\n", mname, r2, delta, interp))
  prev_r2 <- r2
}

# Residual = unexplained by cell composition
residual_r2 <- 1 - summary(m4)$r.squared
cat(sprintf("%-40s %6s %6.3f %s\n", "Unexplained (residual)", "", residual_r2,
            "Cell-intrinsic + technical noise"))

# ---- 4. Key analysis: SC_Score variance partitioning ----
cat("\n========== SC Score Variance Partitioning ==========\n\n")

cat("This addresses: is SC score decline along MGS explained by purity/composition?\n\n")

s0 <- lm(SC_Score ~ 1, data = df)
s1 <- lm(SC_Score ~ MGS, data = df)
s2 <- lm(SC_Score ~ MGS + Purity, data = df)
s3 <- lm(SC_Score ~ MGS + Purity + Immune + Fibroblasts + Monocytes, data = df)

cat(sprintf("SC_Score ~ MGS:                  R²=%.3f (total MGS association)\n",
            summary(s1)$r.squared))
cat(sprintf("SC_Score ~ MGS + Purity:         R²=%.3f (Δ=%.3f controlling purity)\n",
            summary(s2)$r.squared, summary(s2)$r.squared - summary(s1)$r.squared))
cat(sprintf("SC_Score ~ MGS + Purity + cells: R²=%.3f\n",
            summary(s3)$r.squared))

# Partial R² of MGS after controlling everything
ss_mgs_only <- sum(residuals(lm(SC_Score ~ Purity + Immune + Fibroblasts + Monocytes, data = df))^2)
ss_mgs_full <- sum(residuals(s3)^2)
partial_r2_mgs <- 1 - ss_mgs_full / ss_mgs_only
cat(sprintf("\nPartial R² of MGS (after controlling purity + immune + fibroblasts + monocytes): %.3f\n",
            partial_r2_mgs))
cat(sprintf("→ %.1f%% of SC score variance explained by MGS is INDEPENDENT of cell composition\n",
            100 * partial_r2_mgs))

# ---- 5. Proportion of MGS variance from composition vs cell-intrinsic ----
cat("\n========== MGS Decomposition Summary ==========\n\n")

total_r2 <- summary(m4)$r.squared
comp_r2 <- summary(m1)$r.squared  # purity alone
immune_r2 <- summary(m2)$r.squared - comp_r2  # immune beyond purity
fibro_r2 <- summary(m3)$r.squared - summary(m2)$r.squared
other_r2 <- summary(m4)$r.squared - summary(m3)$r.squared

cat(sprintf("MGS variance explained by:\n"))
cat(sprintf("  Tumor purity (cell fraction):     %.1f%%\n", 100 * comp_r2))
cat(sprintf("  Immune infiltration (beyond):     %.1f%%\n", 100 * immune_r2))
cat(sprintf("  Fibroblasts/stroma:               %.1f%%\n", 100 * fibro_r2))
cat(sprintf("  Other cell types:                 %.1f%%\n", 100 * other_r2))
cat(sprintf("  ──────────────────────────────\n"))
cat(sprintf("  TOTAL explained by composition:   %.1f%%\n", 100 * total_r2))
cat(sprintf("  Unexplained (cell-intrinsic + ε): %.1f%%\n", 100 * residual_r2))

# ---- 6. Per-subtype purity-adjusted SC score ----
cat("\n========== Purity-Adjusted SC Score by Subtype ==========\n\n")
cat("SC score after regressing out purity effect:\n")
sc_resid <- residuals(lm(SC_Score ~ Purity, data = df))
names(sc_resid) <- df$Sample
df$SC_purity_adjusted <- sc_resid

adj_summary <- df %>% group_by(Subtype) %>%
  summarise(
    SC_raw = mean(SC_Score),
    SC_purity_adj = mean(SC_purity_adjusted),
    Purity_mean = mean(Purity)
  )
print(adj_summary)

kw_raw <- kruskal.test(SC_Score ~ Subtype, data = df)
kw_adj <- kruskal.test(SC_purity_adjusted ~ Subtype, data = df)
cat(sprintf("\nSC Score (raw):       KW P = %.4f\n", kw_raw$p.value))
cat(sprintf("SC Score (purity-adj): KW P = %.4f\n", kw_adj$p.value))

if (kw_adj$p.value < 0.05) {
  cat("→ SC score differences persist after purity adjustment\n")
  cat("  Consistent with PARTIAL cell-intrinsic Schwann cell state changes\n")
} else {
  cat("→ SC score differences largely explained by purity\n")
  cat("  Consistent with PREDOMINANTLY compositional effects\n")
}

# ---- 7. Correlation matrix of key variables ----
cat("\n========== Key Variable Correlation Matrix ==========\n\n")
vars <- c("MGS","SC_Score","Purity","Immune","Fibroblasts","Monocytes","Tcells","Bcells")
cor_m <- cor(df[, vars], method = "spearman")

cat("Spearman correlations:\n")
cat(sprintf("%-15s", ""))
for (v in vars) cat(sprintf("%8s", substr(v,1,8)))
cat("\n")
for (v1 in vars) {
  cat(sprintf("%-15s", v1))
  for (v2 in vars) {
    cat(sprintf("%8.3f", cor_m[v1, v2]))
  }
  cat("\n")
}

# ---- Save ----
saveRDS(list(
  variance_partitioning = list(
    mgs_r2_decomposition = c(purity = comp_r2, immune = immune_r2, fibro = fibro_r2,
                             other = other_r2, residual = residual_r2),
    sc_partial_r2_mgs = partial_r2_mgs,
    sc_kw_raw = kw_raw$p.value,
    sc_kw_adj = kw_adj$p.value
  ),
  df = df
), file.path(work_dir, "variance_partitioning_results.rds"))

cat("\n========== Complete ==========\n")
