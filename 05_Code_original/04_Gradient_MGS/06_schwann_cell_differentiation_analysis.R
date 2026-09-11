# =============================================================================
# Revision Analysis 7: Schwann Cell Differentiation vs Composition Analysis
# Addresses reviewer concern: "de-differentiation" vs compositional shift
# =============================================================================

library(dplyr)
library(readxl)
library(ggplot2)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"

# ---- Load data ----
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

# ---- 1. Schwann Cell Differentiation Gene Scores ----
cat("========== Schwann Cell Differentiation Markers ==========\n")

# Known Schwann cell differentiation markers (from literature: Jessen & Mirsky, Barrett et al.)
sc_diff_genes <- c("SOX10", "PLP1", "EGR2", "MPZ", "GAP43", "NGFR", "MBP", "PMP22",
                    "PRX", "MAL", "CDH19", "L1CAM", "GFAP", "S100B")

sc_avail <- intersect(sc_diff_genes, rownames(tpm_mat))
cat(sprintf("Schwann cell genes available: %d/%d\n", length(sc_avail), length(sc_diff_genes)))
cat(sprintf("Missing: %s\n", paste(setdiff(sc_diff_genes, sc_avail), collapse=", ")))

# Compute per-gene expression by subtype
cat("\nSchwann cell marker expression by subtype:\n")
for (g in sc_avail) {
  vals <- log_tpm[g, ]
  c1_mean <- mean(vals[subtypes == "C1"])
  c2_mean <- mean(vals[subtypes == "C2"])
  c3_mean <- mean(vals[subtypes == "C3"])
  kw <- kruskal.test(vals ~ subtypes)
  cat(sprintf("  %s: C1=%.2f, C2=%.2f, C3=%.2f, P=%.4f\n",
              g, c1_mean, c2_mean, c3_mean, kw$p.value))
}

# Schwann cell differentiation score (mean of all markers)
sc_score <- colMeans(log_tpm[sc_avail, ])
sc_by_subtype <- data.frame(Sample = names(sc_score), SC_Score = sc_score, Subtype = subtypes)
cat("\nSchwann Cell Differentiation Score by Subtype:\n")
print(sc_by_subtype %>% group_by(Subtype) %>% summarise(mean=mean(SC_Score), sd=sd(SC_Score)))
kw_sc <- kruskal.test(SC_Score ~ Subtype, data = sc_by_subtype)
cat(sprintf("Kruskal-Wallis P = %.6f\n", kw_sc$p.value))

# ---- 2. Injury-like Schwann Cell Genes ----
cat("\n========== Injury-like Schwann Cell Markers ==========\n")

injury_genes <- c("NGFR", "RUNX2", "SPP1", "FN1", "CD44", "VIM", "NES",
                   "SOX2", "CDH2", "SNAI2", "TWIST1", "ZEB1", "CTNNB1")

injury_avail <- intersect(injury_genes, rownames(tpm_mat))
cat(sprintf("Injury genes available: %d/%d\n", length(injury_avail), length(injury_genes)))

injury_score <- colMeans(log_tpm[injury_avail, ])
injury_by_subtype <- data.frame(Sample = names(injury_score), Injury_Score = injury_score, Subtype = subtypes)
cat("\nInjury-like Score by Subtype:\n")
print(injury_by_subtype %>% group_by(Subtype) %>% summarise(mean=mean(Injury_Score), sd=sd(Injury_Score)))
kw_injury <- kruskal.test(Injury_Score ~ Subtype, data = injury_by_subtype)
cat(sprintf("Kruskal-Wallis P = %.6f\n", kw_injury$p.value))

# ---- 3. Compositional vs Cell-Intrinsic Analysis ----
cat("\n========== Compositional vs Cell-Intrinsic Analysis ==========\n")

# Load ESTIMATE and MCPcounter
est_df <- readRDS(file.path(work_dir, "estimate_results.rds"))
mcp <- readRDS(file.path(work_dir, "mcp_counter_results.rds"))
mgs_full <- readRDS(file.path(work_dir, "mgs_full_results.rds"))
mgs_d <- mgs_full$mgs_discovery

# Put everything together
comp_df <- data.frame(
  Sample = our_samples,
  Subtype = subtypes,
  MGS = mgs_d,
  SC_Score = sc_score[our_samples],
  Injury_Score = injury_score[our_samples],
  Purity = est_df$Purity,
  Immune_Score = est_df$Immune,
  Stromal_Score = est_df$Stromal,
  Tcells_MCP = as.numeric(mcp["T cells", our_samples]),
  Bcells_MCP = as.numeric(mcp["B lineage", our_samples]),
  Monocytes_MCP = as.numeric(mcp["Monocytic lineage", our_samples]),
  Fibroblasts_MCP = as.numeric(mcp["Fibroblasts", our_samples]),
  NK_MCP = as.numeric(mcp["NK cells", our_samples])
)

# Key question: Is SC score decline explained by purity alone?
cat("\n--- Partial correlation: SC_Score vs MGS, controlling for Purity ---\n")
# Simple approach: linear models
m1 <- lm(SC_Score ~ MGS, data = comp_df)
m2 <- lm(SC_Score ~ MGS + Purity, data = comp_df)
cat(sprintf("SC_Score ~ MGS: R2=%.3f, MGS coef=%.3f, P=%.4f\n",
            summary(m1)$r.squared, coef(m1)[2], summary(m1)$coefficients[2,4]))
cat(sprintf("SC_Score ~ MGS + Purity: R2=%.3f, MGS coef=%.3f, Purity coef=%.3f\n",
            summary(m2)$r.squared, coef(m2)[2], coef(m2)[3]))
cat(sprintf("Partial R2 of MGS: %.3f\n",
            summary(m2)$r.squared - summary(lm(SC_Score ~ Purity, data=comp_df))$r.squared))

# ---- 4. MKI67 expression validation ----
cat("\n========== MKI67 Expression ==========\n")
if ("MKI67" %in% rownames(tpm_mat)) {
  mki67_vals <- log_tpm["MKI67", ]
  mki67_df <- data.frame(Sample = names(mki67_vals), MKI67 = mki67_vals, Subtype = subtypes)
  cat("MKI67 by subtype:\n")
  print(mki67_df %>% group_by(Subtype) %>% summarise(mean=mean(MKI67), sd=sd(MKI67)))
  kw_mki67 <- kruskal.test(MKI67 ~ Subtype, data = mki67_df)
  cat(sprintf("Kruskal-Wallis P = %.4f\n", kw_mki67$p.value))
}

# ---- 5. Proliferation gene set score ----
cat("\n========== Proliferation Score ==========\n")
prolif_genes <- c("MKI67", "PCNA", "TOP2A", "CCNB1", "CCNB2", "CDK1", "BIRC5",
                   "AURKA", "AURKB", "CDC20", "BUB1", "PLK1", "MCM2", "MCM3", "MCM4")

prolif_avail <- intersect(prolif_genes, rownames(tpm_mat))
prolif_score <- colMeans(log_tpm[prolif_avail, ])
prolif_df <- data.frame(Sample = names(prolif_score), Prolif = prolif_score, Subtype = subtypes)
cat("Proliferation score by subtype:\n")
print(prolif_df %>% group_by(Subtype) %>% summarise(mean=mean(Prolif), sd=sd(Prolif)))
kw_prolif <- kruskal.test(Prolif ~ Subtype, data = prolif_df)
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw_prolif$p.value))

# ---- 6. ECM/Fibrosis gene score ----
cat("\n========== ECM/Fibrosis Score ==========\n")
ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", "COL5A2", "COL6A1",
               "COL6A2", "COL6A3", "FN1", "LOX", "LOXL1", "LOXL2", "MMP2", "MMP9",
               "MMP14", "TIMP1", "TIMP3", "TGFB1", "CTGF", "POSTN", "SPARC", "THBS2")

ecm_avail <- intersect(ecm_genes, rownames(tpm_mat))
ecm_score <- colMeans(log_tpm[ecm_avail, ])
ecm_df <- data.frame(Sample = names(ecm_score), ECM = ecm_score, Subtype = subtypes)
cat("ECM score by subtype:\n")
print(ecm_df %>% group_by(Subtype) %>% summarise(mean=mean(ECM), sd=sd(ECM)))
kw_ecm <- kruskal.test(ECM ~ Subtype, data = ecm_df)
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw_ecm$p.value))

# ---- 7. Immune checkpoint genes ----
cat("\n========== Immune Checkpoint Expression ==========\n")
checkpoint_genes <- c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "CD274", "PDCD1LG2")
check_avail <- intersect(checkpoint_genes, rownames(tpm_mat))
cat(sprintf("Checkpoint genes available: %d/%d\n", length(check_avail), length(checkpoint_genes)))

for (g in check_avail) {
  vals <- log_tpm[g, ]
  c1m <- mean(vals[subtypes == "C1"])
  c2m <- mean(vals[subtypes == "C2"])
  c3m <- mean(vals[subtypes == "C3"])
  kw <- kruskal.test(vals ~ subtypes)
  cat(sprintf("  %s: C1=%.2f C2=%.2f C3=%.2f P=%.4f\n", g, c1m, c2m, c3m, kw$p.value))
}

# ---- 8. Correlation matrix: key variables ----
cat("\n========== Key Correlations ==========\n")
cor_vars <- comp_df[, c("MGS", "SC_Score", "Injury_Score", "Purity", "Immune_Score",
                         "Tcells_MCP", "Bcells_MCP", "Monocytes_MCP", "Fibroblasts_MCP")]
cor_mat <- cor(cor_vars, method = "spearman", use = "complete.obs")

cat("\nSpearman correlations with MGS:\n")
for (v in setdiff(colnames(cor_mat), "MGS")) {
  cat(sprintf("  %s: rho=%.3f\n", v, cor_mat["MGS", v]))
}

# ---- Save ----
saveRDS(comp_df, file.path(work_dir, "composition_analysis.rds"))
cat("\n========== Complete ==========\n")
