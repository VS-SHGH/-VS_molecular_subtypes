# =============================================================================
# Revision Analysis 6: MGS External Validation (using TPM.xlsx gene names)
# =============================================================================

library(dplyr)
library(readxl)
library(princurve)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"

# ---- Load discovery TPM ----
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

tpm_disc <- as.matrix(tpm_raw[, our_samples])
rownames(tpm_disc) <- genes
log_tpm_disc <- log2(tpm_disc + 1)
cat(sprintf("Discovery: %d genes x %d samples\n", nrow(log_tpm_disc), ncol(log_tpm_disc)))

# ---- Load external cohorts ----
load_external <- function(path) {
  dat <- read.csv(path, check.names = FALSE)
  genes <- as.character(dat[, 1])
  keep <- genes != "" & !is.na(genes)
  dat <- dat[keep, ]; genes <- genes[keep]
  dup <- duplicated(genes)
  dat <- dat[!dup, ]; genes <- genes[!dup]
  mat <- as.matrix(dat[, -1, drop = FALSE])
  rownames(mat) <- genes
  return(mat)
}

gse141801 <- load_external("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE141801_Gugel/GSE141801_Expression_Log2.csv")
gse39645 <- load_external("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE39645_Torres/GSE39645_Expression_Log2.csv")
cat(sprintf("GSE141801: %d genes x %d samples\n", nrow(gse141801), ncol(gse141801)))
cat(sprintf("GSE39645: %d genes x %d samples\n", nrow(gse39645), ncol(gse39645)))

# ---- 1. Build MGS on discovery ----
mads <- apply(log_tpm_disc, 1, mad)
top_3000 <- names(sort(mads, decreasing = TRUE))[1:3000]
cat(sprintf("Top 3000 MAD genes selected\n"))

pca_d <- prcomp(t(log_tpm_disc[top_3000, ]), center = TRUE, scale. = TRUE)
fit_d <- principal_curve(as.matrix(pca_d$x[, 1:5]), stretch = 0)
mgs_d <- fit_d$lambda
mgs_d <- (mgs_d - min(mgs_d)) / (max(mgs_d) - min(mgs_d))

# Var explained
var_d <- summary(pca_d)$importance[2, 1:10]
cat(sprintf("Discovery PC1=%.1f%%, PC2=%.1f%%, PC1-5=%.1f%%\n",
            100*var_d[1], 100*var_d[2], 100*sum(var_d[1:5])))

# ---- 2. MGS by subtype ----
cl <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]
sample_ids_clean <- gsub("^tpm\\.", "", our_samples)
subtypes <- setNames(cl$Subtype, cl$SampleID)[sample_ids_clean]

mgs_subtype <- data.frame(Sample = our_samples, MGS = mgs_d, Subtype = subtypes)
cat("\nMGS by subtype (discovery):\n")
print(mgs_subtype %>% group_by(Subtype) %>% summarise(mean=mean(MGS), sd=sd(MGS)))

kw_mgs <- kruskal.test(MGS ~ Subtype, data = mgs_subtype)
cat(sprintf("Kruskal-Wallis P = %.6f\n", kw_mgs$p.value))

# ---- 3. External projection function ----
project_mgs <- function(ext_mat, disc_mat, top_genes, disc_pca, disc_fit, disc_mgs) {
  common <- intersect(rownames(disc_mat), rownames(ext_mat))
  top_avail <- intersect(top_genes, common)
  cat(sprintf("  Common genes: %d, top avail: %d\n", length(common), length(top_avail)))

  disc_sub <- disc_mat[top_avail, ]
  ext_sub <- ext_mat[top_avail, ]

  # Center and scale using discovery parameters
  gmeans <- rowMeans(disc_sub)
  gsds <- apply(disc_sub, 1, sd)
  gsds[gsds < 1e-10] <- 1

  if (any(is.na(gmeans)) || any(is.na(gsds))) {
    cat("  WARNING: NA in scaling params\n")
    return(rep(NA, ncol(ext_mat)))
  }

  # Build PCA on discovery subset
  disc_scaled <- (disc_sub - gmeans) / gsds
  pca_m <- prcomp(t(disc_scaled), center = FALSE, scale. = FALSE)

  # Project external
  ext_scaled <- (ext_sub - gmeans) / gsds
  ext_pc <- as.matrix(t(ext_scaled)) %*% pca_m$rotation[, 1:min(5, ncol(pca_m$rotation))]

  # Build new principal curve on discovery PCs
  fit_m <- principal_curve(as.matrix(pca_m$x[, 1:min(5, ncol(pca_m$rotation))]), stretch = 0)
  mgs_m <- fit_m$lambda
  mgs_m <- (mgs_m - min(mgs_m)) / (max(mgs_m) - min(mgs_m))

  # Find nearest point on curve for each external sample
  curve_pts <- fit_m$s[fit_m$ord, 1:min(5, ncol(ext_pc))]
  ext_mgs <- numeric(nrow(ext_pc))
  for (i in seq_len(nrow(ext_pc))) {
    dists <- sqrt(colSums((t(curve_pts) - ext_pc[i, ])^2))
    ext_mgs[i] <- mgs_m[fit_m$ord[which.min(dists)]]
  }
  return(ext_mgs)
}

# ---- 4. Project to GSE141801 ----
cat("\n--- GSE141801 MGS Projection ---\n")
mgs_141801 <- project_mgs(gse141801, log_tpm_disc, top_3000, pca_d, fit_d, mgs_d)
cat(sprintf("GSE141801 MGS: n=%d, range=[%.3f,%.3f], mean=%.3f, sd=%.3f\n",
            sum(!is.na(mgs_141801)), min(mgs_141801,na.rm=TRUE), max(mgs_141801,na.rm=TRUE),
            mean(mgs_141801,na.rm=TRUE), sd(mgs_141801,na.rm=TRUE)))

# ---- 5. Project to GSE39645 ----
cat("\n--- GSE39645 MGS Projection ---\n")
mgs_39645 <- project_mgs(gse39645, log_tpm_disc, top_3000, pca_d, fit_d, mgs_d)
cat(sprintf("GSE39645 MGS: n=%d, range=[%.3f,%.3f], mean=%.3f, sd=%.3f\n",
            sum(!is.na(mgs_39645)), min(mgs_39645,na.rm=TRUE), max(mgs_39645,na.rm=TRUE),
            mean(mgs_39645,na.rm=TRUE), sd(mgs_39645,na.rm=TRUE)))

# ---- 6. Load external clinical data for MGS validation ----
clin_141801 <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE141801_Gugel/Clinical_Data.csv")
clin_39645 <- read.csv("/Users/yanchen/Desktop/VSbulk+单细胞/Bulk VS 110/GSE39645_Torres/Clinical_Data.csv")

# GSE141801: MGS by Barrett Class
names(mgs_141801) <- colnames(gse141801)
clin_141801$MGS <- mgs_141801[clin_141801$Sample_geo_accession]
cat("\nMGS by Barrett Class (GSE141801):\n")
print(clin_141801 %>% filter(!is.na(MGS)) %>% group_by(Class) %>%
  summarise(mean_MGS=mean(MGS), sd_MGS=sd(MGS), n=n()))
kw_141801 <- kruskal.test(MGS ~ Class, data = clin_141801[!is.na(clin_141801$MGS), ])
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw_141801$p.value))

# GSE39645: MGS by Barrett Class
names(mgs_39645) <- colnames(gse39645)
clin_39645$MGS <- mgs_39645[clin_39645$Sample_geo_accession]
cat("\nMGS by Barrett Class (GSE39645):\n")
print(clin_39645 %>% filter(!is.na(MGS)) %>% group_by(Class) %>%
  summarise(mean_MGS=mean(MGS), sd_MGS=sd(MGS), n=n()))
kw_39645 <- kruskal.test(MGS ~ Class, data = clin_39645[!is.na(clin_39645$MGS), ])
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw_39645$p.value))

# ---- 7. MGS vs NF2 status ----
cat("\n--- MGS by NF2 Status ---\n")
clin_141801$NF2_bin <- ifelse(clin_141801$NF2 == "Y", "Mutant", "WT")
nf2_mgs <- clin_141801[!is.na(clin_141801$MGS) & clin_141801$NF2 %in% c("Y", "N"), ]
cat(sprintf("GSE141801: Mutant MGS=%.3f, WT MGS=%.3f\n",
            mean(nf2_mgs$MGS[nf2_mgs$NF2=="Y"]), mean(nf2_mgs$MGS[nf2_mgs$NF2=="N"])))
wt <- wilcox.test(MGS ~ NF2, data = nf2_mgs)
cat(sprintf("Wilcoxon P = %.4f\n", wt$p.value))

clin_39645$NF2_bin <- ifelse(clin_39645$NF2 == "Y", "Mutant", "WT")
nf2_mgs_39 <- clin_39645[!is.na(clin_39645$MGS) & clin_39645$NF2 %in% c("Y", "N"), ]
if (nrow(nf2_mgs_39) > 5) {
  cat(sprintf("GSE39645: Mutant MGS=%.3f, WT MGS=%.3f\n",
              mean(nf2_mgs_39$MGS[nf2_mgs_39$NF2=="Y"]),
              mean(nf2_mgs_39$MGS[nf2_mgs_39$NF2=="N"])))
}

# ---- 8. MGS vs Tumor Size ----
cat("\n--- MGS by Tumor Size (GSE141801) ---\n")
mgs_size <- clin_141801[!is.na(clin_141801$MGS) & !is.na(clin_141801$Size), ]
cat("MGS by Size:\n")
print(mgs_size %>% group_by(Size) %>% summarise(mean_MGS=mean(MGS), sd=sd(MGS), n=n()))
kw_size <- kruskal.test(MGS ~ Size, data = mgs_size)
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw_size$p.value))

# ---- 9. Comparison of MGS distributions across cohorts ----
cat("\n========== Cross-Cohort MGS Comparison ==========\n")
cat(sprintf("Discovery MGS:   mean=%.3f, sd=%.3f, IQR=[%.3f, %.3f]\n",
            mean(mgs_d), sd(mgs_d), quantile(mgs_d,0.25), quantile(mgs_d,0.75)))
cat(sprintf("GSE141801 MGS:   mean=%.3f, sd=%.3f, IQR=[%.3f, %.3f]\n",
            mean(mgs_141801,na.rm=TRUE), sd(mgs_141801,na.rm=TRUE),
            quantile(mgs_141801,0.25,na.rm=TRUE), quantile(mgs_141801,0.75,na.rm=TRUE)))
cat(sprintf("GSE39645 MGS:    mean=%.3f, sd=%.3f, IQR=[%.3f, %.3f]\n",
            mean(mgs_39645,na.rm=TRUE), sd(mgs_39645,na.rm=TRUE),
            quantile(mgs_39645,0.25,na.rm=TRUE), quantile(mgs_39645,0.75,na.rm=TRUE)))

# ---- 10. Spearman correlation: MGS vs immune/stromal scores ----
cat("\n========== MGS vs Microenvironment ==========\n")
# Load ESTIMATE results
est_df <- readRDS(file.path(work_dir, "estimate_results.rds"))
est_df$MGS <- mgs_d

cat("Spearman correlations with MGS:\n")
cat(sprintf("  Stromal score: rho=%.3f, P=%.4f\n",
            cor(est_df$MGS, est_df$Stromal, method="spearman"),
            cor.test(est_df$MGS, est_df$Stromal, method="spearman")$p.value))
cat(sprintf("  Immune score:  rho=%.3f, P=%.4f\n",
            cor(est_df$MGS, est_df$Immune, method="spearman"),
            cor.test(est_df$MGS, est_df$Immune, method="spearman")$p.value))
cat(sprintf("  Purity:        rho=%.3f, P=%.4f\n",
            cor(est_df$MGS, est_df$Purity, method="spearman"),
            cor.test(est_df$MGS, est_df$Purity, method="spearman")$p.value))

# ---- Save ----
saveRDS(list(
  mgs_discovery = mgs_d,
  mgs_gse141801 = mgs_141801,
  mgs_gse39645 = mgs_39645,
  pca_var = var_d,
  mgs_subtype = mgs_subtype
), file.path(work_dir, "mgs_full_results.rds"))

cat("\n========== Complete ==========\n")
