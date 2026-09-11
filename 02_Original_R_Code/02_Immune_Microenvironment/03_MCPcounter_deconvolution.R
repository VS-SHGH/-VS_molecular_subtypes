# =============================================================================
# Revision Analysis 5: MCP-counter Deconvolution + ESTIMATE Purity
# Compares with original CIBERSORT results and addresses reviewer concern
# =============================================================================

library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)

work_dir <- "/Users/yanchen/Desktop/VS_revision/output"

# ---- Load TPM data with gene names ----
tpm_raw <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/TPM.xlsx")
genes <- tpm_raw$gene_name
# Remove duplicates
keep <- !duplicated(genes)
tpm_raw <- tpm_raw[keep, ]
genes <- genes[keep]

# Keep only our 38 discovery samples
our_samples <- c('tpm.P1','tpm.P10','tpm.P11','tpm.P12','tpm.P13','tpm.P14','tpm.P15',
                 'tpm.P16','tpm.P17','tpm.P18','tpm.P19','tpm.P2','tpm.P20','tpm.P21',
                 'tpm.P22','tpm.P23','tpm.P24','tpm.P25','tpm.P26','tpm.P27','tpm.P28',
                 'tpm.P3','tpm.P30','tpm.P32','tpm.P34','tpm.P35','tpm.P36','tpm.P37',
                 'tpm.P38','tpm.P39','tpm.P4','tpm.P40','tpm.P41','tpm.P5','tpm.P6',
                 'tpm.P7','tpm.S5','tpm.S8')

tpm_mat <- as.matrix(tpm_raw[, our_samples])
rownames(tpm_mat) <- genes
cat(sprintf("TPM matrix: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

# ---- 1. ESTIMATE Tumor Purity Analysis ----
cat("\n========== ESTIMATE Analysis ==========\n")

# ESTIMATE requires specific approach for custom data
# We'll compute stromal and immune scores using ESTIMATE gene signatures
# Common genes for ESTIMATE
estimate_stromal <- c("ABI3BP","ADAM12","ADAM33","ADAMTS2","ADAMTS5","ANGPT1","ANGPTL1",
                       "ANTXR1","ASPN","BGN","C11orf96","CCDC80","CD248","CDH11","CLEC11A",
                       "COL10A1","COL11A1","COL12A1","COL1A1","COL1A2","COL3A1","COL4A1",
                       "COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL6A3","COL8A1",
                       "COMP","CRISPLD2","CSF1","CYR61","DCN","DPT","EFEMP2","ELN",
                       "FAP","FBN1","FBLN1","FBLN2","FBLN5","FNDC1","FN1","GLT8D2",
                       "HTRA1","HTRA3","ISLR","LAMA4","LAMB1","LAMC1","LOX","LOXL1",
                       "LOXL2","LRRC15","LUM","MATN2","MFAP5","MGP","MMP14","MXRA5",
                       "MXRA8","NBL1","NTM","OLFML2B","OMD","PCOLCE","PDGFRB","PDLIM3",
                       "PLAT","PRRX1","PTN","RCN3","RIN2","SFRP4","SPARC","SPOCK1",
                       "SRPX2","SULF1","TAGLN","TGFB1","THBS2","THY1","TIMP3","TNFSF4",
                       "VCAN","WISP1","ZFHX4")

estimate_immune <- c("APC","ARHGAP25","BTK","C1QA","C1QB","C2","C3AR","C5AR1","CCL13",
                      "CCL19","CCL22","CCL5","CCR2","CCR5","CD14","CD163","CD1C","CD2",
                      "CD209","CD28","CD3D","CD3E","CD4","CD40","CD40LG","CD48","CD5",
                      "CD52","CD53","CD6","CD68","CD7","CD74","CD79A","CD79B","CD8A",
                      "CD8B","CD86","CD96","CIITA","CLEC4A","CORO1A","CRTAM","CSF2RB",
                      "CTLA4","CTSS","CX3CR1","CXCL10","CXCL13","CXCL9","CYBB","DOCK2",
                      "EVI2B","FASLG","FLT3","FPR3","FYN","GIMAP5","GPR18","GZMA",
                      "HAVCR2","HCLS1","HLA-DMB","HLA-E","HLA-F","HMMR","ICOS","IFNG",
                      "IKZF1","IL10RA","IL2RA","IL2RG","IL7R","IRF8","ITGAL","ITK",
                      "KLRB1","KLRK1","LAG3","LAIR1","LCK","LCP2","LY86","LY96",
                      "MARCO","MS4A1","NCKAP1L","NKG7","P2RY13","PDCD1","PDCD1LG2",
                      "PRF1","PSMB9","PTPRC","PYHIN1","RAC2","SELPLG","SLAMF6","SLAMF8",
                      "SPN","ST8SIA4","TARP","TCL1A","TLR10","TLR7","TLR8","TNFRSF13B",
                      "TNFRSF17","TRAT1","TRBC1","TRBC2","TRBV28","TREM2","UBA7","WDFY4")

# Compute ssGSEA-like scores (mean expression of signature genes)
compute_signature <- function(expr_mat, sig_genes) {
  avail <- intersect(sig_genes, rownames(expr_mat))
  if (length(avail) < 5) return(rep(NA, ncol(expr_mat)))
  colMeans(expr_mat[avail, , drop = FALSE])
}

stromal_score <- compute_signature(tpm_mat, estimate_stromal)
immune_score <- compute_signature(tpm_mat, estimate_immune)
purity_estimate <- 1 - (stromal_score + immune_score) / (max(stromal_score + immune_score, na.rm = TRUE))

# Load subtype info
cl <- read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]
sample_to_subtype <- setNames(cl$Subtype, cl$SampleID)

# Map scores to subtypes (strip "tpm." prefix)
sample_ids_clean <- gsub("^tpm\\.", "", our_samples)
subtypes <- sample_to_subtype[sample_ids_clean]
names(subtypes) <- our_samples
est_df <- data.frame(
  Sample = our_samples,
  Subtype = subtypes,
  Stromal = stromal_score,
  Immune = immune_score,
  Purity = purity_estimate
)

cat("\nESTIMATE scores by subtype:\n")
print(est_df %>% group_by(Subtype) %>%
  summarise(Stromal=mean(Stromal,na.rm=TRUE), Immune=mean(Immune,na.rm=TRUE),
            Purity=mean(Purity,na.rm=TRUE)))

# Kruskal-Wallis test
for (var in c("Stromal", "Immune", "Purity")) {
  kw <- kruskal.test(as.formula(paste(var, "~ Subtype")), data = est_df)
  cat(sprintf("%s: Kruskal-Wallis P = %.4f\n", var, kw$p.value))
}

# ---- 2. MCP-counter Deconvolution ----
cat("\n========== MCP-counter Deconvolution ==========\n")

# Try loading MCPcounter
mcp_available <- requireNamespace("MCPcounter", quietly = TRUE)
cat(sprintf("MCPcounter available: %s\n", mcp_available))

if (!mcp_available) {
  cat("Installing MCPcounter...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source", upgrade = "never")
  library(MCPcounter)
} else {
  library(MCPcounter)
}

# MCPcounter needs non-log-transformed TPM (or counts)
# Our TPM data is already in linear scale
cat("Running MCPcounter...\n")
mcp_results <- tryCatch({
  MCPcounter.estimate(tpm_mat, featuresType = "HUGO_symbols")
}, error = function(e) {
  cat(sprintf("MCPcounter error: %s\n", e$message))
  cat("Trying with gene symbols...\n")
  # Try alternative approach
  MCPcounter::MCPcounter.estimate(tpm_mat, featuresType = "HUGO_symbols")
})

if (!is.null(mcp_results) && !inherits(mcp_results, "try-error")) {
  cat(sprintf("MCPcounter: %d cell types x %d samples\n", nrow(mcp_results), ncol(mcp_results)))
  cat("Cell types:\n")
  print(rownames(mcp_results))

  # Compare MCPcounter cell types across subtypes
  for (ct in rownames(mcp_results)) {
    ct_vals <- as.numeric(mcp_results[ct, ])
    names(ct_vals) <- colnames(mcp_results)
    ct_df <- data.frame(Sample = names(ct_vals), Value = ct_vals,
                        Subtype = subtypes[names(ct_vals)])
    kw <- kruskal.test(Value ~ Subtype, data = ct_df)
    if (kw$p.value < 0.1) {
      cat(sprintf("%s: KW P = %.4f", ct, kw$p.value))
      means <- ct_df %>% group_by(Subtype) %>% summarise(m = mean(Value, na.rm = TRUE))
      cat(sprintf(" | C1=%.2f C2=%.2f C3=%.2f\n", means$m[1], means$m[2], means$m[3]))
    }
  }

  # Save MCPcounter results
  saveRDS(mcp_results, file.path(work_dir, "mcp_counter_results.rds"))
} else {
  cat("MCPcounter failed. Will document this limitation.\n")
}

# ---- Save ESTIMATE results ----
saveRDS(est_df, file.path(work_dir, "estimate_results.rds"))

cat("\n========== Complete ==========\n")
