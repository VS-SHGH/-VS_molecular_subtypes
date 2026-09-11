# =============================================================================
# VS (Vestibular Schwannoma) Immune Microenvironment Analysis
# Description: Optimized immune scoring, subtype-specific functional scoring,
#              and publication-quality visualization for VS C1/C2/C3
#
# Refactored Version - Key Improvements:
#   [R1] Removed duplicate SECTION 30 code blocks
#   [R2] Unified all annotations to English
#   [R3] Defined magic numbers as constants
#   [R4] Consistent variable naming
#   [R5] Removed residual debug code (score_rows in SECTION 13)
#   [R6] Consolidated ssGSEA calls for performance
#   [R7] Improved memory management
#   [R8] PDF-only output (no PNG)
#   [R9] Fixed CYT score formula (geometric mean)
#   [R10] Enhanced statistical testing with sample size threshold
# =============================================================================

# =============================================================================
# SECTION 0: ENVIRONMENT SETUP & PACKAGE LOADING
# =============================================================================

# --- Constants ---
ESTIMATE_INTERCEPT <- 0.6049872018
ESTIMATE_SLOPE     <- 0.0001467884
MIN_GENES_SSGSEA   <- 3
MIN_SAMPLES_PARAM  <- 20
EPSILON            <- 1e-9

# --- Package Installation ---
if (!requireNamespace("estimate", quietly = TRUE)) {
  message("Installing estimate from R-Forge...")
  install.packages("estimate", repos = "http://r-forge.r-project.org",
                   dependencies = TRUE)
}

if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  message("Installing MCPcounter from GitHub...")
  remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")
}

if (!requireNamespace("GSVA", quietly = TRUE)) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GSVA")
}

GSVA_NEW_API <- tryCatch({
  packageVersion("GSVA") >= "1.44.0"
}, error = function(e) FALSE)
message(sprintf("GSVA version: %s | Using %s API",
                as.character(packageVersion("GSVA")),
                ifelse(GSVA_NEW_API, "new (>=1.44)", "old (<1.44)")))

suppressPackageStartupMessages({
  library(estimate)
  library(MCPcounter)
  library(GSVA)
  library(ggplot2)
  library(ggdist)
  library(ggpubr)
  library(ggsignif)
  library(ComplexHeatmap)
  library(circlize)
  library(fmsb)
  library(corrplot)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(FSA)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
  library(pheatmap)
})

# Optional packages with graceful fallback
HAS_GGRASTR <- requireNamespace("ggrastr", quietly = TRUE)
if (HAS_GGRASTR) library(ggrastr)

HAS_CAIRO <- requireNamespace("Cairo", quietly = TRUE)
if (HAS_CAIRO) library(Cairo)

# --- Global Paths ---
RDATA_PATH  <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/RData/Step2_Clustering.RData"
OUTPUT_DIR  <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/step8_Immune"
TMPDIR      <- file.path(OUTPUT_DIR, "tmp_estimate")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TMPDIR,     recursive = TRUE, showWarnings = FALSE)

# --- Unified Color Palette ---
# C1 = green (Proliferative), C2 = blue (Mesenchymal/Myelination), C3 = red (Immune)
SUBTYPE_COLORS <- c(
  C1 = "#00A087",
  C2 = "#4DBBD5",
  C3 = "#E64B35"
)
SUBTYPE_LEVELS <- c("C1", "C2", "C3")

# --- Font Configuration ---
BASE_FAMILY <- "Arial"
if (HAS_CAIRO) {
  test_fonts <- tryCatch({
    Cairo::CairoFonts()
    "Arial" %in% names(grDevices::pdfFonts())
  }, error = function(e) FALSE)
  if (!test_fonts) {
    BASE_FAMILY <- "Helvetica"
    message("Arial not found, using Helvetica as fallback")
  }
} else {
  BASE_FAMILY <- "sans"
  message("Cairo not available, using base PDF output")
}

# --- Publication Theme ---
theme_pub <- function(base_size = 11, base_family = BASE_FAMILY) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line        = element_line(linewidth = 0.5, colour = "black"),
      axis.text        = element_text(color = "black", size = base_size - 1),
      axis.title       = element_text(color = "black", size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size, face = "bold"),
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      hjust = 0.5),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size),
      panel.grid       = element_blank()
    )
}
theme_set(theme_pub())

# --- Save Helper: PDF Only ---
save_plot <- function(plot_obj, filename_base, width = 7, height = 5) {
  pdf_path <- file.path(OUTPUT_DIR, paste0(filename_base, ".pdf"))
  
  tryCatch({
    if (HAS_CAIRO) {
      Cairo::CairoPDF(pdf_path, width = width, height = height,
                      family = BASE_FAMILY)
    } else {
      pdf(pdf_path, width = width, height = height)
    }
    print(plot_obj)
    dev.off()
    message("Saved PDF: ", pdf_path)
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("PDF save failed: ", e$message)
  })
}

# Wrapper for base-graphics plots
save_base_plot <- function(expr, filename_base, width = 7, height = 6) {
  pdf_path <- file.path(OUTPUT_DIR, paste0(filename_base, ".pdf"))
  
  tryCatch({
    if (HAS_CAIRO) {
      Cairo::CairoPDF(pdf_path, width = width, height = height,
                      family = BASE_FAMILY)
    } else {
      pdf(pdf_path, width = width, height = height)
    }
    eval(expr)
    dev.off()
    message("Saved PDF: ", pdf_path)
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("PDF save failed: ", e$message)
  })
}

# --- GSVA-version-aware ssGSEA wrapper ---
run_ssgsea <- function(expr_mat, gene_sets, normalize = TRUE) {
  expr_mat <- as.matrix(expr_mat)
  gene_sets <- Filter(function(g) length(g) >= MIN_GENES_SSGSEA, gene_sets)
  if (length(gene_sets) == 0) stop("No gene sets with >= 3 genes found")
  
  if (GSVA_NEW_API) {
    param  <- ssgseaParam(expr_mat, gene_sets, normalize = normalize)
    result <- gsva(param, verbose = FALSE)
  } else {
    result <- gsva(expr_mat, gene_sets, method = "ssgsea",
                   ssgsea.norm = normalize, verbose = FALSE)
  }
  result
}

# --- Statistical Test Helper ---
auto_stat_test <- function(df, score_col, group_col = "Subtype") {
  scores <- df[[score_col]]
  groups <- df[[group_col]]
  test_df <- data.frame(score = scores, group = as.factor(groups))
  test_df <- test_df[!is.na(test_df$score), ]
  
  grp_list <- split(test_df$score, test_df$group)
  n_per_grp <- sapply(grp_list, length)
  
  # Use non-parametric if any group < MIN_SAMPLES_PARAM
  use_nonparam <- any(n_per_grp < MIN_SAMPLES_PARAM)
  
  if (!use_nonparam) {
    normal_flags <- sapply(grp_list, function(x) {
      if (length(x) < 3) return(TRUE)
      shapiro.test(x)$p.value > 0.05
    })
    use_nonparam <- !all(normal_flags)
  }
  
  if (!use_nonparam && length(grp_list) >= 2) {
    method   <- "ANOVA + Tukey HSD"
    fit      <- aov(score ~ group, data = test_df)
    global_p <- summary(fit)[[1]][["Pr(>F)"]][1]
    pw       <- as.data.frame(TukeyHSD(fit)[["group"]])
    pw$comparison <- rownames(pw); rownames(pw) <- NULL
    pw_df <- pw[, c("comparison", "p adj")]
    colnames(pw_df) <- c("comparison", "p.adj")
  } else {
    method   <- "Kruskal-Wallis + Dunn (BH)"
    kt       <- kruskal.test(score ~ group, data = test_df)
    global_p <- kt$p.value
    dt       <- dunnTest(score ~ group, data = test_df, method = "bh")
    pw_df    <- dt$res[, c("Comparison", "P.adj")]
    colnames(pw_df) <- c("comparison", "p.adj")
  }
  
  message(sprintf("[%s] %s | Global p = %.4f", score_col, method, global_p))
  list(method = method, global_p = global_p, pairwise_df = pw_df)
}

# Helper: parse comparison strings
parse_comparison <- function(comp_str) {
  parts <- trimws(strsplit(comp_str, "\\s*-\\s*")[[1]])
  parts[1:2]
}

# Helper: single-gene boxplot with significance
plot_gene_box <- function(gene, expr_mat, meta_df, stat_res = NULL,
                          log_transform = TRUE) {
  if (!gene %in% rownames(expr_mat)) return(NULL)
  vals <- as.numeric(expr_mat[gene, meta_df$SampleID])
  if (log_transform) vals <- log2(vals + 1)
  df_g <- data.frame(Subtype = meta_df$Subtype, Value = vals)
  
  if (is.null(stat_res)) {
    stat_res <- auto_stat_test(df_g, "Value", "Subtype")
  }
  
  pw   <- stat_res$pairwise_df
  comp <- lapply(pw$comparison, parse_comparison)
  
  ggplot(df_g, aes(x = Subtype, y = Value, fill = Subtype)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5,
                 width = 0.55, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.5, color = "black") +
    geom_signif(
      comparisons    = comp,
      map_signif_level = function(p) {
        if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
      },
      step_increase  = 0.10, tip_length = 0.02,
      color = "black", textsize = 3
    ) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = if (log_transform) "log2(TPM+1)" else "Expression",
         title = gene) +
    theme_pub() +
    theme(legend.position = "none",
          plot.title = element_text(face = "italic", size = 10))
}

# =============================================================================
# SECTION 1: DATA LOADING & VALIDATION
# =============================================================================

message("\n=== Loading data ===")
load(RDATA_PATH)

stopifnot(
  "clean_tpm not found" = exists("clean_tpm"),
  "meta_step2 not found" = exists("meta_step2"),
  "SampleID missing"    = "SampleID" %in% colnames(meta_step2),
  "Subtype missing"     = "Subtype"  %in% colnames(meta_step2)
)

meta_step2$Subtype <- factor(meta_step2$Subtype, levels = SUBTYPE_LEVELS)

common_samples <- intersect(colnames(clean_tpm), meta_step2$SampleID)
stopifnot("No overlapping samples" = length(common_samples) > 0)
message(sprintf("Common samples: %d / TPM: %d / Meta: %d",
                length(common_samples), ncol(clean_tpm), nrow(meta_step2)))

clean_tpm  <- clean_tpm[, common_samples]
meta_step2 <- meta_step2[match(common_samples, meta_step2$SampleID), ]

gene_names    <- rownames(clean_tpm)
ensembl_count <- sum(grepl("^ENSG\\d+", gene_names))
if (ensembl_count > 0) {
  warning(sprintf(
    "%d Ensembl IDs in rownames; MCP-counter & ESTIMATE need HGNC symbols.",
    ensembl_count))
}
message(sprintf("Gene matrix: %d genes x %d samples",
                nrow(clean_tpm), ncol(clean_tpm)))

HAS_CLINICAL <- "SizeGrade" %in% colnames(meta_step2) ||
  "TumorSize_cm" %in% colnames(meta_step2)

# =============================================================================
# SECTION 2: ESTIMATE ANALYSIS
# =============================================================================

message("\n=== ESTIMATE Analysis ===")

estimate_input_file  <- file.path(TMPDIR, "tpm_for_estimate.txt")
estimate_output_file <- file.path(TMPDIR, "estimate_scores.gct")
estimate_filtered    <- file.path(TMPDIR, "tpm_filtered_estimate.gct")

# Write TPM matrix
tpm_df_est <- as.data.frame(clean_tpm)
tpm_df_est <- cbind(
  NAME        = rownames(tpm_df_est),
  Description = rownames(tpm_df_est),
  tpm_df_est
)
write.table(tpm_df_est, file = estimate_input_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

# Filter common genes
filterCommonGenes(
  input.f  = estimate_input_file,
  output.f = estimate_filtered,
  id       = "GeneSymbol"
)

# Run ESTIMATE (illumina for TPM data)
estimateScore(
  input.ds  = estimate_filtered,
  output.ds = estimate_output_file,
  platform  = "illumina"
)

# Robust GCT parsing
est_raw <- read.table(
  estimate_output_file,
  header           = TRUE,
  skip             = 2,
  sep              = "\t",
  row.names        = 1,
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

desc_cols  <- grep("^Description", colnames(est_raw), value = TRUE)
est_scores <- est_raw[, !colnames(est_raw) %in% desc_cols, drop = FALSE]

estimate_t <- as.data.frame(t(est_scores))
estimate_t[] <- lapply(estimate_t, function(x) suppressWarnings(as.numeric(x)))
estimate_t$SampleID <- rownames(estimate_t)

# Rename columns
colnames(estimate_t) <- gsub("StromalScore",  "Stromal_Score",  colnames(estimate_t))
colnames(estimate_t) <- gsub("ImmuneScore",   "Immune_Score",   colnames(estimate_t))
colnames(estimate_t) <- gsub("ESTIMATEScore", "ESTIMATE_Score", colnames(estimate_t))

# Calculate Tumor Purity (Yoshihara et al. 2013)
if (!"Tumor_Purity" %in% colnames(estimate_t)) {
  estimate_t$Tumor_Purity <- cos(ESTIMATE_INTERCEPT + ESTIMATE_SLOPE * estimate_t$ESTIMATE_Score)
  estimate_t$Tumor_Purity <- pmin(pmax(estimate_t$Tumor_Purity, 0), 1)
  message("Tumor_Purity calculated. Range: ",
          round(min(estimate_t$Tumor_Purity, na.rm = TRUE), 3), " - ",
          round(max(estimate_t$Tumor_Purity, na.rm = TRUE), 3))
}

# Safe merge
meta_step2 <- meta_step2 %>%
  dplyr::select(-dplyr::any_of(c(
    "Immune_Score", "Stromal_Score", "ESTIMATE_Score", "Tumor_Purity"
  ))) %>%
  dplyr::left_join(estimate_t, by = "SampleID")

estimate_check <- grep("Purity|Immune|Stromal|ESTIMATE", colnames(meta_step2), value = TRUE)
message("ESTIMATE columns added: ", paste(estimate_check, collapse = ", "))
stopifnot("Missing ESTIMATE columns" = length(estimate_check) == 4)

# Cleanup
suppressWarnings(file.remove(estimate_input_file, estimate_filtered, estimate_output_file))
rm(est_raw, tpm_df_est, est_scores, estimate_t)
gc()
message("ESTIMATE complete.")

# =============================================================================
# SECTION 3: CONSOLIDATED ssGSEA SCORING
# =============================================================================

message("\n=== Consolidated ssGSEA Scoring ===")

# Define all gene sets with literature references
# Ref: Proliferation - Tirosh et al. 2016 Science; TIS - Ayers et al. 2017 JCI
# Ref: Checkpoint - Auslander et al. 2018 Nat Med; ECM - Naba et al. 2012 Mol Cell Proteomics

PROLIF_GENES <- c(
  "MKI67","PCNA","TOP2A","CCNB1","CCNB2","CCNA2",
  "CDK1","BUB1","AURKB","PLK1","CENPF","MCM2",
  "MCM6","E2F1","MYBL2","RRM2","TYMS","CDC20",
  "BIRC5","PTTG1"
)

ECM_GENES <- c(
  "COL1A1","COL1A2","COL3A1","COL4A1","COL4A2",
  "COL5A1","COL5A2","COL6A1","COL6A2","COL6A3",
  "FN1","POSTN","LAMA2","LAMA4","LAMB1",
  "LAMC1","TNC","THBS1","THBS2","COMP",
  "SPARC","VCAN","HAPLN1","HSPG2","NID1"
)

MYELIN_GENES <- c(
  "MBP","MPZ","PMP22","PRX","EGR2",
  "EGR1","NCAM1","SOX10","S100B","NGFR",
  "GAP43","CLDN19","KCNA1","MAG","CNP"
)

FIBROSIS_GENES <- c(
  "COL1A1","POSTN","FN1","ACTA2","FAP",
  "PDPN","COL10A1","CTHRC1","LOXL2","PRRX1"
)

TIS_GENES <- c(
  "CD3D","IDO1","CIITA","HLA-DQA1","CD276",
  "LAG3","TIGIT","CD8A","PDCD1LG2","CD274",
  "CXCL9","CXCL10","CMKLR1","NKG7","CCL5",
  "PSMB10","CXCR6","HLA-E"
)

BCELL_GENES <- c(
  "CD19","CD79A","MS4A1","IGHG1","IGHG2",
  "IGHM","CD38","SDC1","JCHAIN","IGKC",
  "IGLC2","BLK","PAX5","BANK1","CR2"
)

CHECKPOINT_GENES <- c(
  "PDCD1","CD274","CTLA4","HAVCR2",
  "LAG3","TIGIT","VSIR","CD96","SIGLEC7"
)

ANGIO_GENES <- c(
  "VEGFA","VEGFB","VEGFC","KDR","FLT1",
  "FLT4","ANGPT1","ANGPT2","TEK","PDGFB",
  "PDGFRB","ENG","CDH5","PECAM1","THBS1",
  "THBS2","COL4A3","NRP1","NOTCH4","HIF1A"
)

HLA1_GENES <- c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP","NLRC5")

CHEMO_PROG_GENES <- c("CXCL9","CXCL10","CCL5","CD86","HAVCR2","LGALS9")

# Build consolidated gene set list
ALL_GENE_SETS <- list(
  Proliferation_Score    = intersect(PROLIF_GENES,     rownames(clean_tpm)),
  ECM_Score              = intersect(ECM_GENES,        rownames(clean_tpm)),
  Myelination_Score      = intersect(MYELIN_GENES,     rownames(clean_tpm)),
  Fibrosis_Score         = intersect(FIBROSIS_GENES,   rownames(clean_tpm)),
  TIS_Score              = intersect(TIS_GENES,        rownames(clean_tpm)),
  BCell_Score            = intersect(BCELL_GENES,      rownames(clean_tpm)),
  Checkpoint_Score       = intersect(CHECKPOINT_GENES, rownames(clean_tpm)),
  Angiogenesis_Score     = intersect(ANGIO_GENES,      rownames(clean_tpm)),
  HLA_ClassI_Score       = intersect(HLA1_GENES,       rownames(clean_tpm)),
  Chemokine_Program_Score = intersect(CHEMO_PROG_GENES, rownames(clean_tpm))
)

# Filter gene sets with sufficient genes
ALL_GENE_SETS <- Filter(function(g) length(g) >= MIN_GENES_SSGSEA, ALL_GENE_SETS)

# Report gene set coverage
for (gs_name in names(ALL_GENE_SETS)) {
  message(sprintf("  %s: %d genes", gs_name, length(ALL_GENE_SETS[[gs_name]])))
}

# Run ssGSEA once for all gene sets
if (length(ALL_GENE_SETS) > 0) {
  ssgsea_results <- run_ssgsea(clean_tpm, ALL_GENE_SETS)
  
  for (score_name in rownames(ssgsea_results)) {
    meta_step2[[score_name]] <- as.numeric(ssgsea_results[score_name, meta_step2$SampleID])
  }
  rm(ssgsea_results)
}

# CYT Score: geometric mean (Rooney et al. 2015 Cell)
CYT_GENES <- c("GZMA", "PRF1")
cyt_present <- intersect(CYT_GENES, rownames(clean_tpm))
if (length(cyt_present) == 2) {
  gzma_vals <- as.numeric(clean_tpm["GZMA", meta_step2$SampleID])
  prf1_vals <- as.numeric(clean_tpm["PRF1", meta_step2$SampleID])
  meta_step2$CYT_Score <- sqrt(gzma_vals * prf1_vals)
  message("CYT Score calculated using geometric mean (Rooney et al. 2015)")
} else if (length(cyt_present) >= 1) {
  meta_step2$CYT_Score <- colMeans(
    log2(clean_tpm[cyt_present, meta_step2$SampleID, drop = FALSE] + 1))
  warning("CYT Score: only one gene available, using log2 mean")
}

gc()
message("ssGSEA scoring complete.")

# =============================================================================
# SECTION 4: MCP-COUNTER DECONVOLUTION
# =============================================================================

message("\n=== MCP-counter ===")

# NOTE: MCP-counter was trained on common solid tumors. Results for
# schwannoma (neural origin) should be validated with marker gene correlations.

mcp_result <- MCPcounter.estimate(
  expression   = clean_tpm,
  featuresType = "HUGO_symbols"
)
mcp_df <- as.data.frame(t(mcp_result))
mcp_df$SampleID <- rownames(mcp_df)
colnames(mcp_df) <- gsub("\\s+", "_", colnames(mcp_df))
colnames(mcp_df) <- gsub("[^A-Za-z0-9_]", "", colnames(mcp_df))

meta_step2 <- left_join(meta_step2, mcp_df, by = "SampleID")
MCP_CELL_TYPES <- colnames(mcp_df)[colnames(mcp_df) != "SampleID"]
message(sprintf("MCP-counter: %d cell types estimated.", length(MCP_CELL_TYPES)))

rm(mcp_result, mcp_df)
gc()

# =============================================================================
# SECTION 5: M1/M2 MACROPHAGE POLARIZATION
# =============================================================================

message("\n=== M1/M2 Polarization ===")

M1_GENES <- c("CD86","IL6","TNF","CXCL10","IL1B","NOS2","CXCL9")
M2_GENES <- c("MRC1","CD163","IL10","ARG1","TGFB1","CCL22","CD200R1")

m1_present <- intersect(M1_GENES, rownames(clean_tpm))
m2_present <- intersect(M2_GENES, rownames(clean_tpm))

if (length(m1_present) >= 2 && length(m2_present) >= 2) {
  m1_score <- colMeans(log2(clean_tpm[m1_present, meta_step2$SampleID, drop=FALSE] + 1))
  m2_score <- colMeans(log2(clean_tpm[m2_present, meta_step2$SampleID, drop=FALSE] + 1))
  meta_step2$M1_Score    <- m1_score
  meta_step2$M2_Score    <- m2_score
  meta_step2$M1_M2_Ratio <- m1_score - m2_score
  message("M1/M2 scores calculated.")
}

# =============================================================================
# SECTION 6: SCORE SUMMARY & STATISTICS
# =============================================================================

SCORE_COLS <- c(
  "Immune_Score","Stromal_Score","ESTIMATE_Score","Tumor_Purity",
  "Proliferation_Score","ECM_Score","Myelination_Score",
  "TIS_Score","BCell_Score","Checkpoint_Score","Angiogenesis_Score",
  "CYT_Score","HLA_ClassI_Score","Fibrosis_Score","Chemokine_Program_Score"
)
SCORE_COLS <- intersect(SCORE_COLS, colnames(meta_step2))

score_summary <- meta_step2 %>%
  select(SampleID, Subtype, all_of(SCORE_COLS)) %>%
  arrange(Subtype, SampleID)
write.csv(score_summary, file.path(OUTPUT_DIR, "VS_score_summary.csv"),
          row.names = FALSE, quote = FALSE)

message("\n=== Statistical Testing ===")
stat_results <- setNames(
  lapply(SCORE_COLS, function(sc) {
    auto_stat_test(meta_step2, sc, "Subtype")
  }),
  SCORE_COLS
)

pw_all <- do.call(rbind, lapply(SCORE_COLS, function(sc) {
  pw <- stat_results[[sc]]$pairwise_df
  pw$Score <- sc; pw$Method <- stat_results[[sc]]$method; pw
}))
write.csv(pw_all, file.path(OUTPUT_DIR, "VS_pairwise_stats.csv"),
          row.names = FALSE, quote = FALSE)

# =============================================================================
# SECTION 7: VISUALIZATION FUNCTIONS
# =============================================================================

make_raincloud <- function(score_col, y_label = NULL) {
  if (is.null(y_label)) y_label <- gsub("_", " ", score_col)
  
  df_plot <- meta_step2 %>%
    select(Subtype, Score = all_of(score_col)) %>%
    filter(!is.na(Score))
  
  pw   <- stat_results[[score_col]]$pairwise_df
  comp <- lapply(pw$comparison, parse_comparison)
  
  p <- ggplot(df_plot, aes(x = Subtype, y = Score, fill = Subtype,
                           color = Subtype)) +
    stat_halfeye(adjust = 1.2, width = 0.5, .width = 0,
                 justification = -0.2, point_colour = NA, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9,
                 color = "black", fill = "white") +
    {
      if (HAS_GGRASTR)
        ggrastr::rasterise(geom_jitter(width = 0.05, size = 0.8, alpha = 0.5),
                           dpi = 300)
      else
        geom_jitter(width = 0.05, size = 0.8, alpha = 0.5)
    } +
    geom_signif(
      comparisons = comp,
      map_signif_level = function(p) {
        if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
      },
      step_increase = 0.08, tip_length = 0.02,
      color = "black", textsize = 3
    ) +
    scale_fill_manual(values  = SUBTYPE_COLORS) +
    scale_color_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = y_label) +
    theme_pub() +
    theme(legend.position = "none")
  p
}

make_violin <- function(score_col) {
  df_v <- meta_step2 %>%
    select(Subtype, Score = all_of(score_col)) %>%
    filter(!is.na(Score))
  pw <- stat_results[[score_col]]$pairwise_df
  comp <- lapply(pw$comparison, parse_comparison)
  
  ggplot(df_v, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_violin(trim = FALSE, alpha = 0.75, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 fill = "white", color = "black") +
    geom_signif(comparisons = comp,
                map_signif_level = function(p) {
                  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
                },
                step_increase = 0.10, tip_length = 0.02,
                color = "black", textsize = 3) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = gsub("_", " ", score_col), title = gsub("_", " ", score_col)) +
    theme_pub() + theme(legend.position = "none")
}

# =============================================================================
# SECTION 8: RAINCLOUD PLOTS (All Scores)
# =============================================================================

message("\n=== Plotting: Raincloud Plots ===")

raincloud_list <- setNames(lapply(SCORE_COLS, make_raincloud), SCORE_COLS)

n_cols <- 4
n_rows <- ceiling(length(SCORE_COLS) / n_cols)
save_plot(wrap_plots(raincloud_list, ncol = n_cols),
          "Fig_Raincloud_AllScores",
          width = n_cols * 3.5, height = n_rows * 3.8)

# =============================================================================
# SECTION 9: Fig4A — ESTIMATE 3-panel Violin
# =============================================================================

message("\n=== Fig4A: ESTIMATE violin ===")

estimate_panel_scores <- intersect(
  c("Immune_Score","Stromal_Score","Tumor_Purity"), colnames(meta_step2))

fig4a_plots <- lapply(estimate_panel_scores, make_violin)

save_plot(wrap_plots(fig4a_plots, nrow = 1) +
            plot_annotation(title = "ESTIMATE Tumor Microenvironment Scores",
                            theme = theme(plot.title = element_text(
                              face = "bold", size = 12, hjust = 0.5))),
          "Fig4A_ESTIMATE_Violin", width = 10.5, height = 4.5)

# =============================================================================
# SECTION 10: Fig4B — MCP-counter Boxplot + Heatmap
# =============================================================================

message("\n=== Fig4B: MCP-counter ===")

mcp_var <- apply(meta_step2[, MCP_CELL_TYPES, drop = FALSE], 2, var, na.rm = TRUE)
TOP5_CELLS <- names(sort(mcp_var, decreasing = TRUE))[1:min(5, length(MCP_CELL_TYPES))]

fig4b_plots <- lapply(TOP5_CELLS, function(ct) {
  df_ct <- meta_step2 %>%
    select(Subtype, Score = all_of(ct)) %>% filter(!is.na(Score))
  sr <- auto_stat_test(df_ct, "Score", "Subtype")
  comp <- lapply(sr$pairwise_df$comparison, parse_comparison)
  
  ggplot(df_ct, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5, width = 0.55,
                 alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.45, color = "black") +
    geom_signif(comparisons = comp,
                map_signif_level = function(p) {
                  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
                },
                step_increase = 0.10, tip_length = 0.02,
                color = "black", textsize = 3) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = "Estimated abundance", title = gsub("_", " ", ct)) +
    theme_pub() + theme(legend.position = "none")
})

save_plot(wrap_plots(fig4b_plots, nrow = 1) +
            plot_annotation(title = "MCP-counter Immune Cell Abundances",
                            theme = theme(plot.title = element_text(
                              face = "bold", size = 12, hjust = 0.5))),
          "Fig4B_MCP_Boxplot", width = length(TOP5_CELLS) * 2.5 + 1, height = 4.5)

# MCP Heatmap
mcp_mat <- meta_step2 %>%
  arrange(Subtype, SampleID) %>%
  select(SampleID, all_of(MCP_CELL_TYPES)) %>%
  column_to_rownames("SampleID") %>% t()
mcp_z <- t(scale(t(mcp_mat)))
mcp_z[is.nan(mcp_z)] <- 0

anno_col <- data.frame(
  Subtype = meta_step2 %>% arrange(Subtype, SampleID) %>% pull(Subtype),
  row.names = meta_step2 %>% arrange(Subtype, SampleID) %>% pull(SampleID)
)
anno_colors <- list(Subtype = SUBTYPE_COLORS)

p_heatmap <- pheatmap(mcp_z,
                      annotation_col  = anno_col,
                      annotation_colors = anno_colors,
                      cluster_rows    = TRUE,
                      cluster_cols    = FALSE,
                      show_colnames   = FALSE,
                      fontsize_row    = 8,
                      color = colorRampPalette(c("#4DBBD5","white","#E64B35"))(100),
                      main = "MCP-counter Z-score",
                      silent = TRUE)
save_base_plot(quote(print(p_heatmap)), "Fig4B_MCP_Heatmap", width = 11, height = 4.5)

# =============================================================================
# SECTION 11: Fig4C — Checkpoint 5-gene Boxplot
# =============================================================================

message("\n=== Fig4C: Checkpoint genes ===")

CP5_GENES <- c("CD274","CTLA4","HAVCR2","LAG3","PDCD1")
cp5_present <- intersect(CP5_GENES, rownames(clean_tpm))

fig4c_plots <- lapply(cp5_present, function(g) {
  plot_gene_box(g, clean_tpm, meta_step2)
})
fig4c_plots <- Filter(Negate(is.null), fig4c_plots)

save_plot(wrap_plots(fig4c_plots, nrow = 1) +
            plot_annotation(title = "Immune Checkpoint Gene Expression",
                            theme = theme(plot.title = element_text(
                              face = "bold", size = 12, hjust = 0.5))),
          "Fig4C_Checkpoints", width = length(fig4c_plots) * 2.3 + 0.5, height = 4.5)

# =============================================================================
# SECTION 12: CYT + HLA Class I Violin Panel
# =============================================================================

message("\n=== CYT & HLA panel ===")

cyt_scores <- intersect(c("CYT_Score","HLA_ClassI_Score"), colnames(meta_step2))
if (length(cyt_scores) > 0) {
  cyt_plots <- lapply(cyt_scores, make_violin)
  save_plot(wrap_plots(cyt_plots, nrow = 1),
            "Fig_CYT_HLA_Violin", width = length(cyt_plots) * 3.5, height = 4.5)
}

# =============================================================================
# SECTION 13: Adaptive Resistance — CYT vs Checkpoint Scatter
# =============================================================================

message("\n=== CYT vs Checkpoint scatter ===")

if ("CYT_Score" %in% colnames(meta_step2) &&
    "Checkpoint_Score" %in% colnames(meta_step2)) {
  
  cor_df <- meta_step2 %>%
    filter(!is.na(CYT_Score), !is.na(Checkpoint_Score))
  
  cor_all <- cor.test(cor_df$CYT_Score, cor_df$Checkpoint_Score,
                      method = "spearman", exact = FALSE)
  
  p_cyt_cp <- ggplot(cor_df,
                     aes(CYT_Score, Checkpoint_Score,
                         color = Subtype, fill = Subtype)) +
    {
      if (HAS_GGRASTR)
        ggrastr::rasterise(geom_point(size = 2.5, alpha = 0.75, shape = 21,
                                      color = "white", stroke = 0.4), dpi = 300)
      else
        geom_point(size = 2.5, alpha = 0.75, shape = 21,
                   color = "white", stroke = 0.4)
    } +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Spearman r = %.2f\np = %.3f",
                             cor_all$estimate, cor_all$p.value),
             hjust = 1.1, vjust = 1.5, size = 3.5, color = "black") +
    scale_color_manual(values = SUBTYPE_COLORS) +
    scale_fill_manual(values  = SUBTYPE_COLORS) +
    labs(x = "CYT Score (Cytolytic Activity)",
         y = "Checkpoint Score",
         title = "Adaptive Immune Resistance: CYT vs Checkpoint") +
    theme_pub()
  
  save_plot(p_cyt_cp, "Fig_Adaptive_Resistance_CYT_vs_Checkpoint", width = 6.5, height = 5)
  
  # Purity-corrected residual correlation
  if ("Tumor_Purity" %in% colnames(cor_df)) {
    lm_cyt <- lm(CYT_Score ~ Tumor_Purity + Subtype, data = cor_df)
    lm_cp  <- lm(Checkpoint_Score ~ Tumor_Purity + Subtype, data = cor_df)
    cor_df$CYT_resid <- residuals(lm_cyt)
    cor_df$CP_resid  <- residuals(lm_cp)
    
    cor_resid <- cor.test(cor_df$CYT_resid, cor_df$CP_resid,
                          method = "spearman", exact = FALSE)
    p_resid <- ggplot(cor_df, aes(CYT_resid, CP_resid, color = Subtype)) +
      {
        if (HAS_GGRASTR)
          ggrastr::rasterise(geom_point(size = 2.5, alpha = 0.8), dpi = 300)
        else
          geom_point(size = 2.5, alpha = 0.8)
      } +
      geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                  linewidth = 0.9, color = "gray30", alpha = 0.15) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("Spearman r = %.2f\np = %.3f",
                               cor_resid$estimate, cor_resid$p.value),
               hjust = 1.1, vjust = 1.5, size = 3.5) +
      scale_color_manual(values = SUBTYPE_COLORS) +
      labs(x = "CYT Score (residual, purity+subtype adjusted)",
           y = "Checkpoint Score (residual)",
           title = "Purity-Corrected Adaptive Resistance") +
      theme_pub()
    
    save_plot(p_resid, "Fig_ResidualCorr_CYT_vs_Checkpoint", width = 6.5, height = 5)
  }
}

# =============================================================================
# SECTION 14: Tumor Size vs Checkpoint Paradox
# =============================================================================

message("\n=== Tumor size vs Checkpoint ===")

size_col <- if ("SizeGrade" %in% colnames(meta_step2)) "SizeGrade" else
  if ("TumorSize_cm" %in% colnames(meta_step2)) "TumorSize_cm" else NULL

if (!is.null(size_col) && "Checkpoint_Score" %in% colnames(meta_step2)) {
  paradox_df <- meta_step2 %>%
    filter(!is.na(.data[[size_col]]), !is.na(Checkpoint_Score))
  
  if (is.numeric(paradox_df[[size_col]])) {
    p_paradox <- ggplot(paradox_df,
                        aes(.data[[size_col]], Checkpoint_Score, color = Subtype)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                  color = "gray30", linewidth = 0.9) +
      scale_color_manual(values = SUBTYPE_COLORS) +
      labs(x = "Tumor Size (cm)", y = "Checkpoint Score",
           title = "Tumor Size vs Immune Checkpoint Expression") +
      theme_pub()
  } else {
    paradox_df[[size_col]] <- factor(paradox_df[[size_col]])
    p_paradox <- ggplot(paradox_df,
                        aes(.data[[size_col]], Checkpoint_Score,
                            fill = .data[[size_col]])) +
      geom_boxplot(outlier.shape = 21, width = 0.55, alpha = 0.85) +
      geom_jitter(width = 0.12, size = 1.5, alpha = 0.6, color = "black") +
      scale_fill_brewer(palette = "Blues") +
      labs(x = "Tumor Size Grade", y = "Checkpoint Score",
           title = "Tumor Size vs Immune Checkpoint (Paradox)") +
      theme_pub() + theme(legend.position = "none")
  }
  save_plot(p_paradox, "Fig_TumorSize_vs_Checkpoint_Paradox", width = 6.5, height = 5)
}

# =============================================================================
# SECTION 15: ImmuneScore vs TumorPurity Scatter
# =============================================================================

message("\n=== ImmuneScore vs TumorPurity ===")

if ("Immune_Score" %in% colnames(meta_step2) &&
    "Tumor_Purity" %in% colnames(meta_step2)) {
  
  pur_df <- meta_step2 %>%
    filter(!is.na(Tumor_Purity), !is.na(Immune_Score))
  
  cor_overall <- cor.test(pur_df$Tumor_Purity, pur_df$Immune_Score,
                          method = "spearman", exact = FALSE)
  
  p_pur <- ggplot(pur_df, aes(Tumor_Purity, Immune_Score,
                              color = Subtype, fill = Subtype)) +
    {
      if (HAS_GGRASTR)
        ggrastr::rasterise(geom_point(size = 2, alpha = 0.7, shape = 21,
                                      color = "white", stroke = 0.3), dpi = 300)
      else
        geom_point(size = 2, alpha = 0.7, shape = 21,
                   color = "white", stroke = 0.3)
    } +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Spearman r = %.2f\np = %.3f",
                             cor_overall$estimate, cor_overall$p.value),
             hjust = 1.1, vjust = 1.5, size = 3.5) +
    scale_fill_manual(values  = SUBTYPE_COLORS) +
    scale_color_manual(values = SUBTYPE_COLORS) +
    labs(x = "ESTIMATE Tumor Purity",
         y = "ESTIMATE Immune Score",
         title = "Tumor Purity vs Immune Infiltration") +
    theme_pub()
  
  save_plot(p_pur, "Fig_ImmuneScore_vs_TumorPurity", width = 6.5, height = 5)
  
  p_pur_facet <- p_pur +
    facet_wrap(~ Subtype, nrow = 1) +
    theme(legend.position = "none")
  save_plot(p_pur_facet, "Fig_ImmuneScore_vs_TumorPurity_BySubtype", width = 10, height = 4.5)
}

# =============================================================================
# SECTION 16: Checkpoint 7-gene Faceted Boxplot
# # =============================================================================
# 
# message("\n=== Checkpoint 7-gene panel ===")
# 
# CP7_GENES <- c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","IDO1")
# cp7_present <- intersect(CP7_GENES, rownames(clean_tpm))
# 
# if (length(cp7_present) >= 3) {
#   cp7_long <- do.call(rbind, lapply(cp7_present, function(g) {
#     data.frame(
#       SampleID = meta_step2$SampleID,
#       Subtype  = meta_step2$Subtype,
#       Gene     = g,
#       Expr     = log2(as.numeric(clean_tpm[g, meta_step2$SampleID]) + 1)
#     )
#   }))
#   
#   p_cp7 <- ggplot(cp7_long, aes(Subtype, Expr, fill = Subtype)) +
#     geom_boxplot(outlier.shape = 21, outlier.size = 1, width = 0.55,
#                  alpha = 0.85) +
#     geom_jitter(width = 0.12, size = 0.6, alpha = 0.4, color = "black") +
#     facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
#     stat_compare_means(
#       comparisons = list(c("C1","C2"), c("C1","C3"), c("C2","C3")),
#       method = "wilcox.test", label = "p.signif",
#       step.increase = 0.10, tip.length = 0.02,
#       size = 3, hide.ns = TRUE
#     ) +
#     scale_fill_manual(values = SUBTYPE_COLORS) +
#     labs(x = NULL, y = "log2(TPM+1)",
#          title = "Immune Checkpoint Gene Expression Panel") +
#     theme_pub() +
#     theme(legend.position = "none",
#           axis.text.x = element_text(size = 8, angle = 0),
#           strip.text  = element_text(face = "italic", size = 9))
#   
#   save_plot(p_cp7, "Fig_Checkpoint7_Panel",
#             width = length(cp7_present) * 4 + 0.5, height = 4)
# }
message("\n=== Checkpoint 7-gene panel ===")
CP7_GENES <- c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","IDO1")
cp7_present <- intersect(CP7_GENES, rownames(clean_tpm))
if (length(cp7_present) >= 3) {
  cp7_long <- do.call(rbind, lapply(cp7_present, function(g) {
    data.frame(
      SampleID = meta_step2$SampleID,
      Subtype  = meta_step2$Subtype,
      Gene     = g,
      Expr     = log2(as.numeric(clean_tpm[g, meta_step2$SampleID]) + 1)
    )
  }))
  
  p_cp7 <- ggplot(cp7_long, aes(Subtype, Expr, fill = Subtype)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1, width = 0.55,
                 alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.6, alpha = 0.4, color = "black") +
    facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
    stat_compare_means(
      comparisons = list(c("C1","C2"), c("C1","C3"), c("C2","C3")),
      method = "wilcox.test", label = "p.signif",
      step.increase = 0.10, tip.length = 0.02,
      size = 3, hide.ns = TRUE
    ) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = "log2(TPM+1)",
         title = "Immune Checkpoint Gene Expression Panel") +
    theme_pub() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 0),
          strip.text  = element_text(face = "italic", size = 9),
          # 新增代码：添加四周全包围边框并移除可能冲突的单侧轴线
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.line    = element_blank() 
    )
  
  save_plot(p_cp7, "Fig_Checkpoint7_Panel",
            width = length(cp7_present) * 3.5 + 0.5, height = 4)
}
# =============================================================================
# SECTION 17: Expanded Checkpoint Panel
# =============================================================================

message("\n=== Expanded checkpoint panel ===")

EXP_CP_GENES <- c("LGALS9","PVR","CD276","VSIR","ENTPD1","NT5E",
                  "ADORA2A","CD47","SIGLEC7")
exp_present <- intersect(EXP_CP_GENES, rownames(clean_tpm))

if (length(exp_present) >= 3) {
  exp_plots <- lapply(exp_present, function(g) {
    plot_gene_box(g, clean_tpm, meta_step2)
  })
  exp_plots <- Filter(Negate(is.null), exp_plots)
  n_exp_cols <- min(5, length(exp_plots))
  save_plot(
    wrap_plots(exp_plots, ncol = n_exp_cols) +
      plot_annotation(title = "Expanded Checkpoint Panel",
                      theme = theme(plot.title = element_text(
                        face = "bold", size = 12, hjust = 0.5))),
    "Fig_ExpandedCheckpoint_Panel",
    width = n_exp_cols * 2.3 + 0.5,
    height = ceiling(length(exp_plots) / n_exp_cols) * 4
  )
}

# =============================================================================
# SECTION 18: C2 Fibrosis Marker Panel
# =============================================================================

message("\n=== C2 Fibrosis panel ===")

FIB_GENES <- c("COL1A1","POSTN","FN1","ACTA2","FAP",
               "PDPN","COL10A1","CTHRC1","LOXL2","PRRX1")
fib_gene_present <- intersect(FIB_GENES, rownames(clean_tpm))

if (length(fib_gene_present) >= 3) {
  fib_plots <- lapply(fib_gene_present, function(g)
    plot_gene_box(g, clean_tpm, meta_step2))
  fib_plots <- Filter(Negate(is.null), fib_plots)
  n_fib <- min(5, length(fib_plots))
  save_plot(
    wrap_plots(fib_plots, ncol = n_fib) +
      plot_annotation(title = "C2 Fibrosis/CAF Markers",
                      theme = theme(plot.title = element_text(
                        face = "bold", size = 12, hjust = 0.5))),
    "Fig_C2_Fibrosis_Markers",
    width = n_fib * 2.3 + 0.5,
    height = ceiling(length(fib_plots) / n_fib) * 4
  )
}

if ("Fibrosis_Score" %in% colnames(meta_step2)) {
  save_plot(make_raincloud("Fibrosis_Score", "Fibrosis Score (ssGSEA)"),
            "Fig_C2_Fibrosis_Score", width = 4.5, height = 4.5)
}

# =============================================================================
# SECTION 19: C2 Angiogenesis Marker Panel
# =============================================================================

message("\n=== C2 Angiogenesis panel ===")

ANGIO_MARKER_GENES <- c("PECAM1","VWF","KDR","TEK","FLT1",
                        "FLT4","CDH5","VEGFA","ANGPT2","NRP1")
angio_m_present <- intersect(ANGIO_MARKER_GENES, rownames(clean_tpm))

if (length(angio_m_present) >= 3) {
  angio_plots <- lapply(angio_m_present, function(g)
    plot_gene_box(g, clean_tpm, meta_step2))
  angio_plots <- Filter(Negate(is.null), angio_plots)
  n_ang <- min(5, length(angio_plots))
  save_plot(
    wrap_plots(angio_plots, ncol = n_ang) +
      plot_annotation(title = "C2 Angiogenesis Markers",
                      theme = theme(plot.title = element_text(
                        face = "bold", size = 12, hjust = 0.5))),
    "Fig_C2_Angiogenesis_Markers",
    width = n_ang * 2.3 + 0.5,
    height = ceiling(length(angio_plots) / n_ang) * 4
  )
}

# =============================================================================
# SECTION 20: C3 Chemokine Mechanism Panel
# =============================================================================

message("\n=== C3 Chemokine panel ===")

CHEMO_GENES <- c("CXCL9","CXCL10","CXCL11","CCL5","CCL19","CCL21",
                 "CD86","MRC1","HAVCR2","LGALS9","ENTPD1")
chemo_present <- intersect(CHEMO_GENES, rownames(clean_tpm))

if (length(chemo_present) >= 3) {
  chemo_plots <- lapply(chemo_present, function(g)
    plot_gene_box(g, clean_tpm, meta_step2))
  chemo_plots <- Filter(Negate(is.null), chemo_plots)
  n_ch <- min(5, length(chemo_plots))
  save_plot(
    wrap_plots(chemo_plots, ncol = n_ch) +
      plot_annotation(title = "C3 Chemokine & Myeloid Regulation Markers",
                      theme = theme(plot.title = element_text(
                        face = "bold", size = 12, hjust = 0.5))),
    "Fig_C3_Chemokine_Mechanism",
    width = n_ch * 2.3 + 0.5,
    height = ceiling(length(chemo_plots) / n_ch) * 4
  )
}

if ("Chemokine_Program_Score" %in% colnames(meta_step2)) {
  save_plot(make_raincloud("Chemokine_Program_Score",
                           "Chemokine Program Score (ssGSEA)"),
            "Fig_C3_Chemokine_Program_Score", width = 4.5, height = 4.5)
}

# =============================================================================
# SECTION 21: M1/M2 Macrophage Polarization Ratio
# =============================================================================

message("\n=== M1/M2 polarization ratio ===")

if ("M1_M2_Ratio" %in% colnames(meta_step2)) {
  stat_m1m2 <- auto_stat_test(
    meta_step2 %>% select(Subtype, M1_M2_Ratio) %>%
      rename(score = M1_M2_Ratio), "score", "Subtype"
  )
  comp_m1m2 <- lapply(stat_m1m2$pairwise_df$comparison, parse_comparison)
  
  p_m1m2 <- ggplot(
    meta_step2 %>% filter(!is.na(M1_M2_Ratio)),
    aes(Subtype, M1_M2_Ratio, fill = Subtype)
  ) +
    geom_boxplot(outlier.shape = 21, width = 0.55, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_signif(comparisons = comp_m1m2,
                map_signif_level = function(p) {
                  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
                },
                step_increase = 0.10, tip_length = 0.02,
                color = "black", textsize = 3) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = "M1 - M2 Score (log2-space)",
         title = "Macrophage Polarization (M1 vs M2)") +
    theme_pub() + theme(legend.position = "none")
  
  save_plot(p_m1m2, "Fig_Macrophage_M1M2_Ratio", width = 4.5, height = 5)
}

# =============================================================================
# SECTION 22: Checkpoint Mini-Heatmap
# =============================================================================

message("\n=== Checkpoint mini-heatmap ===")

CP4_GENES <- c("PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2")
cp4_present <- intersect(CP4_GENES, rownames(clean_tpm))

if (length(cp4_present) >= 3) {
  cp4_mat <- as.matrix(log2(clean_tpm[cp4_present,
                                      meta_step2 %>%
                                        arrange(Subtype, SampleID) %>%
                                        pull(SampleID), drop = FALSE] + 1))
  cp4_z   <- t(scale(t(cp4_mat)))
  cp4_z[is.nan(cp4_z)] <- 0
  
  anno_c <- data.frame(
    Subtype = meta_step2 %>% arrange(Subtype, SampleID) %>% pull(Subtype),
    row.names = meta_step2 %>% arrange(Subtype, SampleID) %>% pull(SampleID)
  )
  
  p_cp_heat <- pheatmap(
    cp4_z,
    annotation_col   = anno_c,
    annotation_colors = list(Subtype = SUBTYPE_COLORS),
    cluster_rows     = TRUE,
    cluster_cols     = FALSE,
    show_colnames    = FALSE,
    fontsize_row     = 9,
    fontsize         = 8,
    color = colorRampPalette(c("#4DBBD5","white","#E64B35"))(100),
    main  = "Checkpoint Gene Expression (Row Z-score)",
    silent = TRUE
  )
  save_base_plot(quote(print(p_cp_heat)),
                 "Fig_Checkpoint_MiniHeatmap", width = 10, height = 3.5)
}

# =============================================================================
# SECTION 23: Checkpoint Point-Range Summary
# =============================================================================

message("\n=== Checkpoint point-range summary ===")

if (length(cp4_present) >= 3) {
  cp4_long <- do.call(rbind, lapply(cp4_present, function(g) {
    data.frame(
      Subtype = meta_step2$Subtype,
      Gene    = g,
      Expr    = log2(as.numeric(clean_tpm[g, meta_step2$SampleID]) + 1)
    )
  })) %>%
    group_by(Subtype, Gene) %>%
    summarise(Median = median(Expr, na.rm = TRUE),
              Q25    = quantile(Expr, 0.25, na.rm = TRUE),
              Q75    = quantile(Expr, 0.75, na.rm = TRUE),
              .groups = "drop")
  
  p_pr <- ggplot(cp4_long,
                 aes(x = Median, y = Gene, color = Subtype,
                     xmin = Q25, xmax = Q75)) +
    geom_linerange(position = position_dodge(width = 0.5),
                   linewidth = 1.2, alpha = 0.8) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    scale_color_manual(values = SUBTYPE_COLORS) +
    labs(x = "log2(TPM+1) — Median ± IQR", y = NULL,
         title = "Checkpoint Expression Summary") +
    theme_pub()
  
  save_plot(p_pr, "Fig_Checkpoint_PointRange_Summary", width = 6.5, height = 3.5)
}

# =============================================================================
# SECTION 24: TME Landscape (Stacked Bar)
# =============================================================================

message("\n=== TME Landscape ===")

mcp_long <- meta_step2 %>%
  select(SampleID, Subtype, all_of(MCP_CELL_TYPES)) %>%
  pivot_longer(cols = all_of(MCP_CELL_TYPES),
               names_to = "CellType", values_to = "Score") %>%
  group_by(SampleID) %>%
  mutate(Score = pmax(Score, 0),
         Prop  = Score / (sum(Score, na.rm = TRUE) + EPSILON)) %>%
  ungroup()
mcp_long$CellType <- gsub("_", " ", mcp_long$CellType)

n_ct     <- length(unique(mcp_long$CellType))
ct_pal   <- colorRampPalette(brewer.pal(min(n_ct, 9), "Set1"))(n_ct)
samp_ord <- meta_step2 %>% arrange(Subtype, SampleID) %>% pull(SampleID)
mcp_long$SampleID <- factor(mcp_long$SampleID, levels = samp_ord)

sub_bounds <- meta_step2 %>% arrange(Subtype, SampleID) %>%
  group_by(Subtype) %>% summarise(n = n()) %>%
  mutate(cumN = cumsum(n))

p_tme <- ggplot(mcp_long, aes(SampleID, Prop, fill = CellType)) +
  {
    if (HAS_GGRASTR)
      ggrastr::rasterise(geom_bar(stat = "identity", width = 1), dpi = 300)
    else
      geom_bar(stat = "identity", width = 1)
  } +
  geom_vline(xintercept = head(sub_bounds$cumN, -1) + 0.5,
             linetype = "dashed", color = "white", linewidth = 0.8) +
  scale_fill_manual(values = ct_pal) +
  scale_y_continuous(labels = percent_format()) +
  annotate("text",
           x = sub_bounds$cumN - sub_bounds$n / 2,
           y = 1.04,
           label = as.character(sub_bounds$Subtype),
           color = SUBTYPE_COLORS[as.character(sub_bounds$Subtype)],
           fontface = "bold", size = 4) +
  labs(x = "Samples", y = "Estimated Proportion",
       title = "Tumor Microenvironment Landscape", fill = "Cell Type") +
  theme_pub() +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm"))

save_plot(p_tme, "Fig_TME_Landscape", width = 12, height = 5)

# =============================================================================
# SECTION 25: Ternary Plot (ggplot2 Manual Projection)
# =============================================================================

message("\n=== Ternary Plot ===")

tern_scores <- c("Proliferation_Score","Myelination_Score","TIS_Score")

if (all(tern_scores %in% colnames(meta_step2))) {
  ternary_df <- meta_step2 %>%
    select(SampleID, Subtype,
           Prolif = Proliferation_Score,
           Myelin = Myelination_Score,
           Immune = TIS_Score) %>%
    filter(!is.na(Prolif), !is.na(Myelin), !is.na(Immune)) %>%
    mutate(
      Prolif = Prolif - min(Prolif) + EPSILON,
      Myelin = Myelin - min(Myelin) + EPSILON,
      Immune = Immune - min(Immune) + EPSILON,
      Total  = Prolif + Myelin + Immune,
      tc1 = Prolif / Total,
      tc2 = Myelin / Total,
      tc3 = Immune / Total,
      x = tc2 + 0.5 * tc3,
      y = tc3 * (sqrt(3) / 2)
    )
  
  triangle_df <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, sqrt(3)/2, 0)
  )
  
  label_df <- data.frame(
    x = c(0, 1, 0.5),
    y = c(-0.04, -0.04, sqrt(3)/2 + 0.04),
    label = c("Proliferation", "Myelination", "Immune"),
    hjust = c(1, 0, 0.5),
    vjust = c(1, 1, 0)
  )
  
  p_tern <- ggplot(ternary_df, aes(x = x, y = y, color = Subtype, fill = Subtype)) +
    geom_polygon(data = triangle_df, aes(x = x, y = y), 
                 fill = NA, color = "black", inherit.aes = FALSE, linewidth = 0.6) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, 
                 linetype = 2, linewidth = 0.4) +
    {
      if (HAS_GGRASTR)
        ggrastr::rasterise(geom_point(size = 2.5, alpha = 0.8, shape = 21,
                                      stroke = 0.3, color = "white"), dpi = 300)
      else
        geom_point(size = 2.5, alpha = 0.8, shape = 21,
                   stroke = 0.3, color = "white")
    } +
    geom_text(data = label_df, aes(x = x, y = y, label = label, 
                                   hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 4.5, fontface = "bold", color = "black") +
    scale_color_manual(values = SUBTYPE_COLORS) +
    scale_fill_manual(values  = SUBTYPE_COLORS) +
    coord_fixed(ratio = 1, clip = "off") + 
    labs(title = "VS Subtype Functional Axes") +
    theme_void() + 
    theme(
      text = element_text(family = BASE_FAMILY),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 15)),
      legend.position = "right",
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20) 
    )
  
  save_plot(p_tern, "Fig_Ternary_FunctionalAxes", width = 7, height = 6)
}

# =============================================================================
# SECTION 26: Radar Chart
# =============================================================================

message("\n=== Radar Chart ===")

RADAR_SCORES <- intersect(
  c("TIS_Score","BCell_Score","Checkpoint_Score","Proliferation_Score",
    "ECM_Score","Myelination_Score","Angiogenesis_Score","Immune_Score",
    "CYT_Score","HLA_ClassI_Score"),
  colnames(meta_step2)
)

if (length(RADAR_SCORES) >= 3) {
  radar_meds <- meta_step2 %>%
    group_by(Subtype) %>%
    summarise(across(all_of(RADAR_SCORES), ~median(.x, na.rm = TRUE))) %>%
    column_to_rownames("Subtype")
  
  radar_sc <- as.data.frame(
    apply(radar_meds, 2, function(x) (x - min(x)) / (max(x) - min(x) + EPSILON)))
  radar_fmsb <- rbind(rep(1, ncol(radar_sc)), rep(0, ncol(radar_sc)), radar_sc)
  colnames(radar_fmsb) <- trimws(gsub("_(Score)?", " ", colnames(radar_sc)))
  
  save_base_plot(quote({
    par(mar = c(1,1,2,1))
    radarchart(radar_fmsb,
               axistype   = 1,
               pcol       = SUBTYPE_COLORS[c("C1","C2","C3")],
               pfcol      = adjustcolor(SUBTYPE_COLORS[c("C1","C2","C3")],
                                        alpha.f = 0.15),
               plwd = 2, cglcol = "grey70", cglty = 1, cglwd = 0.4,
               axislabcol = "grey30", vlcex = 0.85,
               title = "Multi-Dimensional Subtype Profile")
    legend("topright", legend = paste0("C", 1:3),
           col = SUBTYPE_COLORS, lty = 1, lwd = 2, bty = "n", cex = 0.9)
  }), "Fig_Radar_SubtypeProfile", width = 7, height = 6)
}

# =============================================================================
# SECTION 27: Correlation Bubble Matrix
# =============================================================================

message("\n=== Correlation Matrix ===")

score_mat_cor <- meta_step2 %>% select(all_of(SCORE_COLS)) %>% drop_na()

if (nrow(score_mat_cor) >= 10) {
  cor_mat  <- cor(score_mat_cor, method = "spearman")
  cor_pmat <- matrix(NA, nrow = ncol(score_mat_cor), ncol = ncol(score_mat_cor),
                     dimnames = list(colnames(score_mat_cor),
                                     colnames(score_mat_cor)))
  for (i in seq_len(ncol(score_mat_cor))) {
    for (j in seq_len(ncol(score_mat_cor))) {
      if (i != j) {
        cor_pmat[i,j] <- cor.test(score_mat_cor[[i]], score_mat_cor[[j]],
                                  method = "spearman", exact = FALSE)$p.value
      } else {
        cor_pmat[i,j] <- 0
      }
    }
  }
  
  save_base_plot(quote({
    corrplot(cor_mat, method = "circle", type = "upper", order = "hclust",
             tl.col = "black", tl.cex = 0.8, tl.srt = 45,
             addCoef.col = "black", number.cex = 0.6,
             p.mat = cor_pmat, sig.level = 0.05, insig = "blank",
             col = colorRampPalette(c("#4DBBD5","white","#E64B35"))(200),
             cl.cex = 0.7,
             title = "Score Correlation Matrix (Spearman, p<0.05)",
             mar = c(0,0,2,0), tl.offset = 0.5)
  }), "Fig_Correlation_BubbleMatrix", width = 9, height = 8)
}

# =============================================================================
# SECTION 28: MCP-Counter ComplexHeatmap
# =============================================================================

message("\n=== MCP ComplexHeatmap ===")

mcp_mat2 <- meta_step2 %>%
  arrange(Subtype, SampleID) %>%
  select(SampleID, all_of(MCP_CELL_TYPES)) %>%
  column_to_rownames("SampleID") %>% t()
mcp_z2 <- t(scale(t(mcp_mat2)))
mcp_z2[is.nan(mcp_z2)] <- 0

sub_ord2 <- meta_step2 %>% arrange(Subtype, SampleID) %>%
  pull(Subtype) %>% as.character()

top_ann <- HeatmapAnnotation(
  Subtype = sub_ord2,
  col = list(Subtype = SUBTYPE_COLORS),
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm")
)
cf2 <- colorRamp2(c(-2,0,2), c("#4DBBD5","white","#E64B35"))

mcp_ht <- Heatmap(mcp_z2, name = "Z-score", col = cf2,
                  top_annotation  = top_ann,
                  cluster_rows    = TRUE, cluster_columns = FALSE,
                  show_column_names = FALSE,
                  row_names_gp      = gpar(fontsize = 9, fontfamily = BASE_FAMILY),
                  column_split      = sub_ord2,
                  column_title_gp   = gpar(fontsize = 10, fontface = "bold",
                                           fontfamily = BASE_FAMILY),
                  column_gap        = unit(2, "mm"),
                  border            = FALSE,
                  use_raster        = TRUE,
                  raster_quality    = 2)

save_base_plot(quote(draw(mcp_ht)), "Fig_Heatmap_MCPcounter", width = 12, height = 5.5)

# =============================================================================
# SECTION 29: Parallel Coordinates
# =============================================================================

message("\n=== Parallel Coordinates ===")

score_norm <- meta_step2 %>%
  select(SampleID, Subtype, all_of(SCORE_COLS)) %>%
  drop_na() %>%
  mutate(across(all_of(SCORE_COLS),
                ~ (. - min(., na.rm = TRUE)) /
                  (max(., na.rm = TRUE) - min(., na.rm = TRUE) + EPSILON)))

score_long2 <- score_norm %>%
  pivot_longer(all_of(SCORE_COLS), names_to = "Score", values_to = "Value") %>%
  mutate(Score = factor(Score, levels = SCORE_COLS))

p_par <- ggplot(score_long2, aes(Score, Value, group = SampleID,
                                 color = Subtype)) +
  {
    if (HAS_GGRASTR)
      ggrastr::rasterise(geom_line(alpha = 0.2, linewidth = 0.3), dpi = 300)
    else
      geom_line(alpha = 0.2, linewidth = 0.3)
  } +
  stat_summary(aes(group = Subtype), fun = median, geom = "line",
               linewidth = 1.8, alpha = 0.95) +
  scale_color_manual(values = SUBTYPE_COLORS) +
  scale_x_discrete(labels = setNames(
    gsub("_Score|_", "\n", SCORE_COLS), SCORE_COLS)) +
  labs(x = NULL, y = "Normalized Score (0-1)",
       title = "Per-Sample Multi-Score Profile") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 8, hjust = 0.5))

save_plot(p_par, "Fig_Parallel_SampleProfile", width = 11, height = 5)

# =============================================================================
# SECTION 30: Integrated Figure Panels
# =============================================================================

message("\n=== Integrated Figure Panels ===")

est_sc <- intersect(c("Immune_Score","Stromal_Score","ESTIMATE_Score",
                      "Tumor_Purity"), colnames(meta_step2))
if (length(est_sc) > 0) {
  p_est <- wrap_plots(lapply(est_sc, make_raincloud), ncol = 4) +
    plot_annotation(title = "ESTIMATE Scores",
                    theme = theme(plot.title = element_text(
                      face = "bold", size = 12, hjust = 0.5)))
  save_plot(p_est, "Fig_Panel_ESTIMATE", width = 14, height = 4.5)
}

func_sc <- intersect(
  c("Proliferation_Score","ECM_Score","Myelination_Score","TIS_Score",
    "BCell_Score","Checkpoint_Score","Angiogenesis_Score","CYT_Score",
    "HLA_ClassI_Score","Fibrosis_Score"),
  colnames(meta_step2)
)
if (length(func_sc) > 0) {
  p_func <- wrap_plots(lapply(func_sc, make_raincloud), ncol = 4) +
    plot_annotation(title = "Subtype-Specific Functional Scores",
                    theme = theme(plot.title = element_text(
                      face = "bold", size = 12, hjust = 0.5)))
  save_plot(p_func, "Fig_Panel_FunctionalScores",
            width = 14, height = ceiling(length(func_sc) / 4) * 4)
}

# =============================================================================
# SECTION 31: SESSION INFO & CLEANUP
# =============================================================================

message("\n=== Analysis Complete ===")
message("All outputs saved to: ", OUTPUT_DIR)

sink(file.path(OUTPUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()

# # Final cleanup
# rm(list = setdiff(ls(), c("meta_step2", "SUBTYPE_COLORS", "OUTPUT_DIR", "SCORE_COLS")))
# gc(verbose = FALSE)
# 
# message("\n=== Done! ===")
# =============================================================================
# FIXED PLOTS - 修复空白图和 p 值显示问题
# =============================================================================

message("\n=== 修复空白图和 p 值显示 ===")

# --- 辅助函数：格式化 p 值（避免显示 0.00）---
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return(sprintf("p < 0.0001"))
  if (p < 0.001)  return(sprintf("p = %.2e", p))
  return(sprintf("p = %.4f", p))
}

# --- 修复1: Fig_Raincloud_AllScores (用 violin + boxplot 替代 ggdist) ---
message("修复 Fig_Raincloud_AllScores...")

make_raincloud_fixed <- function(score_col, y_label = NULL) {
  if (is.null(y_label)) y_label <- gsub("_", " ", score_col)
  
  df_plot <- meta_step2 %>%
    dplyr::select(Subtype, Score = all_of(score_col)) %>%
    dplyr::filter(!is.na(Score))
  
  if (nrow(df_plot) == 0) return(NULL)
  
  # 计算统计
  dt <- FSA::dunnTest(Score ~ Subtype, data = df_plot, method = "bh")
  comp <- lapply(dt$res$Comparison, function(x) trimws(strsplit(x, "\\s*-\\s*")[[1]])[1:2])
  
  ggplot(df_plot, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_violin(trim = FALSE, alpha = 0.6, scale = "width", color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.95, 
                 color = "black", fill = "white") +
    geom_jitter(width = 0.08, size = 0.6, alpha = 0.5, color = "black") +
    geom_signif(comparisons = comp, 
                map_signif_level = function(p) {
                  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
                }, 
                step_increase = 0.08, tip_length = 0.02, color = "black", textsize = 3) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = y_label) + 
    theme_pub() +
    theme(legend.position = "none")
}

raincloud_fixed <- lapply(SCORE_COLS, function(sc) {
  tryCatch(make_raincloud_fixed(sc), error = function(e) NULL)
})
raincloud_fixed <- Filter(Negate(is.null), raincloud_fixed)

if (length(raincloud_fixed) > 0) {
  n_cols <- 4
  n_rows <- ceiling(length(raincloud_fixed) / n_cols)
  save_plot(wrap_plots(raincloud_fixed, ncol = n_cols),
            "Fig_Raincloud_AllScores",
            width = n_cols * 3.5, height = n_rows * 3.8)
}

# --- 修复2: Fig_C2_Fibrosis_Score ---
message("修复 Fig_C2_Fibrosis_Score...")

if ("Fibrosis_Score" %in% colnames(meta_step2)) {
  p_fib <- make_raincloud_fixed("Fibrosis_Score", "Fibrosis Score (ssGSEA)")
  if (!is.null(p_fib)) {
    save_plot(p_fib, "Fig_C2_Fibrosis_Score", width = 4.5, height = 4.5)
  }
}

# --- 修复3: Fig_C3_Chemokine_Program_Score ---
message("修复 Fig_C3_Chemokine_Program_Score...")

if ("Chemokine_Program_Score" %in% colnames(meta_step2)) {
  p_chemo <- make_raincloud_fixed("Chemokine_Program_Score", "Chemokine Program Score (ssGSEA)")
  if (!is.null(p_chemo)) {
    save_plot(p_chemo, "Fig_C3_Chemokine_Program_Score", width = 4.5, height = 4.5)
  }
}

# --- 修复4: Fig_Panel_ESTIMATE ---
message("修复 Fig_Panel_ESTIMATE...")

est_sc <- intersect(c("Immune_Score","Stromal_Score","ESTIMATE_Score","Tumor_Purity"), 
                    colnames(meta_step2))
if (length(est_sc) > 0) {
  est_plots <- lapply(est_sc, function(sc) {
    tryCatch(make_raincloud_fixed(sc), error = function(e) NULL)
  })
  est_plots <- Filter(Negate(is.null), est_plots)
  if (length(est_plots) > 0) {
    p_est <- wrap_plots(est_plots, ncol = 4) + 
      plot_annotation(title = "ESTIMATE Scores", 
                      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))
    save_plot(p_est, "Fig_Panel_ESTIMATE", width = 14, height = 4.5)
  }
}

# --- 修复5: Fig_Panel_FunctionalScores ---
message("修复 Fig_Panel_FunctionalScores...")

func_sc <- intersect(c("Proliferation_Score","ECM_Score","Myelination_Score","TIS_Score",
                       "BCell_Score","Checkpoint_Score","Angiogenesis_Score","CYT_Score",
                       "HLA_ClassI_Score","Fibrosis_Score"), colnames(meta_step2))
if (length(func_sc) > 0) {
  func_plots <- lapply(func_sc, function(sc) {
    tryCatch(make_raincloud_fixed(sc), error = function(e) NULL)
  })
  func_plots <- Filter(Negate(is.null), func_plots)
  if (length(func_plots) > 0) {
    p_func <- wrap_plots(func_plots, ncol = 4) + 
      plot_annotation(title = "Subtype-Specific Functional Scores", 
                      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))
    save_plot(p_func, "Fig_Panel_FunctionalScores", 
              width = 14, height = ceiling(length(func_plots) / 4) * 4)
  }
}

# --- 修复6: 相关性散点图 p 值显示（避免 p = 0.00）---
message("修复相关性图的 p 值显示...")

# CYT vs Checkpoint
if ("CYT_Score" %in% colnames(meta_step2) && "Checkpoint_Score" %in% colnames(meta_step2)) {
  
  cor_df <- meta_step2 %>%
    dplyr::filter(!is.na(CYT_Score), !is.na(Checkpoint_Score))
  
  cor_all <- cor.test(cor_df$CYT_Score, cor_df$Checkpoint_Score,
                      method = "spearman", exact = FALSE)
  
  p_cyt_cp <- ggplot(cor_df, aes(CYT_Score, Checkpoint_Score, color = Subtype, fill = Subtype)) +
    {
      if (HAS_GGRASTR)
        ggrastr::rasterise(geom_point(size = 2.5, alpha = 0.75, shape = 21,
                                      color = "white", stroke = 0.4), dpi = 300)
      else
        geom_point(size = 2.5, alpha = 0.75, shape = 21, color = "white", stroke = 0.4)
    } +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Spearman r = %.2f\n%s", 
                             cor_all$estimate, format_pval(cor_all$p.value)),
             hjust = 1.1, vjust = 1.5, size = 3.5, color = "black") +
    scale_color_manual(values = SUBTYPE_COLORS) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = "CYT Score (Cytolytic Activity)", y = "Checkpoint Score",
         title = "Adaptive Immune Resistance: CYT vs Checkpoint") +
    theme_pub()
  
  save_plot(p_cyt_cp, "Fig_Adaptive_Resistance_CYT_vs_Checkpoint", width = 6.5, height = 5)
  
  # Purity-corrected residual
  if ("Tumor_Purity" %in% colnames(cor_df)) {
    lm_cyt <- lm(CYT_Score ~ Tumor_Purity + Subtype, data = cor_df)
    lm_cp  <- lm(Checkpoint_Score ~ Tumor_Purity + Subtype, data = cor_df)
    cor_df$CYT_resid <- residuals(lm_cyt)
    cor_df$CP_resid  <- residuals(lm_cp)
    
    cor_resid <- cor.test(cor_df$CYT_resid, cor_df$CP_resid, method = "spearman", exact = FALSE)
    
    p_resid <- ggplot(cor_df, aes(CYT_resid, CP_resid, color = Subtype)) +
      {
        if (HAS_GGRASTR)
          ggrastr::rasterise(geom_point(size = 2.5, alpha = 0.8), dpi = 300)
        else
          geom_point(size = 2.5, alpha = 0.8)
      } +
      geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                  linewidth = 0.9, color = "gray30", alpha = 0.15) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("Spearman r = %.2f\n%s", 
                               cor_resid$estimate, format_pval(cor_resid$p.value)),
               hjust = 1.1, vjust = 1.5, size = 3.5) +
      scale_color_manual(values = SUBTYPE_COLORS) +
      labs(x = "CYT Score (residual, purity+subtype adjusted)",
           y = "Checkpoint Score (residual)",
           title = "Purity-Corrected Adaptive Resistance") +
      theme_pub()
    
    save_plot(p_resid, "Fig_ResidualCorr_CYT_vs_Checkpoint", width = 6.5, height = 5)
  }
}

# ImmuneScore vs TumorPurity
if ("Immune_Score" %in% colnames(meta_step2) && "Tumor_Purity" %in% colnames(meta_step2)) {
  
  pur_df <- meta_step2 %>%
    dplyr::filter(!is.na(Tumor_Purity), !is.na(Immune_Score))
  
  cor_overall <- cor.test(pur_df$Tumor_Purity, pur_df$Immune_Score,
                          method = "spearman", exact = FALSE)
  
  p_pur <- ggplot(pur_df, aes(Tumor_Purity, Immune_Score, color = Subtype, fill = Subtype)) +
    {
      if (HAS_GGRASTR)
        ggrastr::rasterise(geom_point(size = 2, alpha = 0.7, shape = 21,
                                      color = "white", stroke = 0.3), dpi = 300)
      else
        geom_point(size = 2, alpha = 0.7, shape = 21, color = "white", stroke = 0.3)
    } +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Spearman r = %.2f\n%s", 
                             cor_overall$estimate, format_pval(cor_overall$p.value)),
             hjust = 1.1, vjust = 1.5, size = 3.5) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    scale_color_manual(values = SUBTYPE_COLORS) +
    labs(x = "ESTIMATE Tumor Purity", y = "ESTIMATE Immune Score",
         title = "Tumor Purity vs Immune Infiltration") +
    theme_pub()
  
  save_plot(p_pur, "Fig_ImmuneScore_vs_TumorPurity", width = 6.5, height = 5)
  
  # Faceted
  p_pur_facet <- p_pur +
    facet_wrap(~ Subtype, nrow = 1) +
    theme(legend.position = "none")
  save_plot(p_pur_facet, "Fig_ImmuneScore_vs_TumorPurity_BySubtype", width = 10, height = 4.5)
}

message("\n=== 图片修复完成 ===")






# =============================================================================
# SELF-CONTAINED FIX — 修复空白图 + p 值格式
# 不依赖任何之前定义的函数，全部重新定义
# =============================================================================

message("\n========== SELF-CONTAINED FIX START ==========\n")

# --- 1. 加载包 ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
  library(patchwork)
  library(dplyr)
  library(grid)
})

# --- 2. 恢复可能被 rm() 清除的全局变量 ---
if (!exists("OUTPUT_DIR")) {
  OUTPUT_DIR <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/step8_Immune"
}
if (!exists("SUBTYPE_COLORS")) {
  SUBTYPE_COLORS <- c(C1 = "#00A087", C2 = "#4DBBD5", C3 = "#E64B35")
}
stopifnot("meta_step2 not found in environment" = exists("meta_step2"))

# 重建 SCORE_COLS
SCORE_COLS_FIX <- intersect(
  c("Immune_Score","Stromal_Score","ESTIMATE_Score","Tumor_Purity",
    "Proliferation_Score","ECM_Score","Myelination_Score",
    "TIS_Score","BCell_Score","Checkpoint_Score","Angiogenesis_Score",
    "CYT_Score","HLA_ClassI_Score","Fibrosis_Score","Chemokine_Program_Score"),
  colnames(meta_step2)
)

# 确保 Subtype 是正确的 factor
meta_step2$Subtype <- factor(meta_step2$Subtype, levels = c("C1","C2","C3"))

# 诊断信息
message("Samples: ", nrow(meta_step2))
message("Subtypes: ", paste(names(table(meta_step2$Subtype)),
                            table(meta_step2$Subtype), sep = "=", collapse = ", "))
message("Score columns: ", length(SCORE_COLS_FIX))
for (sc in SCORE_COLS_FIX) {
  n_ok <- sum(!is.na(meta_step2[[sc]]))
  message(sprintf("  %-30s n_valid = %d", sc, n_ok))
}

# --- 3. 自包含主题 ---
theme_fix <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line        = element_line(linewidth = 0.5, colour = "black"),
      axis.text        = element_text(color = "black", size = base_size - 1),
      axis.title       = element_text(color = "black", size = base_size),
      legend.position  = "none",
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size),
      panel.grid       = element_blank()
    )
}

# --- 4. 自包含 PDF 保存（带诊断）---
save_pdf_fix <- function(plot_obj, filename, width = 7, height = 5) {
  fpath <- file.path(OUTPUT_DIR, paste0(filename, ".pdf"))
  
  # 预渲染测试：如果 plot 构建失败，这里就会报错而不是生成空白 PDF
  grob <- tryCatch(
    ggplotGrob(plot_obj),
    error = function(e) {
      message("!! PLOT BUILD FAILED for ", filename, ": ", e$message)
      NULL
    }
  )
  if (is.null(grob)) return(invisible(NULL))
  
  pdf(fpath, width = width, height = height)
  grid.draw(grob)
  dev.off()
  
  fsize <- file.size(fpath)
  if (fsize < 2000) {
    message("!! WARNING: ", filename, ".pdf is only ", fsize, " bytes — likely blank!")
  } else {
    message("OK: ", filename, ".pdf (", fsize, " bytes)")
  }
}

# --- 5. p 值格式化函数 ---
format_pval <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  if (p < 0.0001)  return(sprintf("p = %.2e", p))
  if (p < 0.001)   return(sprintf("p = %.1e", p))
  sprintf("p = %.4f", p)
}

# --- 6. 核心绘图函数：单个 score 的 violin+box ---
make_viobox <- function(score_col, ylab = NULL) {
  if (is.null(ylab)) ylab <- gsub("_", " ", score_col)
  
  df <- data.frame(
    Subtype = meta_step2$Subtype,
    Score   = as.numeric(meta_step2[[score_col]])
  )
  df <- df[complete.cases(df), ]
  df$Subtype <- droplevels(df$Subtype)
  
  # 至少需要两个组、每组至少 3 个样本
  grp_n <- table(df$Subtype)
  valid_groups <- names(grp_n[grp_n >= 3])
  
  if (length(valid_groups) < 2) {
    message("  SKIP ", score_col, ": fewer than 2 valid groups")
    return(NULL)
  }
  
  df <- df[df$Subtype %in% valid_groups, ]
  df$Subtype <- factor(df$Subtype, levels = intersect(c("C1","C2","C3"), valid_groups))
  
  # 动态生成比较对
  all_pairs <- list(c("C1","C2"), c("C1","C3"), c("C2","C3"))
  comps <- Filter(function(pair) all(pair %in% valid_groups), all_pairs)
  
  p <- ggplot(df, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_violin(trim = FALSE, alpha = 0.6, scale = "width", color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(width = 0.08, size = 0.6, alpha = 0.5, color = "black") +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(x = NULL, y = ylab) +
    theme_fix()
  
  # 安全地添加显著性标注
  if (length(comps) > 0) {
    p <- p + geom_signif(
      comparisons      = comps,
      test             = "wilcox.test",
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
      step_increase    = 0.08,
      tip_length       = 0.02,
      color            = "black",
      textsize         = 3
    )
  }
  
  return(p)
}

# =============================================================================
# FIX A: Fig_Raincloud_AllScores
# =============================================================================
message("\n--- FIX A: Fig_Raincloud_AllScores ---")

rc_plots <- list()
for (sc in SCORE_COLS_FIX) {
  p <- tryCatch(make_viobox(sc), error = function(e) {
    message("  ERROR in ", sc, ": ", e$message); NULL
  })
  if (!is.null(p)) rc_plots[[sc]] <- p
}

message(sprintf("Successfully built %d / %d score plots", 
                length(rc_plots), length(SCORE_COLS_FIX)))

if (length(rc_plots) > 0) {
  nc <- 4
  nr <- ceiling(length(rc_plots) / nc)
  combined <- wrap_plots(rc_plots, ncol = nc)
  save_pdf_fix(combined, "Fig_Raincloud_AllScores",
               width = nc * 3.5, height = nr * 3.8)
}

# =============================================================================
# FIX B: Fig_C2_Fibrosis_Score
# =============================================================================
message("\n--- FIX B: Fig_C2_Fibrosis_Score ---")

if ("Fibrosis_Score" %in% colnames(meta_step2)) {
  p_fib <- tryCatch(
    make_viobox("Fibrosis_Score", "Fibrosis Score (ssGSEA)"),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  if (!is.null(p_fib)) save_pdf_fix(p_fib, "Fig_C2_Fibrosis_Score", width = 4.5, height = 4.5)
} else {
  message("  Column 'Fibrosis_Score' not found in meta_step2")
}

# =============================================================================
# FIX C: Fig_C3_Chemokine_Program_Score
# =============================================================================
message("\n--- FIX C: Fig_C3_Chemokine_Program_Score ---")

if ("Chemokine_Program_Score" %in% colnames(meta_step2)) {
  p_chemo <- tryCatch(
    make_viobox("Chemokine_Program_Score", "Chemokine Program Score (ssGSEA)"),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  if (!is.null(p_chemo)) save_pdf_fix(p_chemo, "Fig_C3_Chemokine_Program_Score", width = 4.5, height = 4.5)
} else {
  message("  Column 'Chemokine_Program_Score' not found in meta_step2")
}

# =============================================================================
# FIX D: Fig_Panel_ESTIMATE
# =============================================================================
message("\n--- FIX D: Fig_Panel_ESTIMATE ---")

est_cols <- intersect(c("Immune_Score","Stromal_Score","ESTIMATE_Score","Tumor_Purity"),
                      colnames(meta_step2))
message("  ESTIMATE columns found: ", paste(est_cols, collapse = ", "))

if (length(est_cols) > 0) {
  est_plots <- list()
  for (sc in est_cols) {
    p <- tryCatch(make_viobox(sc), error = function(e) {
      message("  ERROR in ", sc, ": ", e$message); NULL
    })
    if (!is.null(p)) est_plots[[sc]] <- p
  }
  message(sprintf("  Built %d / %d ESTIMATE plots", length(est_plots), length(est_cols)))
  
  if (length(est_plots) > 0) {
    p_est <- wrap_plots(est_plots, ncol = length(est_plots)) +
      plot_annotation(
        title = "ESTIMATE Scores",
        theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
      )
    save_pdf_fix(p_est, "Fig_Panel_ESTIMATE",
                 width = length(est_plots) * 3.5, height = 4.5)
  }
}

# =============================================================================
# FIX E: Fig_Panel_FunctionalScores
# =============================================================================
message("\n--- FIX E: Fig_Panel_FunctionalScores ---")

func_cols <- intersect(
  c("Proliferation_Score","ECM_Score","Myelination_Score","TIS_Score",
    "BCell_Score","Checkpoint_Score","Angiogenesis_Score","CYT_Score",
    "HLA_ClassI_Score","Fibrosis_Score"),
  colnames(meta_step2)
)
message("  Functional columns found: ", paste(func_cols, collapse = ", "))

if (length(func_cols) > 0) {
  func_plots <- list()
  for (sc in func_cols) {
    p <- tryCatch(make_viobox(sc), error = function(e) {
      message("  ERROR in ", sc, ": ", e$message); NULL
    })
    if (!is.null(p)) func_plots[[sc]] <- p
  }
  message(sprintf("  Built %d / %d functional plots", length(func_plots), length(func_cols)))
  
  if (length(func_plots) > 0) {
    nc_f <- 4
    nr_f <- ceiling(length(func_plots) / nc_f)
    p_func <- wrap_plots(func_plots, ncol = nc_f) +
      plot_annotation(
        title = "Subtype-Specific Functional Scores",
        theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
      )
    save_pdf_fix(p_func, "Fig_Panel_FunctionalScores",
                 width = 14, height = nr_f * 4)
  }
}

# =============================================================================
# FIX F: 相关性散点图 p 值格式化
# =============================================================================
message("\n--- FIX F: Correlation plots with proper p-value formatting ---")

HAS_GGRASTR_FIX <- requireNamespace("ggrastr", quietly = TRUE)

# F1: CYT vs Checkpoint
if (all(c("CYT_Score", "Checkpoint_Score") %in% colnames(meta_step2))) {
  
  cor_df <- meta_step2 %>%
    filter(!is.na(CYT_Score), !is.na(Checkpoint_Score))
  message("  CYT vs Checkpoint: n = ", nrow(cor_df))
  
  if (nrow(cor_df) >= 10) {
    cor_all <- cor.test(cor_df$CYT_Score, cor_df$Checkpoint_Score,
                        method = "spearman", exact = FALSE)
    
    p_cyt_cp <- ggplot(cor_df, aes(CYT_Score, Checkpoint_Score,
                                   color = Subtype, fill = Subtype)) +
      {
        if (HAS_GGRASTR_FIX)
          ggrastr::rasterise(geom_point(size = 2.5, alpha = 0.75, shape = 21,
                                        color = "white", stroke = 0.4), dpi = 300)
        else
          geom_point(size = 2.5, alpha = 0.75, shape = 21,
                     color = "white", stroke = 0.4)
      } +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("Spearman rho = %.3f\n%s",
                               cor_all$estimate, format_pval(cor_all$p.value)),
               hjust = 1.1, vjust = 1.5, size = 3.5, color = "black") +
      scale_color_manual(values = SUBTYPE_COLORS) +
      scale_fill_manual(values = SUBTYPE_COLORS) +
      labs(x = "CYT Score (Cytolytic Activity)", y = "Checkpoint Score",
           title = "Adaptive Immune Resistance: CYT vs Checkpoint") +
      theme_fix() + theme(legend.position = "right")
    
    save_pdf_fix(p_cyt_cp, "Fig_Adaptive_Resistance_CYT_vs_Checkpoint",
                 width = 6.5, height = 5)
    
    # Purity-corrected residual
    if ("Tumor_Purity" %in% colnames(cor_df)) {
      lm_cyt <- lm(CYT_Score ~ Tumor_Purity + Subtype, data = cor_df)
      lm_cp  <- lm(Checkpoint_Score ~ Tumor_Purity + Subtype, data = cor_df)
      cor_df$CYT_resid <- residuals(lm_cyt)
      cor_df$CP_resid  <- residuals(lm_cp)
      
      cor_resid <- cor.test(cor_df$CYT_resid, cor_df$CP_resid,
                            method = "spearman", exact = FALSE)
      
      p_resid <- ggplot(cor_df, aes(CYT_resid, CP_resid, color = Subtype)) +
        {
          if (HAS_GGRASTR_FIX)
            ggrastr::rasterise(geom_point(size = 2.5, alpha = 0.8), dpi = 300)
          else
            geom_point(size = 2.5, alpha = 0.8)
        } +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                    linewidth = 0.9, color = "gray30", alpha = 0.15) +
        annotate("text", x = Inf, y = Inf,
                 label = sprintf("Spearman rho = %.3f\n%s",
                                 cor_resid$estimate, format_pval(cor_resid$p.value)),
                 hjust = 1.1, vjust = 1.5, size = 3.5) +
        scale_color_manual(values = SUBTYPE_COLORS) +
        labs(x = "CYT Score (residual, purity+subtype adjusted)",
             y = "Checkpoint Score (residual)",
             title = "Purity-Corrected Adaptive Resistance") +
        theme_fix() + theme(legend.position = "right")
      
      save_pdf_fix(p_resid, "Fig_ResidualCorr_CYT_vs_Checkpoint",
                   width = 6.5, height = 5)
    }
  }
}

# F2: ImmuneScore vs TumorPurity
if (all(c("Immune_Score", "Tumor_Purity") %in% colnames(meta_step2))) {
  
  pur_df <- meta_step2 %>% filter(!is.na(Tumor_Purity), !is.na(Immune_Score))
  message("  ImmuneScore vs TumorPurity: n = ", nrow(pur_df))
  
  if (nrow(pur_df) >= 10) {
    cor_ov <- cor.test(pur_df$Tumor_Purity, pur_df$Immune_Score,
                       method = "spearman", exact = FALSE)
    
    p_pur <- ggplot(pur_df, aes(Tumor_Purity, Immune_Score,
                                color = Subtype, fill = Subtype)) +
      {
        if (HAS_GGRASTR_FIX)
          ggrastr::rasterise(geom_point(size = 2, alpha = 0.7, shape = 21,
                                        color = "white", stroke = 0.3), dpi = 300)
        else
          geom_point(size = 2, alpha = 0.7, shape = 21,
                     color = "white", stroke = 0.3)
      } +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("Spearman rho = %.3f\n%s",
                               cor_ov$estimate, format_pval(cor_ov$p.value)),
               hjust = 1.1, vjust = 1.5, size = 3.5) +
      scale_fill_manual(values = SUBTYPE_COLORS) +
      scale_color_manual(values = SUBTYPE_COLORS) +
      labs(x = "ESTIMATE Tumor Purity", y = "ESTIMATE Immune Score",
           title = "Tumor Purity vs Immune Infiltration") +
      theme_fix() + theme(legend.position = "right")
    
    save_pdf_fix(p_pur, "Fig_ImmuneScore_vs_TumorPurity", width = 6.5, height = 5)
    
    # Faceted
    p_pur_facet <- p_pur +
      facet_wrap(~ Subtype, nrow = 1) +
      theme(legend.position = "none")
    save_pdf_fix(p_pur_facet, "Fig_ImmuneScore_vs_TumorPurity_BySubtype",
                 width = 10, height = 4.5)
  }
}

message("\n========== SELF-CONTAINED FIX COMPLETE ==========\n")

# =============================================================================
# 高阶可视化：策略 1 & 策略 2 (全景热图与雷达图)
# =============================================================================

message("\n=== 开始绘制高阶 SCI 组图：全景热图与雷达图 ===")

# --- 0. 加载必要的极客画图包 ---
# ComplexHeatmap 是画顶刊热图的绝对主力，fmsb 用来画雷达图
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!require("circlize")) install.packages("circlize")
if (!require("fmsb")) install.packages("fmsb")

library(ComplexHeatmap)
library(circlize)
library(fmsb)
library(dplyr)
library(tidyr)

# 强制统一颜色字典（防止前面环境变量丢失）
my_subtype_colors <- c("C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35")

# --- 1. 整理 15 个评分，并严格按照 4 个生物学模块排序 ---
macro_cols    <- c("Tumor_Purity", "ESTIMATE_Score", "Immune_Score", "Stromal_Score")
immune_cols   <- c("CYT_Score", "TIS_Score", "BCell_Score", "HLA_ClassI_Score", "Chemokine_Program_Score", "Checkpoint_Score")
stromal_cols  <- c("ECM_Score", "Fibrosis_Score", "Angiogenesis_Score")
intrinsic_cols<- c("Proliferation_Score", "Myelination_Score")

all_15_cols <- c(macro_cols, immune_cols, stromal_cols, intrinsic_cols)

# 定义这 15 个变量所属的模块名称（用于热图顶部分块）
col_split_factor <- factor(
  c(rep("Macro-TME", length(macro_cols)),
    rep("Immune", length(immune_cols)),
    rep("Stromal", length(stromal_cols)),
    rep("Intrinsic", length(intrinsic_cols))),
  levels = c("Macro-TME", "Immune", "Stromal", "Intrinsic")
)

# 清洗数据：提取相关列并去除含有 NA 的样本
tme_data <- meta_step2 %>% 
  select(SampleID, Subtype, all_of(all_15_cols)) %>% 
  drop_na() %>% 
  arrange(Subtype) # 让样本按照 C1, C2, C3 排序，图面更整齐

# 提取纯数字矩阵
score_mat <- as.matrix(tme_data[, all_15_cols])
rownames(score_mat) <- tme_data$SampleID


# =============================================================================
# 🚀 策略 1: 模块化全景热图 (ComplexHeatmap)
# =============================================================================
message("\n--- 正在生成策略1：模块化全景热图 ---")

# 1. 计算 Z-score 标准化（按列计算，使得不同评分可以放在同一种颜色刻度下）
score_mat_z <- scale(score_mat)
# 防止出现极端异常值导致颜色失真，将 Z-score 截断在 -3 到 3 之间
score_mat_z[score_mat_z > 3] <- 3
score_mat_z[score_mat_z < -3] <- -3

# 2. 定义热图的红白蓝连续配色 (遵循您的整体配色调性)
col_fun <- colorRamp2(c(-2, 0, 2), c("#471669", "#238C8D", "#FCE828"))

# 3. 创建左侧样本分组的颜色条注释
row_anno <- rowAnnotation(
  Subtype = tme_data$Subtype,
  col = list(Subtype = my_subtype_colors),
  show_annotation_name = FALSE,
  annotation_width = unit(5, "mm")
)

# 4. 绘制终极 ComplexHeatmap
ht <- Heatmap(
  score_mat_z,
  name = "Z-score",               # 图例名称
  col = col_fun,                  # 颜色映射
  cluster_rows = FALSE,           # 不对行(样本)聚类，因为我们已经按C1/C2/C3排好序了
  cluster_columns = FALSE,        # 不对列(评分)聚类，严格遵循我们定义的4个生物学模块顺序
  row_split = tme_data$Subtype,   # 行按 C1, C2, C3 物理断开
  column_split = col_split_factor,# 列按 4 个生物学模块物理断开
  left_annotation = row_anno,     # 添加左侧颜色条
  show_row_names = FALSE,         # 隐藏样本名，保持图面干净
  show_column_names = TRUE,
  column_names_rot = 45,          # 底部文字倾斜45度
  column_names_gp = gpar(fontsize = 10),
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  border = TRUE                   # 给每个区块加上边框，高级感拉满
)

# 5. 保存 PDF
pdf("Fig_Strategy1_Panorama_Heatmap.pdf", width = 7, height =5)
draw(ht, column_title = "Comprehensive TME Profiling by Subtype", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()



# 
# 
# # =============================================================================
# # 🚀 策略 1: 模块化全景热图 (ComplexHeatmap - 独立方格版)
# # =============================================================================
# message("\n--- 正在生成策略1：模块化全景热图 (方格样式) ---")
# 
# # 1. 计算 Z-score 标准化
# score_mat_z <- scale(score_mat)
# score_mat_z[score_mat_z > 3] <- 3
# score_mat_z[score_mat_z < -3] <- -3
# 
# # 2. 定义热图的配色 (延续您的冷暖撞色)
# col_fun <- colorRamp2(c(-2, 0, 2), c("#471669", "#238C8D", "#FCE828"))
# 
# # 3. 创建左侧样本分组的颜色条注释
# row_anno <- rowAnnotation(
#   Subtype = tme_data$Subtype,
#   col = list(Subtype = my_subtype_colors),
#   show_annotation_name = FALSE,
#   annotation_width = unit(5, "mm")
# )
# 
# # 4. 绘制终极 ComplexHeatmap
# ht <- Heatmap(
#   score_mat_z,
#   name = "Z-score",               
#   col = col_fun,                  
#   cluster_rows = FALSE,           
#   cluster_columns = FALSE,        
#   row_split = tme_data$Subtype,   
#   column_split = col_split_factor,
#   left_annotation = row_anno,     
#   show_row_names = FALSE,         
#   show_column_names = TRUE,
#   column_names_rot = 45,          
#   column_names_gp = gpar(fontsize = 10),
#   row_title_gp = gpar(fontsize = 12, fontface = "bold"),
#   column_title_gp = gpar(fontsize = 12, fontface = "bold"),
#   border = TRUE,                  
#   
#   # ✨ 核心魔法：为每个小格子加上白色边框，制造“独立方块”的视觉效果
#   # col 控制缝隙颜色，lwd 控制缝隙的宽度。您可以调整 lwd（比如 1, 1.5, 2）来改变缝隙大小
#   rect_gp = gpar(col = "white", lwd = 1.5) 
# )
# 
# # 5. 保存 PDF
# # ⚠️ 极度预警：因为有约 200 个样本加上了白缝，高度为 5 绝对会被缝隙吃掉颜色！
# # 建议将 height 调高至 12 或以上。您可以先用 height=12 测试，如果格子依然太扁，就继续加高。
# pdf("Fig_Strategy1_Panorama_Heatmap_Tiles.pdf", width = 7, height = 12) 
# 
# draw(ht, column_title = "Comprehensive TME Profiling by Subtype", 
#      column_title_gp = gpar(fontsize = 16, fontface = "bold"))
# dev.off()
# =============================================================================
# 🚀 策略 1: 模块化全景热图 (横向方格宽屏版)
# =============================================================================
message("\n--- 正在生成策略1：横向模块化全景热图 ---")

# 1. 核心操作：矩阵转置 (Transpose)
# 让行变成 15 个功能评分，列变成您的约 200 个样本
score_mat_z_t <- t(score_mat_z)

# 配色保持您原来的高级撞色不变
col_fun <- colorRamp2(c(-2, 0, 2), c("#471669", "#238C8D", "#FCE828"))

# 2. 注释条位置转换：从左侧 (rowAnnotation) 改为顶部 (HeatmapAnnotation)
top_anno <- HeatmapAnnotation(
  Subtype = tme_data$Subtype,
  col = list(Subtype = my_subtype_colors),
  show_annotation_name = FALSE,
  simple_anno_size = unit(5, "mm") # 控制顶部颜色条的粗细
)

# 3. 绘制横向热图
ht_horizontal <- Heatmap(
  score_mat_z_t,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,           
  cluster_columns = FALSE,        
  
  # ✨ 参数乾坤大挪移
  row_split = col_split_factor,    # 现在的行，用 4 个生物学模块来切分
  column_split = tme_data$Subtype, # 现在的列，用 C1, C2, C3 来切分
  top_annotation = top_anno,       # 分组颜色条顶置
  
  # ✨ 文字标签对齐
  show_column_names = FALSE,       # 隐藏底部密密麻麻的样本名
  show_row_names = TRUE,           # 显示那 15 个评分的名字
  row_names_side = "left",         # 评分名放在左侧，紧挨着色块，最符合阅读习惯
  row_names_gp = gpar(fontsize = 10),
  
  # ✨ 让外侧的分类大标题 (Macro-TME等) 横向显示，提高可读性
  row_title_rot = 0,               
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  border = TRUE,
  # 继续保留高级的白色方格马赛克效果 (适当微调线宽防止太挤)
  rect_gp = gpar(col = "white", lwd = 1.2) 
)

# 4. 保存为宽幅 PDF (注意 width 和 height 已经对调)
# 宽度 14 英寸足够容纳 200 个样本的列，高度 6 或 7 足够容纳 15 行评分
pdf("Fig_Strategy1_Panorama_Heatmap_Horizontal.pdf", width = 12, height =5)

draw(ht_horizontal, 
     column_title = "Comprehensive TME Profiling by Subtype", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold"))

dev.off()










# =============================================================================
# 🚀 策略 2: 微环境特征雷达图 (Radar Chart)
# =============================================================================
message("\n--- 正在生成策略2：微环境特征雷达图 ---")

# 1. 计算每个 Subtype 在 15 个指标上的中位数
radar_meds <- tme_data %>%
  group_by(Subtype) %>%
  summarise(across(all_of(all_15_cols), ~median(.x, na.rm = TRUE))) %>%
  column_to_rownames("Subtype")

# 2. 0-1 标准化：雷达图的各个轴尺度必须一致
# 将每一列的最大值缩放为 1，最小值缩放为 0
radar_sc <- as.data.frame(
  apply(radar_meds, 2, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9))
)

# 3. fmsb 画雷达图需要特定格式：第一行必须是上限(1)，第二行必须是下限(0)
radar_fmsb <- rbind(rep(1, ncol(radar_sc)), rep(0, ncol(radar_sc)), radar_sc)

# 清理列名，去掉烦人的 "_Score" 字符，让雷达图文字更美观
colnames(radar_fmsb) <- gsub("_Score", "", colnames(radar_fmsb))
colnames(radar_fmsb) <- gsub("_", " ", colnames(radar_fmsb))

# 4. 设置带有透明度的雷达图填充颜色 (使用 base R 的 adjustcolor)
fill_colors <- adjustcolor(my_subtype_colors[c("C1", "C2", "C3")], alpha.f = 0.2)
line_colors <- my_subtype_colors[c("C1", "C2", "C3")]

# 5. 绘制并保存雷达图
pdf("Fig_Strategy2_TME_RadarChart.pdf", width = 8, height = 8)
# 设置画布边缘留白，防止标签被切掉
par(mar = c(2, 2, 4, 2)) 

radarchart(
  radar_fmsb,
  axistype = 1,
  pcol = line_colors,          # 边框颜色
  pfcol = fill_colors,         # 填充颜色
  plwd = 2,                    # 线条粗细
  cglcol = "grey70",           # 背景网格颜色
  cglty = 1,                   # 背景网格线型 (实线)
  cglwd = 0.8,                 # 背景网格线粗细
  axislabcol = "grey30",       # 轴标签颜色
  caxislabels = seq(0, 1, 0.25), # 轴刻度 (0, 0.25, 0.5, 0.75, 1)
  vlcex = 0.9,                 # 顶点文字大小
  title = "TME Multi-Dimensional Fingerprint"
)

# 添加图例
legend(
  x = "topright", 
  legend = c("C1 (Proliferative)", "C2 (Stromal)", "C3 (Immune)"),
  col = line_colors, 
  lty = 1, 
  lwd = 3,
  bty = "n", 
  cex = 1.1
)
dev.off()

message("\n=== 恭喜！全景热图与雷达图已成功保存到工作目录！ ===")



















# === 加载R包 ===
library(CIBERSORT)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(vioplot)
library(RColorBrewer)
library(readxl)
# expr <- read.table("clean_tpm.csv", header =TRUE, sep ="\t", row.names =1, check.names =FALSE)
# 
# # === 运行CIBERSORT（需准备好LM22特征矩阵与CIBERSORT.R脚本）===
# source("CIBERSORT.R")
 # res_cibersort <-CIBERSORT("LM22.txt","TPM.txt",
 #                                      perm =1000, QN =TRUE)
# write.csv(res_cibersort,"CIBERSORT_Result.csv")
# 1. 安装并加载必要的包 (如果还没安装，请取消注释第一行)
# install.packages("readxl")
# 
# # === 1. 加载 CIBERSORT 脚本 ===
# 
# source("CIBERSORT.R") 
# source("doPerm.R")    
# source("CoreAlg.R")
# source("utils-pipe.R")
# source("untils.R")

expr_data <- read.csv("clean_tpm name.csv", header = TRUE, check.names = FALSE)
colnames(expr_data)[1] <- "Gene" # 统一第一列名字

# === 3. 数据清洗：去除重复基因并设置行名 ===
expr_clean <- expr_data[!duplicated(expr_data[, 1]), ] # 去除第一列有重复的行
rownames(expr_clean) <- expr_clean[, 1]                # 把第一列变成行名
expr_clean <- expr_clean[, -1]                         # 删掉已经变成行名的第一列

# === 4. 导出为 CIBERSORT 专用的 txt 文件 ===
write.table(expr_clean, file = "TPM_for_Cibersort.txt", sep = "\t", quote = FALSE, col.names = NA)

# 
# # === 1. 安装并加载并行计算缺失的包 ===
# if (!require("parallelly")) install.packages("parallelly")
# library(parallelly)
# 
# # 为了防止接下来报 purrr 的错，顺手也把它加载上
# if (!require("purrr")) install.packages("purrr")
# library(purrr)
# # === 1. 一次性补齐作者使用的现代 R 包全家桶 ===
# # 安装并加载 furrr (解决 current future_map 报错)
# if (!require("furrr")) install.packages("furrr")
# library(furrr)
# 
# # 预防性加载 tidyverse 核心包（防止一会儿又报 %>% 或 mutate 找不到）
# if (!require("dplyr")) install.packages("dplyr")
# if (!require("tidyr")) install.packages("tidyr")
# library(dplyr)
# library(tidyr)
# # === 2. 激动人心的最后一次运行！ ===
# res_cibersort <- cibersort(
#   "LM22.txt", 
#   "TPM_for_Cibersort.txt", 
#   perm = 1000, 
#   QN = FALSE  
# )
# 
# # === 3. 保存结果 ===
# write.csv(res_cibersort, "CIBERSORT_Result.csv")
install.packages("devtools")
devtools::install_github("bmbolstad/preprocessCore")
devtools::install_github("Moonerss/CIBERSORT")
# ==========================================
# 第二步：像使用普通 R 包一样，加载并运行（以后每次只需跑这段）
# ==========================================
library(CIBERSORT)

# 激动人心的决战时刻：替换成咱们刚才千辛万苦洗干净的本地数据！
res_cibersort <- cibersort(
  sig_matrix = "LM22.txt",             # 您的 LM22 特征矩阵
  mixture_file = "TPM_for_Cibersort.txt",  # 您刚才去重处理好的纯净 TPM 数据
  perm = 1000, 
  QN = FALSE                           # 再次强调，TPM 数据务必设为 FALSE
)

# 导出胜利的果实
write.csv(res_cibersort, "CIBERSORT_Result.csv")

# 1. 读取完整的 CIBERSORT 原始结果（包含 P-value 列）
# 使用 check.names = FALSE 防止 R 把 "P-value" 自动改成 "P.value"
res_raw <- read.csv("Figures/step8_Immune/CIBERSORT_Result.csv", row.names = 1, check.names = FALSE)

# 2. 智能抓取 P 值所在的列名（无论是 P-value 还是 P.value 都能抓住）
pval_col <- grep("P-value|P.value", colnames(res_raw), value = TRUE, ignore.case = TRUE)[1]

# 3. 计算 P < 0.05 的样本数和总样本数
n_retained <- sum(res_raw[[pval_col]] < 0.05, na.rm = TRUE)
n_total <- nrow(res_raw)

# 4. 打印最终文章需要的标准结果格式
cat(sprintf("✅ CIBERSORT: number of samples retained after P < 0.05 filter (%d/%d).\n", n_retained, n_total))

getwd()
setwd("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/step8_Immune/")
# === 加载所需 R 包 ===
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(vioplot)

# === 1. 读取并清理数据 ===
res1 <- read.csv("CIBERSORT_Result.csv", row.names = 1)
# 严格提取前22列免疫细胞比例（剔除P值和RMSE列）
res1 <- res1[, 1:22] 

# ==========================================
# 📊 图1：免疫细胞比例堆积条形图 (Stacked Barplot)
# ==========================================
# # 将宽矩阵转换为适合 ggplot 绘图的长矩阵
# res_melt <- melt(as.matrix(res1))
# colnames(res_melt) <- c("Sample", "CellType", "Proportion")
# 
# p1 <- ggplot(res_melt, aes(x = Sample, y = Proportion, fill = CellType)) +
#   geom_bar(stat = "identity", width = 1) +  # width=1 可以让柱子紧凑贴合
#   theme_bw() +
#   theme(
#     axis.text.x = element_blank(),   # 样本太多时隐藏 x 轴名字
#     axis.ticks.x = element_blank(),
#     legend.position = "right",
#     legend.key.size = unit(0.4, "cm")
#   ) +
#   labs(title = "Immune Cell Composition", x = "Samples", y = "Relative Fraction")
# 
# ggsave("CIBERSORT_Barplot.pdf", plot = p1, width = 12, height = 6)
# === 1. 读取并合并数据 ===
# 假设 res1 依然是您的前 22 列细胞比例矩阵
group <- read.csv("group.csv") 
res1$Sample <- rownames(res1)
res_merged <- merge(res1, group, by = "Sample")

# === 2. 宽数据转长数据 ===
library(reshape2)
res_melt <- melt(res_merged, 
                 id.vars = c("Sample", "Group"), 
                 variable.name = "CellType", 
                 value.name = "Proportion")

# === 3. 【重点】自定义 22 种免疫细胞的颜色 ===
# 这里为您精心挑选了 22 种区分度极高、柔和且不刺眼的莫兰迪/马卡龙色系。
# 如果您不喜欢某一种，直接把对应的 "#XXXXXX" 替换成您想要的色号即可。
my_22_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#E31A1C", "#BC80BD", "#CCEBC5", "#FFED6F",
  "#1F78B4", "#A9DAD1", "#D9D9D9","#FB9A99" , "#FCCDE5", "#6A3D9A", "#B15928",
  "#A6CEE3", "#B2DF8A", "#FF7F00", "#FDBF6F"
)
# strip_colors <- c("#00A087", "#4DBBD5", "#E64B35") 

# === 1. 只加载最基础、最稳妥的包 ===
library(ggplot2)
library(grid)      # R语言自带的底层画图引擎，不用安装
library(reshape2)

# (假设前面的 res_melt 和 my_22_colors 依然存在)

# === 2. 画出最纯净的原生 ggplot (白底框) ===
p1 <- ggplot(res_melt, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1, color = NA) + 
  
  # 换回原生绝对不报错的 facet_grid
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  
  scale_fill_manual(values = my_22_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    panel.spacing = unit(0.2, "lines") 
  ) +
  labs(title = "Immune Cell Composition by Group", 
       x = "Grouped Samples", 
       y = "Relative Percent",
       fill = "Immune Cell Type")

# === 3. 【神级操作：进入底层修改颜色】 ===
# 定义您的 3 个组的颜色
# strip_colors <- c("#C1FBB4", "#B4DDFB", "#FBC1B4")
strip_colors <- c("#00A087", "#4DBBD5", "#E64B35")
# 将图转换为底层图形对象 (Grob)
g <- ggplotGrob(p1)

# 自动寻找图中所有属于 "strip-t" (顶部标签框) 的位置
strips <- which(grepl('strip-t', g$layout$name))

# 循环给这 3 个框强行注入颜色
for (i in seq_along(strips)) {
  # 锁定背景矩形块并修改填充色 (fill)
  rect_index <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[rect_index]]$gp$fill <- strip_colors[i]
}

# === 4. 保存这个修改后的底层对象 ===
# ggsave 非常聪明，它可以直接保存 Grob 对象
ggsave("CIBERSORT_Grouped_ColoredStrips_Ultimate.pdf", plot = g, width = 10, height = 5)





# 📊 图2：免疫细胞间的相关性热图 (Heatmap)
# ==========================================
# # 计算细胞之间的斯皮尔曼相关系数
# cor_matrix <- cor(res1, method = "spearman")
# 
# pdf("CIBERSORT_Correlation.pdf", width = 8, height = 8)
# pheatmap(cor_matrix,
#          color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
#          display_numbers = FALSE, # 如果格子太小，不显示具体数字更美观
#          fontsize_row = 10,
#          fontsize_col = 10,
#          main = "Immune Cell Correlation")
# dev.off()
# ==========================================
# 📊 图3：免疫细胞间的相关性热图 (修复版)
# ==========================================

# 1. 计算每一列（每种细胞）的方差，保留方差大于 0 的列
# apply(res1, 2, var) 会计算每一列的方差
res1_filtered <- res1[, apply(res1, 2, var) > 0]

# 2. 使用过滤后的干净数据计算相关性
# 这次就不会再报“标准差为零”的警告了
cor_matrix <- cor(res1_filtered, method = "spearman")

# 3. 绘制热图
pdf("CIBERSORT_Correlation.pdf", width = 8, height = 8)
pheatmap(cor_matrix,
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         display_numbers = FALSE, 
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Immune Cell Correlation")
dev.off()

# # ==========================================
# # 📊 图3：组间免疫细胞差异小提琴图 (Violin Plot)
# # ==========================================
# # 读取分组信息 (需包含列名：Sample 和 Group)
# group <- read.csv("group.csv")
# 
# # 将分组信息与细胞比例合并
# res1$Sample <- rownames(res1)
# res1_grouped <- merge(res1, group, by = "Sample")
# 
# pdf("CIBERSORT_Vioplot.pdf", width = 14, height = 10)
# # 设置画布布局：4行6列，刚好放下 22 个细胞的图
# par(mfrow = c(4, 6), mar = c(3, 3, 2, 1))
# 
# # 循环绘制 22 种细胞
# for (i in 1:22) {
#   cell_name <- colnames(res1)[i]
# 
#   # 提取两组的具体数值（这里以 Large 和 Small 分组为例）
#   # 请务必将 "Large" 和 "Small" 替换为您 group.csv 里真实的组名
#   group_A <- res1_grouped[[cell_name]][res1_grouped$Group == "Large"]
#   group_B <- res1_grouped[[cell_name]][res1_grouped$Group == "Small"]
# 
#   # 只有当两组都有数据时才画图，防止因数据缺失报错
#   if(length(group_A) > 0 & length(group_B) > 0) {
#     vioplot(group_A, group_B,
#             names = c("Large", "Small"),
#             col = c("#80B1D3", "#FB8072"),
#             main = cell_name,
#             cex.main = 0.8) # 缩小标题字号防止重叠
#   }
# }
# dev.off()
# res2<- res1 / rowSums(res1)
# res2_melt<- melt(as.matrix(res2))
# colnames(res2_melt) <- c("Sample","CellType","Fraction")
# p2<- ggplot(res2_melt, aes(x = Sample, y = Fraction, fill = CellType)) +geom_bar(stat ="identity") +theme_minimal() +labs(title ="Normalized Immune Cell Fractions")
# ggsave("CIBERSORT2.pdf", plot = p2, width =10, height =5)
# 
# 
# cor_matrix<- cor(res1, method ="spearman")
# pdf("cor.pdf", width =8, height =8)
# pheatmap(cor_matrix,color= colorRampPalette(brewer.pal(9,"RdBu"))(100),main="Immune Cell Correlation Heatmap")
# dev.off()

# 
# # === 1. 严格过滤掉标准差为 0 的细胞列 ===
# # 使用 sd() 过滤比 var() 更稳妥
# res1_filtered <- res1[, apply(res1, 2, sd) > 0]
# 
# # === 2. 计算相关性矩阵 ===
# cor_matrix <- cor(res1_filtered, method = "spearman")
# 
# # === 3. 【关键补丁】将矩阵中所有残留的 NA 替换为 0 ===
# cor_matrix[is.na(cor_matrix)] <- 0
# 
# # === 4. 放心大胆地画图 ===
# pdf("cor.pdf", width = 8, height = 8)
# pheatmap(cor_matrix,
#          color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
#          display_numbers = FALSE, # 隐藏数字让热图更清爽
#          fontsize_row = 10,
#          fontsize_col = 10,
#          main = "Immune Cell Correlation Heatmap")
# dev.off()



# === 1. 加载学术绘图神器 ===
library(ggplot2)
library(reshape2)
library(ggpubr)

# === 2. 【缺失的关键步】读取 CIBERSORTx 网页版的结果文件 ===
# 请把下方引号里的名字，替换成您真实下载的那个 csv 文件名！
res_raw <- read.csv("CIBERSORT_Result.csv", row.names = 1)

# 提取前 22 列真正的免疫细胞比例（剔除最后面的 P-value、Correlation、RMSE 列）
res1 <- res_raw[, 1:22]

# === 3. 读取分组信息并合并 ===
# 确保您的 group.csv 和 R 工作路径在同一个文件夹下
group <- read.csv("group.csv") 

# 把行名（样本名）变成一列，方便与 group 表合并
res1$Sample <- rownames(res1)
res_merged <- merge(res1, group, by = "Sample")

# === 4. 数据宽转长（ggplot 的标准格式要求） ===
res_melt <- melt(res_merged, 
                 id.vars = c("Sample", "Group"), 
                 variable.name = "CellType", 
                 value.name = "Proportion")

# # === 5. 绘制高颜值多组并排箱线图并自动计算 P 值 ===
# p3 <- ggplot(res_melt, aes(x = CellType, y = Proportion, fill = Group)) +
#   geom_boxplot(outlier.shape = 21, outlier.size = 0.5, alpha = 0.8) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), # 细胞名倾斜防重叠
#     legend.position = "top",                                            # 图例放顶部更美观
#     panel.grid.major.x = element_blank()                                # 去除多余网格线
#   ) +
#   labs(title = "Immune Infiltration across Three Groups",
#        x = "Immune Cell Types", 
#        y = "Relative Fraction",
#        fill = "Tumor Group") +
#   # 自动执行 Kruskal-Wallis 多组检验，并在有显著差异的细胞上方打星号
#   stat_compare_means(aes(group = Group), 
#                      label = "p.signif",    # 只显示星号 (*, **, ***)
#                      method = "kruskal.test", 
#                      hide.ns = TRUE)        # 隐藏没有显著差异 (ns) 的标记
# print(p3)
# # === 6. 输出宽幅 PDF ===
# # 3组并排需要更宽的画布，width 设为 16
# ggsave("CIBERSORT_3Groups_Boxplot.pdf", plot = p3, width = 10, height = 5)
# === 自定义您的分组颜色 ===
# 这里我为您预设了 3 个非常经典的 SCI 柔和配色（分别对应组1、组2、组3）
# 您可以随意替换双引号里的十六进制颜色码或英文颜色名
# my_colors <- c("#00A087", "#4DBBD5", "#E64B35") 
my_colors <- c("#7BCEC1", "#9BD9E8", "#F09C90") 
# === 绘制并应用自定义颜色 ===
p3 <- ggplot(res_melt, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 0.5, alpha = 0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
    legend.position = "top",                                            
    panel.grid.major.x = element_blank()                                
  ) +
  labs(title = "Immune Infiltration across Three Groups",
       x = "Immune Cell Types", 
       y = "Relative Fraction",
       fill = "Tumor Group") +
  stat_compare_means(aes(group = Group), 
                     label = "p.signif",    
                     method = "kruskal.test", 
                     hide.ns = TRUE) +
  # ✨ 重点在这里：将您定义的颜色向量传给这行代码
  scale_fill_manual(values = my_colors) 

# 输出图片
ggsave("CIBERSORT_3Groups_CustomColor.pdf", plot = p3, width = 8, height = 6)






# === 1. 加载画图包 ===
library(pheatmap)
# === 2. 【缺失的关键步】读取 CIBERSORTx 网页版的结果文件 ===
# 请把下方引号里的名字，替换成您真实下载的那个 csv 文件名！
res_raw <- read.csv("CIBERSORT_Result.csv", row.names = 1)

# # 提取前 22 列真正的免疫细胞比例（剔除最后面的 P-value、Correlation、RMSE 列）
 res1 <- res_raw[, 1:22]
# === 2. 提取数据（假设您的结果保存在 res1 中） ===
cibersort_result <- as.data.frame(res1[, 1:22])

# === 3. 读取并严格对齐分组信息 ===
group_info <- read.csv('group.csv', header = TRUE)
# 强制将列名重命名为标准格式，彻底杜绝大小写报错
colnames(group_info) <- c("Sample", "Group") 

# 找出结果矩阵和分组表里共有的样本，防止有些样本没有临床信息导致报错
common_samples <- intersect(rownames(cibersort_result), group_info$Sample)
group_info <- group_info[match(common_samples, group_info$Sample), ]

# 【美观核心】将样本按照 Group 分类重新排序
group_info <- group_info[order(group_info$Group), ]
# 让表达矩阵的行顺序也跟着改变，完美对齐
cibersort_result <- cibersort_result[group_info$Sample, ]

# === 4. 准备热图的高级注释条 ===
annotation_row <- data.frame(Group = group_info$Group)
rownames(annotation_row) <- group_info$Sample

# （可选）如果您想自定义右侧分组条的颜色，取消下面这几行的注释：
# ann_colors <- list(
#   Group = c(Group1 = "#E64B35", Group2 = "#4DBBD5", Group3 = "#00A087") # 替换为您的真实组名
# )

# === 1. 【核心修复】自动识别并剔除标准差为 0 的“死水”细胞 ===
# apply 函数会计算每一列的 sd（标准差），保留那些有数值波动的细胞
valid_cells <- apply(cibersort_result, 2, sd) > 0
cibersort_result_clean <- cibersort_result[, valid_cells]

# === 2. 重新开启 PDF 画布 ===
pdf("Immune_Cell_Abundance_Heatmap.pdf", width = 10, height = 12) 

# === 3. 用清洗后的干净数据画图 ===
pheatmap(as.matrix(cibersort_result_clean), 
         scale = "column",       
         cluster_rows = FALSE,   
         cluster_cols = TRUE,    
         color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100), 
         show_rownames = FALSE,  
         show_colnames = TRUE,   
         main = "Immune Cell Abundance Landscape",
         annotation_row = annotation_row, 
         annotation_names_row = FALSE,
         fontsize_col = 10,
         angle_col = "45"        
)

# === 4. 保存文件 ===
dev.off()












# =============================================================================
# SECTION 32: C3 immune recruitment and suppression
# 目标：
#   Step 1. 判断 C3 的 immune recruitment 是增强还是减弱
#   Step 2. 判断 C3 内部 recruitment 是否与 suppression / checkpoint 相关
# =============================================================================

message("\n=== SECTION 32: C3 immune recruitment and suppression ===")

# -------------------------------------------------------------------------
# 0. 基因集：尽量和你之前代码保持一致的 strict 版本
# -------------------------------------------------------------------------
IMMUNE_RECRUITMENT_STRICT <- c("CCL4","CCL5","CXCL9","CXCL10")
# IMMUNE_SUPPRESSION_STRICT <- c("CD274","CTLA4","HAVCR2","LAG3","TIGIT","PDCD1LG2","IDO1")
IMMUNE_SUPPRESSION_STRICT <- c("PDCD1","CD274","CTLA4","HAVCR2",
"LAG3","TIGIT","VSIR","CD96","SIGLEC7")
recruitment_present <- intersect(IMMUNE_RECRUITMENT_STRICT, rownames(clean_tpm))
suppression_present <- intersect(IMMUNE_SUPPRESSION_STRICT, rownames(clean_tpm))

message("Recruitment genes present: ", paste(recruitment_present, collapse = ", "))
message("Suppression genes present: ", paste(suppression_present, collapse = ", "))

if (length(recruitment_present) < 3) {
  stop("Too few recruitment genes found in clean_tpm (<3).")
}
if (length(suppression_present) < 3) {
  stop("Too few suppression genes found in clean_tpm (<3).")
}

# -------------------------------------------------------------------------
# 1. 计算 Recruitment / Suppression score
#    沿用你前面的 run_ssgsea 和 meta_step2
# -------------------------------------------------------------------------
RS_GENESETS <- list(
  Immune_Recruitment_Score = recruitment_present,
  Immune_Suppression_Score = suppression_present
)

rs_scores <- run_ssgsea(clean_tpm, RS_GENESETS)

meta_step2$Immune_Recruitment_Score <- as.numeric(
  rs_scores["Immune_Recruitment_Score", meta_step2$SampleID]
)
meta_step2$Immune_Suppression_Score <- as.numeric(
  rs_scores["Immune_Suppression_Score", meta_step2$SampleID]
)

rm(rs_scores)
gc()

# -------------------------------------------------------------------------
# 2. 辅助函数：p 值格式
# -------------------------------------------------------------------------
format_pval2 <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  if (p < 1e-4) return(sprintf("p = %.2e", p))
  sprintf("p = %.4f", p)
}

pair_comp <- list(c("C1","C2"), c("C1","C3"), c("C2","C3"))

# -------------------------------------------------------------------------
# STEP 1：先判断 C3 的 immune recruitment 是增强还是减弱
# -------------------------------------------------------------------------

# 2.1 summary table
recruitment_summary <- meta_step2 %>%
  dplyr::select(SampleID, Subtype, Immune_Recruitment_Score, Immune_Suppression_Score) %>%
  dplyr::group_by(Subtype) %>%
  dplyr::summarise(
    n = dplyr::n(),
    Recruitment_mean   = mean(Immune_Recruitment_Score, na.rm = TRUE),
    Recruitment_median = median(Immune_Recruitment_Score, na.rm = TRUE),
    Suppression_mean   = mean(Immune_Suppression_Score, na.rm = TRUE),
    Suppression_median = median(Immune_Suppression_Score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  recruitment_summary,
  file.path(OUTPUT_DIR, "C3_ImmuneRecruitment_Suppression_Summary.csv"),
  row.names = FALSE, quote = FALSE
)

# 2.2 recruitment 组间统计
rec_stat <- auto_stat_test(meta_step2, "Immune_Recruitment_Score", "Subtype")

rec_pw <- rec_stat$pairwise_df
rec_comp <- lapply(rec_pw$comparison, parse_comparison)

df_rec <- meta_step2 %>%
  dplyr::select(Subtype, Score = Immune_Recruitment_Score) %>%
  dplyr::filter(!is.na(Score))

p_rec <- ggplot(df_rec, aes(x = Subtype, y = Score, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.70, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.45, color = "black") +
  geom_signif(
    comparisons = rec_comp,
    map_signif_level = function(p) {
      if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
    },
    step_increase = 0.08,
    tip_length = 0.02,
    color = "black",
    textsize = 3
  ) +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  labs(
    x = NULL,
    y = "Immune Recruitment Score",
    title = "Immune Recruitment across VS Subtypes"
  ) +
  theme_pub() +
  theme(legend.position = "none")

save_plot(p_rec, "Fig_C3_ImmuneRecruitment_Violin", width = 5, height = 5)

# 2.3 suppression 组间统计（辅助图，建议保留）
sup_stat <- auto_stat_test(meta_step2, "Immune_Suppression_Score", "Subtype")

sup_pw <- sup_stat$pairwise_df
sup_comp <- lapply(sup_pw$comparison, parse_comparison)

df_sup <- meta_step2 %>%
  dplyr::select(Subtype, Score = Immune_Suppression_Score) %>%
  dplyr::filter(!is.na(Score))

p_sup <- ggplot(df_sup, aes(x = Subtype, y = Score, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.70, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.45, color = "black") +
  geom_signif(
    comparisons = sup_comp,
    map_signif_level = function(p) {
      if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
    },
    step_increase = 0.08,
    tip_length = 0.02,
    color = "black",
    textsize = 3
  ) +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  labs(
    x = NULL,
    y = "Immune Suppression Score",
    title = "Immune Suppression across VS Subtypes"
  ) +
  theme_pub() +
  theme(legend.position = "none")

save_plot(p_sup, "Fig_C3_ImmuneSuppression_Violin", width = 5, height = 5)

# -------------------------------------------------------------------------
# STEP 2：再看 C3 内部招募和抑制 / checkpoint 是否相关
# -------------------------------------------------------------------------

c3_df <- meta_step2 %>%
  dplyr::filter(Subtype == "C3") %>%
  dplyr::select(
    SampleID, Subtype,
    Immune_Recruitment_Score,
    Immune_Suppression_Score,
    dplyr::any_of(c("Checkpoint_Score", "Tumor_Purity"))
  ) %>%
  dplyr::filter(!is.na(Immune_Recruitment_Score), !is.na(Immune_Suppression_Score))

# 3.1 C3 内部：Recruitment vs Suppression
if (nrow(c3_df) >= 6) {
  cor_c3 <- cor.test(
    c3_df$Immune_Recruitment_Score,
    c3_df$Immune_Suppression_Score,
    method = "spearman",
    exact = FALSE
  )
  
  p_cor_c3 <- ggplot(
    c3_df,
    aes(Immune_Recruitment_Score, Immune_Suppression_Score)
  ) +
    geom_point(
      size = 2.8, alpha = 0.9, shape = 21,
      fill = SUBTYPE_COLORS["C3"], color = "white", stroke = 0.3
    ) +
    geom_smooth(
      method = "lm", se = TRUE,
      color = SUBTYPE_COLORS["C3"],
      linewidth = 1, alpha = 0.15
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = sprintf(
        "C3 only\nSpearman rho = %.3f\n%s",
        unname(cor_c3$estimate), format_pval2(cor_c3$p.value)
      ),
      hjust = 1.08, vjust = 1.4, size = 3.6, color = "black"
    ) +
    labs(
      x = "Immune Recruitment Score",
      y = "Immune Suppression Score",
      title = "C3: Recruitment vs Suppression"
    ) +
    theme_pub()
  
  save_plot(p_cor_c3, "Fig_C3_Recruitment_vs_Suppression", width = 6, height = 5)
  
  # 3.2 C3 四象限图
  med_rec <- median(c3_df$Immune_Recruitment_Score, na.rm = TRUE)
  med_sup <- median(c3_df$Immune_Suppression_Score, na.rm = TRUE)
  
  c3_df <- c3_df %>%
    dplyr::mutate(
      Quadrant = dplyr::case_when(
        Immune_Recruitment_Score >= med_rec & Immune_Suppression_Score >= med_sup ~ "High recruit / High suppress",
        Immune_Recruitment_Score >= med_rec & Immune_Suppression_Score <  med_sup ~ "High recruit / Low suppress",
        Immune_Recruitment_Score <  med_rec & Immune_Suppression_Score >= med_sup ~ "Low recruit / High suppress",
        TRUE ~ "Low recruit / Low suppress"
      )
    )
  
  quad_count <- c3_df %>%
    dplyr::count(Quadrant)
  
  write.csv(
    quad_count,
    file.path(OUTPUT_DIR, "C3_Recruitment_Suppression_Quadrant_Counts.csv"),
    row.names = FALSE, quote = FALSE
  )
  
  p_quad <- ggplot(c3_df, aes(Immune_Recruitment_Score, Immune_Suppression_Score)) +
    geom_vline(xintercept = med_rec, linetype = "dashed", color = "gray45") +
    geom_hline(yintercept = med_sup, linetype = "dashed", color = "gray45") +
    geom_point(
      size = 2.8, alpha = 0.9, shape = 21,
      fill = SUBTYPE_COLORS["C3"], color = "white", stroke = 0.3
    ) +
    annotate("text", x = med_rec, y = Inf, label = "median recruitment",
             vjust = 1.4, hjust = -0.05, size = 3.2, color = "gray30") +
    annotate("text", x = Inf, y = med_sup, label = "median suppression",
             vjust = -0.3, hjust = 1.05, size = 3.2, color = "gray30") +
    labs(
      x = "Immune Recruitment Score",
      y = "Immune Suppression Score",
      title = "C3 Quadrant: Recruitment and Suppression"
    ) +
    theme_pub()
  
  save_plot(p_quad, "Fig_C3_Recruitment_Suppression_Quadrants", width = 6, height = 5)
  
  # 3.3 可选：校正 Tumor_Purity 后再看 C3 内部相关性
  if ("Tumor_Purity" %in% colnames(c3_df) &&
      sum(complete.cases(c3_df[, c("Immune_Recruitment_Score",
                                   "Immune_Suppression_Score",
                                   "Tumor_Purity")])) >= 10) {
    
    lm_rec <- lm(Immune_Recruitment_Score ~ Tumor_Purity, data = c3_df)
    lm_sup <- lm(Immune_Suppression_Score ~ Tumor_Purity, data = c3_df)
    
    c3_df$Recruitment_resid <- residuals(lm_rec)
    c3_df$Suppression_resid <- residuals(lm_sup)
    
    cor_c3_resid <- cor.test(
      c3_df$Recruitment_resid,
      c3_df$Suppression_resid,
      method = "spearman",
      exact = FALSE
    )
    
    p_cor_c3_resid <- ggplot(
      c3_df,
      aes(Recruitment_resid, Suppression_resid)
    ) +
      geom_point(
        size = 2.8, alpha = 0.9, shape = 21,
        fill = SUBTYPE_COLORS["C3"], color = "white", stroke = 0.3
      ) +
      geom_smooth(
        method = "lm", se = TRUE,
        color = SUBTYPE_COLORS["C3"],
        linewidth = 1, alpha = 0.15
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = sprintf(
          "C3 only, purity-adjusted\nSpearman rho = %.3f\n%s",
          unname(cor_c3_resid$estimate), format_pval2(cor_c3_resid$p.value)
        ),
        hjust = 1.08, vjust = 1.4, size = 3.5, color = "black"
      ) +
      labs(
        x = "Recruitment residual (adjusted for purity)",
        y = "Suppression residual (adjusted for purity)",
        title = "C3: Recruitment vs Suppression after Purity Adjustment"
      ) +
      theme_pub()
    
    save_plot(
      p_cor_c3_resid,
      "Fig_C3_Recruitment_vs_Suppression_PurityAdjusted",
      width = 6.2, height = 5
    )
  }
}

# 3.4 C3 内部：Recruitment vs Checkpoint_Score
if (all(c("Immune_Recruitment_Score", "Checkpoint_Score") %in% colnames(meta_step2))) {
  
  c3_cp_df <- meta_step2 %>%
    dplyr::filter(Subtype == "C3") %>%
    dplyr::select(SampleID, Immune_Recruitment_Score, Checkpoint_Score) %>%
    dplyr::filter(!is.na(Immune_Recruitment_Score), !is.na(Checkpoint_Score))
  
  if (nrow(c3_cp_df) >= 6) {
    cor_cp <- cor.test(
      c3_cp_df$Immune_Recruitment_Score,
      c3_cp_df$Checkpoint_Score,
      method = "spearman",
      exact = FALSE
    )
    
    p_cp <- ggplot(c3_cp_df, aes(Immune_Recruitment_Score, Checkpoint_Score)) +
      geom_point(
        size = 2.8, alpha = 0.9, shape = 21,
        fill = SUBTYPE_COLORS["C3"], color = "white", stroke = 0.3
      ) +
      geom_smooth(
        method = "lm", se = TRUE,
        color = SUBTYPE_COLORS["C3"],
        linewidth = 1, alpha = 0.15
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = sprintf(
          "C3 only\nSpearman rho = %.3f\n%s",
          unname(cor_cp$estimate), format_pval2(cor_cp$p.value)
        ),
        hjust = 1.08, vjust = 1.4, size = 3.6, color = "black"
      ) +
      labs(
        x = "Immune Recruitment Score",
        y = "Checkpoint Score",
        title = "C3: Recruitment vs Checkpoint Program"
      ) +
      theme_pub()
    
    save_plot(p_cp, "Fig_C3_Recruitment_vs_Checkpoint", width = 5, height = 5)
  }
}

message("SECTION 32 complete.")




getwd()
setwd("/media/desk16/iy5111/VS.sc/vs.paper/VS.4")
# ==========================================================================
# 导出 CIBERSORT 22 种免疫细胞在 C1/C2/C3 三组间的完整统计表
# ==========================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(FSA)      # Dunn's test
})

# ---- 1. 读取 CIBERSORT 结果 + 分组 ----
res_raw <- read.csv("Figures/step8_Immune/CIBERSORT_Result.csv", row.names = 1, check.names = FALSE)

# 抓 P-value 列（兼容 P-value / P.value）
pval_col <- grep("P-value|P.value", colnames(res_raw), value = TRUE, ignore.case = TRUE)[1]

# 两套样本集：全部 vs P<0.05 过滤后
n_total    <- nrow(res_raw)
keep_idx   <- which(res_raw[[pval_col]] < 0.05)
n_retained <- length(keep_idx)
cat(sprintf("CIBERSORT deconvolution: %d/%d samples passed P < 0.05 filter.\n",
            n_retained, n_total))

# 22 种细胞比例
cell_mat <- res_raw[, 1:22]
cell_mat$Sample <- rownames(cell_mat)

# 分组信息
group <- read.csv("Figures/step8_Immune/group.csv")
colnames(group) <- c("Sample", "Group")

# ---- 2. 两种样本集分别分析 ----
run_stats <- function(df, label) {
  long <- df %>%
    pivot_longer(-c(Sample, Group), names_to = "CellType", values_to = "Fraction")
  
  # 每种细胞：中位数 IQR + KW + Dunn
  out <- long %>%
    group_by(CellType) %>%
    group_modify(~ {
      dat <- .x
      # 每组中位数 [IQR]
      stats_per_group <- dat %>%
        group_by(Group) %>%
        summarise(
          n      = n(),
          median = median(Fraction, na.rm = TRUE),
          q25    = quantile(Fraction, 0.25, na.rm = TRUE),
          q75    = quantile(Fraction, 0.75, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(summary = sprintf("%.4f [%.4f–%.4f] (n=%d)", median, q25, q75, n))
      
      wide_summary <- stats_per_group %>%
        select(Group, summary) %>%
        pivot_wider(names_from = Group, values_from = summary,
                    names_prefix = "Median_IQR_")
      
      # Kruskal-Wallis
      kw_p <- tryCatch(
        kruskal.test(Fraction ~ factor(Group), data = dat)$p.value,
        error = function(e) NA_real_
      )
      
      # Dunn pairwise
      dunn_res <- tryCatch({
        d <- dunnTest(Fraction ~ factor(Group, levels = c("C1","C2","C3")),
                      data = dat, method = "bh")$res
        setNames(d$P.adj, d$Comparison)
      }, error = function(e) c())
      
      data.frame(
        wide_summary,
        KW_p      = kw_p,
        Dunn_C1_C2 = ifelse("C1 - C2" %in% names(dunn_res), dunn_res["C1 - C2"], NA),
        Dunn_C1_C3 = ifelse("C1 - C3" %in% names(dunn_res), dunn_res["C1 - C3"], NA),
        Dunn_C2_C3 = ifelse("C2 - C3" %in% names(dunn_res), dunn_res["C2 - C3"], NA),
        check.names = FALSE
      )
    }) %>%
    ungroup() %>%
    mutate(
      KW_signif = case_when(
        is.na(KW_p)      ~ "NA",
        KW_p < 0.001     ~ "***",
        KW_p < 0.01      ~ "**",
        KW_p < 0.05      ~ "*",
        TRUE             ~ "ns"
      ),
      Sample_set = label
    ) %>%
    arrange(KW_p)
  
  return(out)
}

# 全部 38 样本
df_all <- merge(cell_mat, group, by = "Sample")
res_all <- run_stats(df_all, sprintf("All samples (n=%d)", n_total))

# P<0.05 过滤后样本
if (n_retained > 0) {
  cell_mat_f <- cell_mat[rownames(res_raw)[keep_idx], ]
  df_filt    <- merge(cell_mat_f, group, by = "Sample")
  res_filt   <- run_stats(df_filt, sprintf("Filtered P<0.05 (n=%d)", n_retained))
} else {
  res_filt <- NULL
}

# ---- 3. 合并并导出 ----
res_final <- bind_rows(res_all, res_filt)

write.csv(res_final,
          "CIBERSORT_Subtype_Comparison_Table.csv",
          row.names = FALSE)

cat("\n✅ 结果已导出: CIBERSORT_Subtype_Comparison_Table.csv\n")
cat("\n=== 全部 38 样本，KW P<0.05 的细胞 ===\n")
print(res_all %>% filter(KW_p < 0.05) %>%
        select(CellType, starts_with("Median_IQR"), KW_p, KW_signif,
               Dunn_C1_C2, Dunn_C1_C3, Dunn_C2_C3))

if (!is.null(res_filt)) {
  cat("\n=== 过滤后样本，KW P<0.05 的细胞 ===\n")
  print(res_filt %>% filter(KW_p < 0.05) %>%
          select(CellType, starts_with("Median_IQR"), KW_p, KW_signif,
                 Dunn_C1_C2, Dunn_C1_C3, Dunn_C2_C3))
}

\









# ==============================================================================
# CIBERSORT 优质数据 (P < 0.05) 全景可视化图集
# 包含：数据清洗、堆积条形图、组间差异箱线图、细胞相关性热图、丰度全景热图
# ==============================================================================
getwd()
# === 0. 加载必要的极客绘图包 ===
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

# === 1. 数据读取与核心过滤 (只保留 P < 0.05) ===
message("--- 1. 数据加载与过滤 ---")
res_raw <- read.csv("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/step8_Immune/CIBERSORT_Result.csv", row.names = 1, check.names = FALSE)
pval_col <- grep("P-value|P.value", colnames(res_raw), value = TRUE, ignore.case = TRUE)[1]

# 核心：过滤并提取前 22 列免疫细胞比例
res_clean <- res_raw[res_raw[[pval_col]] < 0.05, 1:22]
res_clean$Sample <- rownames(res_clean)

message(sprintf("过滤完成：共保留 %d 个优质样本 (P < 0.05)。", nrow(res_clean)))

# === 2. 读取分组信息并合并 ===
# 假设你的分组文件里有 Sample 和 Subtype 两列 (对应 C1, C2, C3)
group_info <- read.csv("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/step8_Immune/group.csv") 
# 如果你的列名叫 Group，请把它改成 Subtype 保持统一
if("Group" %in% colnames(group_info)) colnames(group_info)[colnames(group_info)=="Group"] <- "Subtype"

# 只保留在优质数据里的分组信息，并合并
res_merged <- merge(res_clean, group_info, by = "Sample")

# 数据宽转长 (用于 ggplot)
res_melt <- melt(res_merged, 
                 id.vars = c("Sample", "Subtype"), 
                 variable.name = "CellType", 
                 value.name = "Proportion")

# === 3. 设定高颜值统一配色 ===
# 亚型颜色
SUBTYPE_COLORS <- c("C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35")

# 22 种细胞专属马卡龙配色
my_22_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#E31A1C", "#BC80BD", "#CCEBC5", "#FFED6F",
  "#1F78B4", "#A9DAD1", "#D9D9D9", "#FB9A99", "#FCCDE5", "#6A3D9A", "#B15928",
  "#A6CEE3", "#B2DF8A", "#FF7F00", "#FDBF6F"
)


# ==============================================================================
# 📊 图 A: 免疫细胞组分堆积条形图 (Landscape)
# ==============================================================================
message("--- 绘制图 A：堆积条形图 ---")

p_bar <- ggplot(res_melt, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1, color = NA) + 
  facet_grid(~ Subtype, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = my_22_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(), 
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black")
  ) +
  labs(title = "Immune Cell Composition (P < 0.05 Validated Samples)", 
       x = "Samples grouped by Subtype", y = "Relative Fraction")

ggsave("Fig_CIBERSORT_Filtered_Barplot.pdf", plot = p_bar, width = 12, height = 5)


# ==============================================================================
# 📊 图 B: 三组别免疫浸润差异箱线图 (带显著性星号)
# ==============================================================================
message("--- 绘制图 B：差异箱线图 ---")

p_box <- ggplot(res_melt, aes(x = CellType, y = Proportion, fill = Subtype)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 0.5, alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
    legend.position = "top",                                            
    panel.grid.major.x = element_blank()                                
  ) +
  labs(title = "Differential Immune Infiltration Across Subtypes",
       x = NULL, y = "Estimated Proportion", fill = "Subtype") +
  # 自动执行 Kruskal-Wallis 或 Wilcoxon 并在上方打星号
  stat_compare_means(aes(group = Subtype), 
                     label = "p.signif",    
                     method = "kruskal.test", 
                     hide.ns = TRUE,
                     size = 4) 

ggsave("Fig_CIBERSORT_Filtered_Boxplot.pdf", plot = p_box, width = 12, height = 5)


# ==============================================================================
# 📊 图 C: 免疫细胞共表达/排斥相关性热图 (Correlation)
# ==============================================================================
message("--- 绘制图 C：相关性热图 ---")

# 剔除由于样本少导致的方差为 0 的列（否则算相关性会报错）
res_matrix_only <- res_clean[, 1:22]
valid_cells <- apply(res_matrix_only, 2, sd) > 0
res_var_filtered <- res_matrix_only[, valid_cells]

cor_matrix <- cor(res_var_filtered, method = "spearman")
cor_matrix[is.na(cor_matrix)] <- 0 # 兜底防止残存 NA

pdf("Fig_CIBERSORT_Filtered_Correlation.pdf", width = 8, height = 8)
pheatmap(cor_matrix,
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         display_numbers = FALSE, 
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Immune Cell Spearman Correlation")
dev.off()


# ==============================================================================
# 📊 图 D: 丰度聚类全景热图 (Abundance Heatmap)
# ==============================================================================
message("--- 绘制图 D：丰度全景热图 ---")

# 准备热图注释条
annotation_col <- data.frame(Subtype = res_merged$Subtype)
rownames(annotation_col) <- res_merged$Sample

# 设置注释条颜色
ann_colors <- list(Subtype = SUBTYPE_COLORS)

# 为了让热图显示颜色差异更明显，我们画行(基因/细胞)的 Z-score
res_mat_z <- t(scale(res_var_filtered)) 
# 截断极端值防止颜色失真
res_mat_z[res_mat_z > 3] <- 3
res_mat_z[res_mat_z < -3] <- -3

pdf("Fig_CIBERSORT_Filtered_Abundance_Heatmap.pdf", width = 10, height = 6)
pheatmap(res_mat_z, 
         cluster_rows = TRUE,    # 细胞类型聚类
         cluster_cols = TRUE,    # 样本聚类
         color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100), 
         show_colnames = FALSE,  # 隐藏底部样本名让画面干净
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         main = "Immune Cell Abundance (Z-score)",
         fontsize_row = 10)
dev.off()

message("✅ 全部可视化已成功完成，请在当前工作目录检查生成的 4 份 PDF 文件！")





# ==============================================================================
# CIBERSORT 过滤质控分析 (Sample Attrition Analysis)
# 目标：统计并可视化 C1/C2/C3 各组中被过滤掉 (P >= 0.05) 的样本比例
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

message("--- 开始样本流失质控分析 ---")

# === 1. 读取原始结果与分组信息 ===
res_raw <- read.csv("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/step8_Immune/CIBERSORT_Result.csv", row.names = 1, check.names = FALSE)
group_info <- read.csv("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/step8_Immune/group.csv")

# 统一列名防呆
colnames(group_info) <- c("Sample", "Subtype")

# === 2. 提取 P 值并打标签 (Retained vs Filtered) ===
pval_col <- grep("P-value|P.value", colnames(res_raw), value = TRUE, ignore.case = TRUE)[1]

status_df <- data.frame(
  Sample = rownames(res_raw),
  P_value = res_raw[[pval_col]]
) %>%
  mutate(Status = ifelse(P_value < 0.05, "Retained (P < 0.05)", "Filtered Out (P >= 0.05)"))

# 合并分组信息
merged_status <- merge(status_df, group_info, by = "Sample")

# === 3. 生成统计汇总表 ===
summary_table <- merged_status %>%
  group_by(Subtype, Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(
    Total_in_Subtype = sum(Count),
    Percentage = round(Count / Total_in_Subtype * 100, 1),
    Label = sprintf("%d (%.1f%%)", Count, Percentage)
  )

cat("\n=== 各亚型样本过滤情况汇总 ===\n")
print(summary_table %>% select(Subtype, Status, Count, Percentage))

# === 4. 统计学检验 (Fisher's Exact Test) ===
# 检验“被过滤的概率”是否与“肿瘤亚型”显著相关
table_for_test <- table(merged_status$Subtype, merged_status$Status)
fisher_res <- fisher.test(table_for_test)

cat(sprintf("\n[统计检验] 亚型与过滤率的相关性 Fisher's exact test P-value = %.4f\n", fisher_res$p.value))

if(fisher_res$p.value < 0.05) {
  cat("💡 结论：不同亚型被过滤的比例存在【显著差异】。这可能暗示某些亚型(如高比例过滤组)属于免疫冷肿瘤！\n")
} else {
  cat("💡 结论：不同亚型被过滤的比例没有显著差异，样本流失是相对均匀的。\n")
}

# === 5. 绘制可视化堆叠条形图 ===
# 定义颜色：保留的用沉稳的灰色，过滤掉的用醒目的红色
status_colors <- c("Retained (P < 0.05)" = "#A0A0A0", "Filtered Out (P >= 0.05)" = "#E64B35")

p_attrition <- ggplot(summary_table, aes(x = Subtype, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.5) +
  # 在柱子上添加具体的数量和百分比标签
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), 
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = status_colors) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "CIBERSORT Sample Attrition by Subtype",
    subtitle = sprintf("Fisher's Exact Test, P = %.4f", fisher_res$p.value),
    x = "Predicted Subtype",
    y = "Proportion of Samples (%)"
  )

# 保存图片
ggsave("Fig_CIBERSORT_Sample_Attrition.pdf", plot = p_attrition, width = 6, height = 6)
message("\n✅ 质控分析完成，已保存统计图：Fig_CIBERSORT_Sample_Attrition.pdf")
