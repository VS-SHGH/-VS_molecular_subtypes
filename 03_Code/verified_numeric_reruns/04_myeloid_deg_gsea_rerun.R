# Reproduce Seurat 5.0.1's default Wilcoxon/limma rank-sum branch for the
# C3-high versus C3-low myeloid comparison. The Rcpp helper implements the
# same continuity-corrected rank-sum calculation used by
# Seurat:::WilcoxDETest when limma is installed. No figures are generated.

options(stringsAsFactors = FALSE)
source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
ROOT <- package_root
OUT <- file.path(AUDIT_ROOT, "singlecell_cellchat_rerun")
SC_QS <- SC_QS_INPUT
MYELOID_META <- file.path(OUT, "myeloid_scored_metadata.rds")
OLD_DEG <- file.path(LEGACY_FIGURE_DIR, "Step07_C3Myeloid", "DEG_C3high_vs_C3low.csv")

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
  library(dplyr)
  library(Rcpp)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

message("Loading expression object and fresh myeloid labels...")
sc <- qread(SC_QS)
meta <- readRDS(MYELOID_META)
cells <- intersect(colnames(sc), rownames(meta))
sc <- subset(sc, cells = cells)
sc@meta.data <- meta[cells, , drop = FALSE]
Idents(sc) <- "C3_status"
cells1 <- WhichCells(sc, idents = "C3-high")
cells2 <- WhichCells(sc, idents = "C3-low")
n1 <- length(cells1)
n2 <- length(cells2)
message("C3-high=", n1, "; C3-low=", n2)

data_mat <- GetAssayData(sc, assay = "RNA", slot = "data")
# The saved manuscript-era DEG table uses the older Seurat fold-change
# convention log2(mean(expm1(data)) + pseudocount), as confirmed by direct
# comparison to the cached values.  Supply that mean function explicitly so
# the current Seurat 5.0.1 runtime does not silently use the newer
# pseudocount-before-division convention.
legacy_log1p_mean <- function(x) {
  log(x = rowSums(expm1(x = x)) / NCOL(x) + 1, base = 2)
}
fc <- Seurat:::FoldChange.Assay(
  sc[["RNA"]], cells.1 = cells1, cells.2 = cells2,
  features = rownames(data_mat), slot = "data",
  pseudocount.use = 1, base = 2, mean.fxn = legacy_log1p_mean
)
# FoldChange.Assay can return a data.frame whose row names are partially
# corrupted on this Seurat/R combination.  The requested feature order is
# authoritative and is preserved explicitly before filtering.
fc <- as.data.frame(fc)
stopifnot(nrow(fc) == nrow(data_mat))
fc$gene <- rownames(data_mat)
rownames(fc) <- rownames(data_mat)
keep <- pmax(fc$pct.1, fc$pct.2) >= 0.1 &
  abs(fc$avg_log2FC) >= 0.25
fc_test <- fc[keep, , drop = FALSE]
tested <- fc_test$gene
tested <- tested[!is.na(tested) & tested != "" & tested != "NA"]
fc_test <- fc_test[match(tested, fc_test$gene), , drop = FALSE]
message("Genes passing Seurat min.pct and logfc.threshold: ", length(tested))

# R's rankSumTestWithCorrelation(index, statistics, correlation=0) formula.
# The input matrix is columns=cells1 followed by cells2, so group 1 is the
# first n1 columns. Zero values are handled as one tied rank group.
if (identical(unname(Sys.info()["sysname"]), "Darwin") && nzchar(Sys.which("xcrun"))) {
  sdk <- system2("xcrun", "--show-sdk-path", stdout = TRUE)
  Sys.setenv(PKG_CXXFLAGS = paste(
    "-arch x86_64 -isysroot", sdk,
    "-isystem", file.path(sdk, "usr/include/c++/v1")
  ))
}
cppFunction(
  code = '
  NumericVector fast_rank_sum_p(NumericMatrix X, int n1) {
    int nr = X.nrow();
    int n = X.ncol();
    int n2 = n - n1;
    NumericVector out(nr);
    for (int i = 0; i < nr; ++i) {
      std::vector< std::pair<double, int> > positive;
      positive.reserve(n);
      for (int j = 0; j < n; ++j) {
        double v = X(i, j);
        if (v != 0.0) positive.push_back(std::make_pair(v, j));
      }
      std::sort(positive.begin(), positive.end(),
                [](const std::pair<double, int>& a,
                   const std::pair<double, int>& b) {
                  return a.first < b.first;
                });
      int nzero = n - (int)positive.size();
      int nzero_group1 = n1;
      for (size_t k = 0; k < positive.size(); ++k)
        if (positive[k].second < n1) --nzero_group1;
      double sum_r1 = nzero_group1 * (nzero + 1) / 2.0;
      // The zero-expression values form one large tied rank group and must
      // contribute to the limma tie correction as well.
      double tie_term = (double)nzero * nzero * nzero - nzero;
      size_t k = 0;
      while (k < positive.size()) {
        size_t e = k + 1;
        while (e < positive.size() && positive[e].first == positive[k].first)
          ++e;
        double t = (double)(e - k);
        double avg_rank = nzero + ((double)(k + 1) + (double)e) / 2.0;
        for (size_t q = k; q < e; ++q)
          if (positive[q].second < n1) sum_r1 += avg_rank;
        tie_term += t * t * t - t;
        k = e;
      }
      double U = (double)n1 * n2 + (double)n1 * (n1 + 1) / 2.0 - sum_r1;
      double mu = (double)n1 * n2 / 2.0;
      double sigma2 = (double)n1 * n2 * (n + 1) / 12.0;
      if (n > 1)
        sigma2 *= 1.0 - tie_term / ((double)n * (n + 1) * (n - 1));
      if (!R_FINITE(sigma2) || sigma2 <= 0.0) {
        out[i] = 1.0;
      } else {
        double zlow = (U + 0.5 - mu) / std::sqrt(sigma2);
        double zup = (U - 0.5 - mu) / std::sqrt(sigma2);
        double p_less = R::pnorm5(zup, 0.0, 1.0, false, false);
        double p_greater = R::pnorm5(zlow, 0.0, 1.0, true, false);
        out[i] = std::min(1.0, 2.0 * std::min(p_less, p_greater));
      }
    }
    return out;
  }
'
)

cell_idx <- match(c(cells1, cells2), colnames(data_mat))
gene_idx <- match(tested, rownames(data_mat))
stopifnot(!anyNA(cell_idx), !anyNA(gene_idx))
X <- as.matrix(data_mat[gene_idx, cell_idx, drop = FALSE])
rownames(X) <- tested
message("Computing tied-rank p-values in C++...")
pvals <- fast_rank_sum_p(X, n1)
deg <- fc_test
deg$p_val <- pvals
# FindMarkers.default uses n=nrow(object), not the number of tested genes.
deg$p_val_adj <- p.adjust(pvals, method = "bonferroni", n = nrow(data_mat))
deg <- deg[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene")]
rownames(deg) <- deg$gene
write.csv(deg, file.path(OUT, "DEG_C3high_vs_C3low_rcpp_all.csv"), row.names = FALSE)
deg_sig <- deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)
write.csv(deg_sig, file.path(OUT, "DEG_C3high_vs_C3low_rcpp_sig.csv"), row.names = FALSE)

comparison <- data.frame(
  metric = c("genes_passing_Seurat_filters", "significant_genes_rcpp"),
  rerun = c(nrow(deg), nrow(deg_sig)),
  cached = c(NA, if (file.exists(OLD_DEG)) nrow(read.csv(OLD_DEG, check.names = FALSE)) else NA)
)
if (file.exists(OLD_DEG)) {
  old <- read.csv(OLD_DEG, check.names = FALSE)
  if ("gene" %in% colnames(old)) {
    common_genes <- intersect(deg_sig$gene, old$gene)
    comparison <- rbind(
      comparison,
      data.frame(
        metric = c("sig_gene_overlap", "sig_gene_jaccard"),
        rerun = c(length(common_genes),
                  length(common_genes) / length(union(deg_sig$gene, old$gene))),
        cached = c(nrow(old), NA)
      )
    )
    if (length(common_genes)) {
      a <- deg_sig[match(common_genes, deg_sig$gene), ]
      b <- old[match(common_genes, old$gene), ]
      comparison <- rbind(
        comparison,
        data.frame(
          metric = c("sig_avg_log2FC_max_abs_diff", "sig_pct1_max_abs_diff",
                     "sig_pct2_max_abs_diff", "sig_p_adj_max_abs_diff"),
          rerun = c(max(abs(a$avg_log2FC - b$avg_log2FC), na.rm = TRUE),
                    max(abs(a$pct.1 - b$pct.1), na.rm = TRUE),
                    max(abs(a$pct.2 - b$pct.2), na.rm = TRUE),
                    max(abs(a$p_val_adj - b$p_val_adj), na.rm = TRUE)),
          cached = NA
        )
      )
    }
  }
}
write.csv(comparison, file.path(OUT, "DEG_rcpp_cache_comparison.csv"), row.names = FALSE)

message("Running GO BP GSEA from Rcpp DEG output...")
gene_list <- deg$avg_log2FC
names(gene_list) <- deg$gene
gene_list <- sort(gene_list, decreasing = TRUE)
# Match the original script's explicit GSEA cleaning step: genes with
# NA/Inf statistics are excluded before calling clusterProfiler::gseGO.
gene_list <- gene_list[is.finite(gene_list) & !is.na(names(gene_list)) &
                         names(gene_list) != ""]
gene_list <- gene_list[!duplicated(names(gene_list))]
gsea_go <- gseGO(
  geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "BP", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500,
  verbose = FALSE
)
gsea_df <- as.data.frame(gsea_go)
write.csv(gsea_df, file.path(OUT, "GSEA_GO_BP_C3high_vs_C3low_rcpp.csv"), row.names = FALSE)
if (nrow(gsea_df)) {
  write.csv(bind_rows(
    gsea_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = 5) %>%
      mutate(Direction = "Up in C3-high"),
    gsea_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = 10) %>%
      mutate(Direction = "Down in C3-high")
  ), file.path(OUT, "GSEA_GO_BP_top_directional_rcpp.csv"), row.names = FALSE)
}
message("Rcpp myeloid DEG/GSEA rerun finished.")
