source("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/00_config.R")
suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(pheatmap)
  library(dplyr)
  library(ggplot2)
  library(factoextra)
  library(Rtsne)
})

# ==============================================================================
# 1. 加载清洗后的数据
# ==============================================================================
data_file <- file.path(RDATA_DIR, "Step1_Clean_Data_VS3.RData")
load(data_file)


log_tpm <- log2(clean_tpm + 1)

mads <- apply(log_tpm, 1, mad)
n_select <- min(PARAMS$n_top_mad_genes, sum(mads > 0))
top_genes <- names(sort(mads, decreasing = TRUE))[1:n_select]


clust_input <- log_tpm[top_genes, ]

clust_input <- sweep(clust_input, 1, apply(clust_input, 1, median, na.rm = TRUE))



# ==============================================================================
# 3
# ==============================================================================
step2_fig_dir <- file.path(FIG_DIR, "Step2")

message("onsensusClusterPlus (reps=", PARAMS$ccp_reps, ") ...")

ccp_results <- ConsensusClusterPlus(
  as.matrix(clust_input),
  maxK         = PARAMS$ccp_maxK,
  reps         = PARAMS$ccp_reps,
  pItem        = PARAMS$ccp_pItem,
  pFeature     = 1,
  title        = step2_fig_dir,
  clusterAlg   = "hc",
  distance     = "pearson",
  seed         = PARAMS$ccp_seed,
  plot         = "pdf"        
)

message("✅ 共识聚类完成")

# ==============================================================================
# 4. k 值选择的客观评估
# ==============================================================================
calc_pac <- function(ccp_res, k, lower = 0.1, upper = 0.9) {
  consensus_mat <- ccp_res[[k]]$consensusMatrix
  vals <- consensus_mat[lower.tri(consensus_mat)]
  pac <- mean(vals > lower & vals < upper)
  return(pac)
}

pac_values <- sapply(2:PARAMS$ccp_maxK, function(k) calc_pac(ccp_results, k))
names(pac_values) <- paste0("k=", 2:PARAMS$ccp_maxK)

message("📊 PAC 值 (越小越好):")
print(round(pac_values, 4))

# 绘制 PAC 图
pac_df <- data.frame(k = 2:PARAMS$ccp_maxK, PAC = pac_values)

p_pac <- ggplot(pac_df, aes(x = k, y = PAC)) +
  geom_line(linewidth = 0.8, colour = "grey30") +
  geom_point(size = 3, colour = "firebrick") +
  geom_vline(xintercept = PARAMS$k_final, linetype = "dashed", colour = "#E64B35") +
  annotate("text", x = PARAMS$k_final + 0.3, y = max(pac_df$PAC) * 0.9,
           label = paste0("k = ", PARAMS$k_final), colour = "#E64B35",
           fontface = "bold", size = 4.5) +
  scale_x_continuous(breaks = 2:PARAMS$ccp_maxK) +
  labs(x = "Number of Clusters (k)", y = "PAC Score",
       title = "Optimal k Selection (PAC)") +
  theme_publication(base_size = 13)

save_pdf(p_pac, file.path(step2_fig_dir, "PAC_k_Selection.pdf"), width = 5, height = 4)

# ==============================================================================
# 5. 
# ==============================================================================
k <- PARAMS$k_final
cluster_labels <- ccp_results[[k]]$consensusClass

meta_step2 <- meta_qc
meta_step2$Subtype_Num <- cluster_labels[meta_step2$SampleID]
meta_step2$Subtype     <- paste0("C", meta_step2$Subtype_Num)

message("📊 各亚型样本数:")
print(table(meta_step2$Subtype))
# 
# ==============================================================================
# 6
# ==============================================================================
# 按亚型排序
meta_ordered <- meta_step2[order(meta_step2$Subtype), ]
clust_ordered <- clust_input[, meta_ordered$SampleID]

plot_mat <- clust_ordered[top500, ]

# 1. 检查 NA / NaN / Inf
message("NA 数量: ", sum(is.na(plot_mat)))
message("NaN 数量: ", sum(is.nan(plot_mat)))
message("Inf 数量: ", sum(is.infinite(plot_mat)))

# 2. 
row_sds <- apply(plot_mat, 1, sd, na.rm = TRUE)
message("零方差基因数: ", sum(row_sds == 0 | is.na(row_sds)))

# 3. 查看哪些基因有问题
bad_genes <- names(which(row_sds == 0 | is.na(row_sds)))
if (length(bad_genes) > 0) print(head(bad_genes, 20))
ann_col <- data.frame(Subtype = meta_ordered$Subtype)
rownames(ann_col) <- meta_ordered$SampleID

if ("Tumor_Score" %in% colnames(meta_ordered)) {
  ann_col$TumorPurity <- meta_ordered$Tumor_Score
}

ann_colors <- list(
  Subtype     = SUBTYPE_COLORS[sort(unique(meta_ordered$Subtype))],
  TumorPurity = c("white", "firebrick")
)

# Top 500 
top500 <- head(rownames(clust_ordered), 500)
break_list <- seq(-2, 2, length.out = 101)

# 
gap_positions <- cumsum(table(meta_ordered$Subtype))
gap_positions <- gap_positions[-length(gap_positions)]

pdf(file.path(step2_fig_dir, "Heatmap_Sorted_Publication.pdf"),
    width = 9, height = 11)
pheatmap(clust_ordered[top500, ],
         scale          = "row",
         breaks         = break_list,
         color          = HEATMAP_COLORS,
         cluster_cols   = FALSE,
         cluster_rows   = TRUE,
         clustering_method = "ward.D2",
         gaps_col       = gap_positions,
         show_rownames  = FALSE,
         show_colnames  = FALSE,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         main           = "Molecular Subtypes (Consensus Clustering)",
         border_color   = NA)
dev.off()
message("✅ 排序热图已保存")




# ==============================================================================
# 7.  PCA 
# ==============================================================================
pca_res <- prcomp(t(log_tpm[top_genes, ]), scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x[, 1:2])
pca_df$SampleID <- rownames(pca_df)
pca_df  <- left_join(pca_df, meta_step2[, c("SampleID", "Subtype")], by = "SampleID")

var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Subtype)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 0.6, linetype = "dashed") +
  scale_colour_manual(values = SUBTYPE_COLORS) +
  labs(
    x     = paste0("PC1 (", var_exp[1], "%)"),
    y     = paste0("PC2 (", var_exp[2], "%)"),
    title = "PCA of Molecular Subtypes"
  ) +
  theme_publication(base_size = 13)

save_half_vector(p_pca,
                 file.path(step2_fig_dir, "PCA_Subtypes.pdf"),
                 width = 6.5, height = 5.5)


# 4. 
ann_col <- data.frame(Subtype = meta_ordered$Subtype,
                      row.names = meta_ordered$SampleID)
ann_colors <- list(Subtype = SUBTYPE_COLORS[sort(unique(meta_ordered$Subtype))])

# 5
pdf(file.path(step2_fig_dir, "Heatmap_Sorted_Publication.pdf"),
    width = 9, height = 11)
pheatmap(plot_mat_scaled,
         scale             = "none",
         breaks            = break_list,
         color             = HEATMAP_COLORS,
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         clustering_method = "ward.D2",
         gaps_col          = gap_positions,
         annotation_col    = ann_col,
         annotation_colors = ann_colors,
         show_rownames     = FALSE,
         show_colnames     = FALSE,
         main              = "Molecular Subtypes (Consensus Clustering)",
         border_color      = NA)
dev.off()






# ==============================================================================
# Step 3


suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
})


# ==============================================================================
load(file.path(RDATA_DIR, "Step2_Clustering.RData"))

# 确保对齐
clean_tpm <- clean_tpm[, meta_step2$SampleID]
stopifnot(all(colnames(clean_tpm) == meta_step2$SampleID))


# ==============================================================================
# 2. limma-trend 
# ==============================================================================
log_tpm <- log2(clean_tpm + 1)
group <- factor(meta_step2$Subtype, levels = c("C1", "C2", "C3"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- meta_step2$SampleID

fit <- lmFit(log_tpm, design)

contrasts_mat <- makeContrasts(
  C1_vs_Others = C1 - (C2 + C3) / 2,
  C2_vs_Others = C2 - (C1 + C3) / 2,
  C3_vs_Others = C3 - (C1 + C2) / 2,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2, trend = TRUE) 

# ==============================================================================
# 3
# ==============================================================================
dea_dir <- file.path(TABLE_DIR, "Step3_DEA")
if (!dir.exists(dea_dir)) dir.create(dea_dir, recursive = TRUE)

extract_deg <- function(fit_obj, coef_name) {
  res <- topTable(fit_obj, coef = coef_name, number = Inf, adjust.method = "BH")
  res$Gene <- rownames(res)
  res$Change <- case_when(
    res$adj.P.Val < PARAMS$dea_fdr_cutoff & res$logFC >  PARAMS$dea_logfc_cutoff ~ "Up",
    res$adj.P.Val < PARAMS$dea_fdr_cutoff & res$logFC < -PARAMS$dea_logfc_cutoff ~ "Down",
    TRUE ~ "Stable"
  )
  return(res)
}

deg_c1 <- extract_deg(fit2, "C1_vs_Others")
deg_c2 <- extract_deg(fit2, "C2_vs_Others")
deg_c3 <- extract_deg(fit2, "C3_vs_Others")

# 
for (nm in c("C1", "C2", "C3")) {
  deg <- get(paste0("deg_c", gsub("C", "", nm)))
  message(nm, ": Up=", sum(deg$Change == "Up"),
          " Down=", sum(deg$Change == "Down"),
          " Stable=", sum(deg$Change == "Stable"))
}

# 
write.csv(deg_c1, file.path(dea_dir, "DEG_C1_vs_Others.csv"), row.names = FALSE)
write.csv(deg_c2, file.path(dea_dir, "DEG_C2_vs_Others.csv"), row.names = FALSE)
write.csv(deg_c3, file.path(dea_dir, "DEG_C3_vs_Others.csv"), row.names = FALSE)

# ==============================================================================
# 4. 
n <- PARAMS$dea_top_n_markers

top_c1 <- deg_c1 %>% filter(Change == "Up") %>%
  slice_max(order_by = logFC, n = n) %>% pull(Gene)
top_c2 <- deg_c2 %>% filter(Change == "Up") %>%
  slice_max(order_by = logFC, n = n) %>% pull(Gene)
top_c3 <- deg_c3 %>% filter(Change == "Up") %>%
  slice_max(order_by = logFC, n = n) %>% pull(Gene)

signature_genes <- c(top_c1, top_c2, top_c3)
message("🧬 签名基因: C1=", length(top_c1),
        " C2=", length(top_c2), " C3=", length(top_c3))

# ==============================================================================
# 5.
# ==============================================================================
fig3_dir <- file.path(FIG_DIR, "Step3_DEA")

meta_ordered <- meta_step2[order(meta_step2$Subtype), ]
expr_plot    <- log_tpm[signature_genes, meta_ordered$SampleID]

ann_col <- data.frame(Subtype = meta_ordered$Subtype)
rownames(ann_col) <- meta_ordered$SampleID
ann_colors <- list(Subtype = SUBTYPE_COLORS[sort(unique(meta_ordered$Subtype))])

gap_cols <- cumsum(table(meta_ordered$Subtype))
gap_cols <- gap_cols[-length(gap_cols)]
gap_rows <- c(length(top_c1), length(top_c1) + length(top_c2))

break_list <- seq(-2, 2, length.out = 101)

pdf(file.path(fig3_dir, "Signature_Heatmap.pdf"), width = 8, height = 10)
pheatmap(expr_plot,
         scale            = "row",
         breaks           = break_list,
         color            = HEATMAP_COLORS,
         cluster_cols     = FALSE,
         cluster_rows     = FALSE,
         gaps_col         = gap_cols,
         gaps_row         = gap_rows,
         show_colnames    = FALSE,
         show_rownames    = TRUE,
         fontsize_row     = 7,
         annotation_col   = ann_col,
         annotation_colors = ann_colors,
         main             = paste0("Top ", n, " Marker Genes per Subtype"),
         border_color     = NA)
dev.off()

# ==============================================================================
# 6. 可视化
# ==============================================================================
plot_volcano <- function(deg_df, subtype_name, colour) {
  ggplot(deg_df, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(colour = Change), size = 0.8, alpha = 0.6) +
    scale_colour_manual(
      values = c("Up" = colour, "Down" = "steelblue", "Stable" = "grey80")
    ) +
    geom_hline(yintercept = -log10(PARAMS$dea_fdr_cutoff),
               linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-PARAMS$dea_logfc_cutoff, PARAMS$dea_logfc_cutoff),
               linetype = "dashed", colour = "grey40") +
    labs(
      x     = expression(log[2]~Fold~Change),
      y     = expression(-log[10]~FDR),
      title = paste0(subtype_name, " vs Others")
    ) +
    theme_publication(base_size = 12) +
    theme(legend.position = "none")
}

p_vol_c1 <- plot_volcano(deg_c1, "C1", SUBTYPE_COLORS["C1"])
p_vol_c2 <- plot_volcano(deg_c2, "C2", SUBTYPE_COLORS["C2"])
p_vol_c3 <- plot_volcano(deg_c3, "C3", SUBTYPE_COLORS["C3"])

p_vol_all <- cowplot::plot_grid(p_vol_c1, p_vol_c2, p_vol_c3,
                                ncol = 3, labels = c("C1", "C2", "C3"))

save_half_vector(p_vol_all,
                 file.path(fig3_dir, "Volcano_All_Subtypes.pdf"),
                 width = 14, height = 4.5)

# ==============================================================================
# 7. 
# ==============================================================================
save(deg_c1, deg_c2, deg_c3, signature_genes, top_c1, top_c2, top_c3,
     file = file.path(RDATA_DIR, "Step3_DEA_Results.RData"))

message("🎉 Step 3 完成!")







suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(GSVA)
  library(AnnotationDbi)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(pheatmap)
  library(ggpubr)
})

load(file.path(RDATA_DIR, "Step3_DEA_Results.RData"))
load(file.path(RDATA_DIR, "Step2_Clustering.RData"))

message("✅ 数据加载成功")


# --- Top N 基因提取 (ORA 用) ---
prepare_topn_genes <- function(deg_table, group_name,
                               top_n = PARAMS$enrich_top_n) {
  sig <- deg_table %>%
    filter(adj.P.Val < PARAMS$dea_fdr_cutoff & logFC > 0) %>%
    arrange(desc(logFC)) %>%
    head(top_n)
  
  ids <- bitr(sig$Gene, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  message(group_name, ": ", nrow(ids), " 个基因成功转换 ID")
  return(ids$ENTREZID)
}

# --- GSEA 
build_gsea_vector <- function(deg_table) {
  df <- deg_table %>%
    arrange(desc(logFC)) %>%
    dplyr::select(Gene, logFC)
  ids <- bitr(df$Gene, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  df_m <- merge(df, ids, by.x = "Gene", by.y = "SYMBOL")
  df_m <- df_m[!duplicated(df_m$ENTREZID), ]
  vec  <- setNames(df_m$logFC, df_m$ENTREZID)
  sort(vec, decreasing = TRUE)
}

# --- GO  ---
run_go_safe <- function(gene_ids) {
  if (length(gene_ids) < 10) return(NULL)
  enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "BP",
           pAdjustMethod = "BH", pvalueCutoff = 0.05,
           qvalueCutoff = 0.2, readable = TRUE)
}

# ---  KEGG--
run_kegg_safe <- function(gene_ids) {
  if (length(gene_ids) < 10) return(NULL)
  tryCatch(
    enrichKEGG(gene = gene_ids, organism = "hsa", pvalueCutoff = 0.05),
    error = function(e) {
      warning("KEGG 查询失败 (网络问题?): ", e$message)
      return(NULL)
    }
  )
}

# ---  GSEA ---
run_gsea_safe <- function(gsea_vec, label) {
  message("🚀 GSEA: ", label)
  gseGO(geneList = gsea_vec, OrgDb = org.Hs.eg.db, ont = "BP",
        minGSSize = 10, maxGSSize = 500,
        pvalueCutoff = PARAMS$gsea_pval_cutoff,
        verbose = FALSE, seed = PARAMS$global_seed)
}

# ==============================================================================
# 3. ORA: GO + KEGG (Top N 策略)
# ==============================================================================
genes_c1 <- prepare_topn_genes(deg_c1, "C1")
genes_c2 <- prepare_topn_genes(deg_c2, "C2")
genes_c3 <- prepare_topn_genes(deg_c3, "C3")

go_c1 <- run_go_safe(genes_c1)
go_c2 <- run_go_safe(genes_c2)
go_c3 <- run_go_safe(genes_c3)

kegg_c1 <- run_kegg_safe(genes_c1)
kegg_c2 <- run_kegg_safe(genes_c2)
kegg_c3 <- run_kegg_safe(genes_c3)

# 绘制并保存 dotplot
enrich_fig_dir <- file.path(FIG_DIR, "Step4_Enrichment")

plot_save_dotplot <- function(enrich_obj, filename, title_txt, colour) {
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) return(invisible(NULL))
  p <- dotplot(enrich_obj, showCategory = 10, title = title_txt) +
    scale_colour_gradient(low = colour, high = "grey20") +
    theme_publication(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
  save_pdf(p, file.path(enrich_fig_dir, filename), width = 8, height = 6)
}

plot_save_dotplot(go_c1, "GO_C1.pdf", "C1: GO Biological Process", SUBTYPE_COLORS["C1"])
plot_save_dotplot(go_c2, "GO_C2.pdf", "C2: GO Biological Process", SUBTYPE_COLORS["C2"])
plot_save_dotplot(go_c3, "GO_C3.pdf", "C3: GO Biological Process", SUBTYPE_COLORS["C3"])

plot_save_dotplot(kegg_c1, "KEGG_C1.pdf", "C1: KEGG Pathway", SUBTYPE_COLORS["C1"])
plot_save_dotplot(kegg_c2, "KEGG_C2.pdf", "C2: KEGG Pathway", SUBTYPE_COLORS["C2"])
plot_save_dotplot(kegg_c3, "KEGG_C3.pdf", "C3: KEGG Pathway", SUBTYPE_COLORS["C3"])

# ==============================================================================
# 4. GSEA 分析 (全部亚型)
# ==============================================================================
vec_c1 <- build_gsea_vector(deg_c1)
vec_c2 <- build_gsea_vector(deg_c2)
vec_c3 <- build_gsea_vector(deg_c3)

gse_c1 <- run_gsea_safe(vec_c1, "C1")
gse_c2 <- run_gsea_safe(vec_c2, "C2")
gse_c3 <- run_gsea_safe(vec_c3, "C3")

# --- GSEA ---
plot_gsea_mountain <- function(gse_obj, subtype_key, top_k = 3) {
  if (is.null(gse_obj) || nrow(gse_obj@result) == 0) return(invisible(NULL))
  
  pos_df <- gse_obj@result %>% filter(NES > 0) %>% arrange(p.adjust, desc(NES))
  if (nrow(pos_df) == 0) {
    message("⚠️ ", subtype_key, " 无 NES>0 的正向通路")
    return(invisible(NULL))
  }
  
  ids <- head(pos_df$ID, top_k)
  base_col <- SUBTYPE_COLORS[subtype_key]
  cols <- scales::alpha(base_col, c(1, 0.7, 0.45))[seq_along(ids)]
  
  p <- gseaplot2(gse_obj, geneSetID = ids,
                 title = paste0(subtype_key, ": Top Activated Pathways"),
                 color = cols, pvalue_table = TRUE,
                 ES_geom = "line", subplots = 1:2)
  
  # 、
  p[[1]] <- p[[1]] + theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 14),
    axis.text    = element_text(size = 11),
    legend.text  = element_text(size = 11)
  )
  p[[2]] <- p[[2]] + theme(
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 11)
  )
  
  outfile <- file.path(enrich_fig_dir, paste0("GSEA_Mountain_", subtype_key, ".pdf"))
  ggsave(outfile, p, width = 10, height = 6, device = "pdf")
  message("✅ ", outfile)
}

plot_gsea_mountain(gse_c1, "C1")
plot_gsea_mountain(gse_c2, "C2")
plot_gsea_mountain(gse_c3, "C3")

# --- GSEA  ---
merge_gsea_results <- function(gse_list, names_list) {
  dfs <- lapply(seq_along(gse_list), function(i) {
    obj <- gse_list[[i]]
    if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
    df <- obj@result
    df$Subtype <- names_list[i]
    df
  })
  bind_rows(dfs)
}

all_gsea <- merge_gsea_results(
  list(gse_c1, gse_c2, gse_c3),
  c("C1: Proliferative", "C2: Mesenchymal", "C3: Immune")
)

if (nrow(all_gsea) > 0) {
  top_paths <- all_gsea %>%
    filter(NES > 0) %>%
    group_by(Subtype) %>%
    slice_max(order_by = NES, n = 5) %>%
    pull(Description) %>%
    unique()
  
  plot_gsea_df <- all_gsea %>%
    filter(Description %in% top_paths) %>%
    mutate(
      Subtype = factor(Subtype, levels = c("C1: Proliferative",
                                           "C2: Mesenchymal",
                                           "C3: Immune")),
      Description = str_wrap(Description, width = 45)
    )
  
  # 
  y_order <- plot_gsea_df %>%
    arrange(Subtype, NES) %>%
    pull(Description) %>%
    unique()
  plot_gsea_df$Description <- factor(plot_gsea_df$Description, levels = y_order)
  
  p_summary <- ggplot(plot_gsea_df, aes(x = Subtype, y = Description)) +
    geom_point(aes(size = NES, colour = p.adjust)) +
    scale_colour_gradientn(
      colours = c("#E64B35", "#F39B7F", "grey80"),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size(range = c(3, 8)) +
    labs(x = NULL, y = NULL, colour = "FDR", size = "NES",
         title = "Functional Landscape of VS Subtypes") +
    theme_publication(base_size = 12) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  
  save_pdf(p_summary,
           file.path(enrich_fig_dir, "GSEA_Summary_Dotplot.pdf"),
           width = 9, height = 8)
}



# ============================================================
# GSEA Ridgeplot
# ============================================================
library(ggridges)
library(patchwork)
library(stringr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggnewscale)

TOP_N  <- 3    
WRAP_W <- 30   

# ── 1. 选 Top N ──────────────────────
pick_top <- function(gse, grp, n = TOP_N) {
  gse@result %>%
    filter(!is.na(p.adjust), !is.na(NES), NES > 0) %>%
    arrange(p.adjust, desc(NES)) %>%
    slice_head(n = n) %>%
    mutate(
      Group     = grp,
      neglog10p = -log10(p.adjust + 1e-10),
      Desc_wrap = str_wrap(Description, WRAP_W)
    ) %>%
    dplyr::select(ID, Description, Desc_wrap, p.adjust, neglog10p, NES, Group)
}

top_all <- bind_rows(
  pick_top(gse_c1, "C1"),
  pick_top(gse_c2, "C2"),
  pick_top(gse_c3, "C3")
)
cat("各组通路数:\n"); print(table(top_all$Group))

# ── 2.  ────────────────────────────────────
get_density <- function(gse, grp) {
  sel <- top_all %>% filter(Group == grp)
  gl  <- gse@geneList
  lapply(seq_len(nrow(sel)), function(i) {
    members <- gse@geneSets[[sel$ID[i]]]
    if (is.null(members)) return(NULL)
    vals <- gl[names(gl) %in% members]
    if (!length(vals)) return(NULL)
    data.frame(Description = sel$Description[i], Group = grp,
               x = as.numeric(vals))
  }) %>% bind_rows()
}

density_df <- bind_rows(
  get_density(gse_c1, "C1"),
  get_density(gse_c2, "C2"),
  get_density(gse_c3, "C3")
)

# ── 3. Y ─────────────────────────
path_order <- top_all %>%
  arrange(factor(Group, levels = c("C3","C2","C1")), p.adjust) %>%
  pull(Desc_wrap) %>% unique()

# 
density_df <- density_df %>%
  left_join(
    top_all %>%
      dplyr::select(Description, Desc_wrap, neglog10p, NES, Group) %>%
      distinct(),
    by = c("Description", "Group")
  ) %>%
  mutate(
    Desc_wrap = factor(Desc_wrap, levels = path_order),
    Group     = factor(Group,     levels = c("C1","C2","C3"))
  ) %>%
  filter(!is.na(Desc_wrap))

cat("密度数据各组通路数:\n")
print(table(density_df$Group, density_df$Desc_wrap))

# ── 4. 
x_lo  <- quantile(density_df$x, 0.005, na.rm = TRUE)
x_hi  <- quantile(density_df$x, 0.998, na.rm = TRUE)
NES_X <- x_lo - (x_hi - x_lo) * 0.18   # NES点固定在密度分布左侧

padj_max <- max(density_df$neglog10p, na.rm = TRUE)
nes_lim  <- max(abs(top_all$NES), na.rm = TRUE)

meta_df <- density_df %>%
  distinct(Desc_wrap, Group, NES, neglog10p) %>%
  mutate(x_nes = NES_X)

# ── 5. 分组颜色 ──────────────────────────────────────────────
grp_col <- c(
  C1 = unname(SUBTYPE_COLORS["C1"]),
  C2 = unname(SUBTYPE_COLORS["C2"]),
  C3 = unname(SUBTYPE_COLORS["C3"])
)

# ── 6. ─────────────────────
p_main <- ggplot(density_df, aes(x = x, y = Desc_wrap)) +
  
  # --- Ridge 密度：
  geom_density_ridges(
    aes(fill = neglog10p),
    scale = 1.5, alpha = 0.90,
    rel_min_height = 0.01,
    color = "white", linewidth = 0.25
  ) +
  scale_fill_gradient2(
    low      = "#998ec3",      
    mid      = "#f7f7f7",       
    high     = "#DE712F",       #
    midpoint = padj_max * 0.4,  
    limits   = c(0, padj_max),
    oob      = squish,
    name     = expression(-log[10](p.adjust))
  ) +
  
  # --- 切换 fill scale ---
  new_scale_fill() +
  
  # --- NES 点：
  geom_point(
    data  = meta_df,
    aes(x = x_nes, y = Desc_wrap, fill = NES, size = neglog10p),
    shape = 21, color = "grey20", stroke = 0.4
  ) +
  scale_fill_gradient2(
    low      = "#4575b4",   # 蓝（NES<0）
    mid      = "#f7f7f7",   # 白（NES≈0）
    high     = "#d73027",   # 红（NES>0）
    midpoint = 0,
    limits   = c(-nes_lim, nes_lim),
    name     = "NES"
  ) +
  scale_size(range = c(3, 9), guide = "none") +
  
  # --- 坐标轴（截短 x 轴）---
  coord_cartesian(xlim = c(NES_X - 0.3, x_hi)) +
  
  labs(
    x     = "Ranked Gene Score (logFC)",
    y     = NULL,
    title = "GSEA Ridgeplot"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.border     = element_rect(color = "#000000", linewidth = 1),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text.y      = element_text(size = 8.5, color = "#000000", hjust = 1),
    axis.text.x      = element_text(size = 9,   color = "#000000"),
    axis.title.x     = element_text(size = 11),
    legend.position  = "right",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(4, 6, 4, 0)
  )

# ── 7. ───────────────────────────────
strip_df <- meta_df %>%
  distinct(Desc_wrap, Group) %>%
  mutate(x = 1)

p_strip <- ggplot(strip_df, aes(x = x, y = Desc_wrap, fill = Group)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = grp_col, name = "Subtype") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin     = margin(4, 2, 4, 4)
  )

# ── 8. 单独提取 Subtype 图例 ─────────────────────────────────
p_leg_src <- ggplot(strip_df, aes(x = x, y = Desc_wrap, fill = Group)) +
  geom_tile() +
  scale_fill_manual(values = grp_col, name = "Subtype") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.5, "cm")
  )
grp_legend <- cowplot::get_legend(p_leg_src)

# ── 9. 拼图 ────────────────────────────────────────────────
p_combined <- p_strip + p_main +
  plot_layout(widths = c(1, 20), guides = "collect")

p_final <- cowplot::plot_grid(
  p_combined,
  grp_legend,
  nrow       = 1,
  rel_widths = c(1, 0.1)
)

# ── 10. 保存 ────────────────────────────────────────────────
n_paths <- nrow(top_all)
ggsave(
  file.path(enrich_fig_dir, "GSEA_Ridgeplot_final.pdf"),
  p_final,
  width  = 10,
  height = 1.8 + n_paths * 0.72,
  device = "pdf"
)
message("✅ GSEA_Ridgeplot_final.pdf 已保存")



# ==============================================================================
# 0. 
# ==============================================================================
library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(grid)

select <- dplyr::select

# ==============================================================================
# 1. 准备绘图数据 (保持不变)
# ==============================================================================
go_res <- bind_rows(
  if(exists("go_c1") && !is.null(go_c1) && nrow(go_c1@result) > 0) go_c1@result %>% mutate(Cluster = "C1"),
  if(exists("go_c2") && !is.null(go_c2) && nrow(go_c2@result) > 0) go_c2@result %>% mutate(Cluster = "C2"),
  if(exists("go_c3") && !is.null(go_c3) && nrow(go_c3@result) > 0) go_c3@result %>% mutate(Cluster = "C3")
) %>% select(Cluster, Description, Count, p.adjust)

kegg_res <- bind_rows(
  if(exists("kegg_c1_relaxed") && !is.null(kegg_c1_relaxed) && nrow(kegg_c1_relaxed@result) > 0) {
    kegg_c1_relaxed@result %>% mutate(Cluster = "C1")
  } else if (exists("kegg_c1") && !is.null(kegg_c1) && nrow(kegg_c1) > 0) {
    kegg_c1@result %>% mutate(Cluster = "C1")
  },
  if(exists("kegg_c2") && !is.null(kegg_c2) && nrow(kegg_c2) > 0) kegg_c2@result %>% mutate(Cluster = "C2"),
  if(exists("kegg_c3") && !is.null(kegg_c3) && nrow(kegg_c3) > 0) kegg_c3@result %>% mutate(Cluster = "C3")
) %>% select(Cluster, Description, Count, p.adjust)

my_colors <- c("C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35")

# ==============================================================================
# 2. 彻底修复图例重叠 + 进一步优化排版的 SCI 级和弦图函数
# ==============================================================================
plot_custom_chord <- function(enrich_df, output_path, top_n = 5, 
                              cluster_colors = c("C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35")) {
  
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  
  # 2.1 强制过滤并排序
  plot_data <- enrich_df %>%
    filter(p.adjust < 0.1) %>%
    mutate(neg_log10_p = -log10(p.adjust + 1e-10)) %>%
    group_by(Cluster) %>%
    arrange(desc(Count)) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  if(nrow(plot_data) == 0) {
    message("⚠️ 警告: 数据过滤后为空，无法绘制: ", basename(output_path))
    return(invisible(NULL))
  }
  
  # 2.2 建立唯一的 Pathway 图例字典 (含 Count & adjp)
  term_info <- plot_data %>%
    arrange(desc(Count)) %>% 
    distinct(Description, .keep_all = TRUE) %>%
    mutate(
      Short_ID = paste0("P", row_number()),
      Label_For_Legend = paste0(Short_ID, " : ", Description, 
                                " (Count=", Count, ", adjp=", sprintf("%.2e", p.adjust), ")")
    ) %>%
    select(Original_Name = Description, Short_ID, Label_For_Legend)
  
  # 2.3 安全合并回原始绘图数据
  plot_data <- plot_data %>%
    left_join(term_info, by = c("Description" = "Original_Name")) %>%
    rename(Original_Name = Description) %>% 
    mutate(Description = Short_ID)          
  
  csv_path <- sub("\\.pdf$", "_Mapping_Table.csv", output_path)
  write.csv(term_info, csv_path, row.names = FALSE)
  
  # ====================== 关键修复：画布与图例空间计算 ======================
  
  actual_clusters <- intersect(names(cluster_colors), unique(plot_data$Cluster))
  all_terms <- unique(plot_data$Description)
  
  # 设置颜色
  term_colors <- rep("grey85", length(all_terms))
  names(term_colors) <- all_terms
  grid_col <- c(cluster_colors[actual_clusters], term_colors)
  order_sectors <- c(actual_clusters, all_terms)
  
  min_p <- min(plot_data$neg_log10_p)
  max_p <- max(plot_data$neg_log10_p)
  col_fun <- colorRamp2(c(min_p, (min_p + max_p)/2, max_p), 
                        c("#85C6EA", "#F3BDB7", "#F2722A"))
  link_colors <- col_fun(plot_data$neg_log10_p)
  
  # ====================== 动态间隙（防止扇形标签重叠） ======================
  gap_sizes <- c(
    rep(1.5, length(actual_clusters) - 1), 15,   # Cluster 之间大间隔
    rep(1.5, length(all_terms) - 1), 15           # Pathway 之间大间隔
  )
  
  # ====================== 绘图 ======================
  # 关键修复：增大画布宽度，圆盘左移，右侧留出充足图例空间
  pdf(output_path, width = 18, height = 13)
  
  circos.clear()
  circos.par(
    start.degree = 180, 
    gap.after = gap_sizes,
    track.margin = c(0.01, 0.01),
    # 核心修复：canvas 左右不对称，左侧正常，右侧收窄，把圆盘推向左边，为图例让路
    canvas.xlim = c(-1.2, 0.7),
    canvas.ylim = c(-1.15, 1.15)
  )
  
  chordDiagram(plot_data[, c("Cluster", "Description", "Count")], 
               order = order_sectors, 
               grid.col = grid_col, 
               col = link_colors,
               transparency = 0.35, 
               annotationTrack = c("grid"), 
               preAllocateTracks = list(track.height = 0.18))
  
  # 扇形标签绘制
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    
    is_cluster <- sector.name %in% actual_clusters
    text_cex  <- ifelse(is_cluster, 1.25, 0.72)
    text_font <- ifelse(is_cluster, 2, 1)
    
    circos.text(mean(xlim), ylim[1] + 0.25, sector.name, 
                facing = "clockwise", niceFacing = TRUE, 
                adj = c(0, 0.5), cex = text_cex, font = text_font)
  }, bg.border = NA)
  
  # ====================== 右侧图例（关键修复） ======================
  lgd_p <- Legend(
    col_fun = col_fun, 
    title = expression(-log[10](p.adjust)), 
    direction = "horizontal",
    legend_width = unit(6, "cm"),
    title_position = "topcenter"
  )
  
  lgd_pathways <- Legend(
    labels = term_info$Label_For_Legend,   
    title = "Enriched Pathways (Count & adjp)",
    type = "points", 
    pch = 15,
    legend_gp = gpar(col = "grey85"),
    labels_gp = gpar(fontsize = 9.5, lineheight = 1.35),
    title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  lgd_combined <- packLegend(lgd_p, lgd_pathways, 
                             direction = "vertical", 
                             gap = unit(1.5, "cm"))
  
  # 图例绘制在右侧，与圆盘不重叠
  draw(lgd_combined, 
       x = unit(0.01, "npc"), 
       y = unit(0.5, "npc"), 
       just = c("left", "center"))
  
  dev.off()
  circos.clear()
  
  cat("✅ 成功保存修复版和弦图至:", output_path, "\n")
  cat("📁 映射表格已更新至:", csv_path, "\n\n")
}

# ==============================================================================
# 3. 执行绘图
# ==============================================================================
plot_custom_chord(
  enrich_df = go_res, 
  output_path = "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/Step4_Enrichment/GO_Enrichment_Chord_V2.pdf", 
  top_n = 5,       # GO 通路名称较长，可改为 4 
  cluster_colors = my_colors
)

plot_custom_chord(
  enrich_df = kegg_res, 
  output_path = "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/Step4_Enrichment/KEGG_Enrichment_Chord_V2.pdf", 
  top_n = 5, 
  cluster_colors = my_colors
)









# ==============================================================================
# Step 2-Ext: 外部数据集验证
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(GSVA)
  library(pheatmap)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})


# ==============================================================================
# 🎨 0.1 全局颜色配置区 (确保所有热图颜色绝对一致)
# ==============================================================================
GLOBAL_SUBTYPE_COL <- c("C1" = "#00A087", "C2" = "#00BFC4", "C3" = "#E64B35")

GLOBAL_SCORE_COL <- colorRamp2(c(-2, -1, 0, 1, 2), 
                               c("navy", "#4575B4", "white", "#D73027", "firebrick"))

GLOBAL_CLIN_COL <- list(
  "Sex"  = c("F" = "#FB8072", "M" = "#80B1D3"),
  "Size" = c("large" = "#FDB462", "small" = "#B3DE69"),
  "NF2"  = c("N" = "#D9D9D9", "Y" = "#BC80BD")
)

# ==============================================================================
# 🌟 0.2
# ==============================================================================
TOP_N_GENES <- 100   # 

# 👇 
DATASETS <- list(
  list(
    id        = "GSE141801",
    label     = "Gugel et al 2020",
    expr_path = file.path(BULK_CLEAN_DIR, "GSE141801_Gugel/GSE141801_Expression_Log2.csv"),
    clin_path = file.path(BULK_CLEAN_DIR, "GSE141801_Gugel/Clinical_Data.csv"),
    sample_col = "ID",   
    vars = list(
      list(key = "Sex",         label = "Sex",        type = "discrete"),
      list(key = "Age",         label = "Age",        type = "continuous"),
      list(key = "NF2",         label = "NF2",        type = "discrete"),
      list(key = "Size",        label = "Size",       type = "discrete"),
      list(key = "Grade",       label = "Grade",      type = "discrete"),
      list(key = "Recurrence",  label = "Recur",      type = "discrete")
    )
  ),
  list(
    id        = "GSE39645",
    label     = "Torres-Martin et al 2013",
    expr_path = file.path(BULK_CLEAN_DIR, "GSE39645_Torres/GSE39645_Expression_Log2.csv"),
    clin_path = file.path(BULK_CLEAN_DIR, "GSE39645_Torres/Clinical_Data.csv"),
    sample_col = "ID",   
    vars = list(
      list(key = "Sex",         label = "Sex",        type = "discrete"),
      list(key = "Age",         label = "Age",        type = "continuous"),
      list(key = "NF2",         label = "NF2",        type = "discrete"),
      list(key = "Size",        label = "Size",       type = "discrete"),
      list(key = "Grade",       label = "Grade",      type = "discrete"),
      list(key = "Recurrence",  label = "Recur",      type = "discrete")
    )
  )
)


# ==============================================================================
load(file.path(RDATA_DIR, "Step3_DEA_Results.RData"))

get_top_sig <- function(deg_df, n) {
  deg_df %>%
    filter(adj.P.Val < 0.05 & logFC > 1) %>% 
    arrange(desc(logFC)) %>%
    head(n) %>%
    pull(Gene)
}

gene_sets_raw <- list(
  C1_Proliferative = get_top_sig(deg_c1, TOP_N_GENES),
  C2_Mesenchymal   = get_top_sig(deg_c2, TOP_N_GENES),
  C3_Immune        = get_top_sig(deg_c3, TOP_N_GENES)
)

overlap_12 <- intersect(gene_sets_raw$C1_Proliferative, gene_sets_raw$C2_Mesenchymal)
overlap_13 <- intersect(gene_sets_raw$C1_Proliferative, gene_sets_raw$C3_Immune)
overlap_23 <- intersect(gene_sets_raw$C2_Mesenchymal,   gene_sets_raw$C3_Immune)
all_overlap <- unique(c(overlap_12, overlap_13, overlap_23))

if (length(all_overlap) > 0) {
  message("⚠️ 发现 ", length(all_overlap), " 个重叠基因，正在移除...")
  gene_sets <- lapply(gene_sets_raw, function(x) setdiff(x, all_overlap))
} else {
  gene_sets <- gene_sets_raw
}

message(sprintf("✅ 基因签名准备完毕 (目标Top %d):", TOP_N_GENES))
for (nm in names(gene_sets)) {
  message("  ", nm, ": ", length(gene_sets[[nm]]), " 个特异性基因")
}

# ==============================================================================
# 2. 辅助函数定义
# ==============================================================================
load_expr <- function(filepath) {
  raw <- read.csv(filepath, check.names = FALSE)
  if (is.character(raw[, 1])) {
    mat <- as.matrix(raw[, -1])
    rownames(mat) <- raw[, 1]
  } else {
    mat <- as.matrix(raw)
  }
  if (any(duplicated(rownames(mat)))) mat <- limma::avereps(mat)
  mat
}

classify_samples <- function(expr_mat, gene_sets) {
  if (packageVersion("GSVA") >= "1.50.0") {
    raw_scores <- gsva(ssgseaParam(expr_mat, gene_sets), verbose = FALSE)
  } else {
    raw_scores <- gsva(expr_mat, gene_sets, method = "ssgsea", kcdf = "Gaussian")
  }
  scaled    <- t(scale(t(raw_scores)))
  predicted <- apply(scaled, 2, function(x) names(x)[which.max(x)])
  list(scaled_scores = scaled, predicted = predicted)
}

load_clin <- function(clin_path, sample_col) {
  if (!file.exists(clin_path)) return(NULL)
  df <- read.csv(clin_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (sample_col == "rownames") {
    rownames(df) <- df[, 1]
    df <- df[, -1, drop = FALSE]
  } else if (sample_col %in% colnames(df)) {
    rownames(df) <- df[[sample_col]]
    df <- df[, colnames(df) != sample_col, drop = FALSE]
  }
  df
}

short_subtype <- function(x) {
  x <- gsub("C1_Proliferative", "C1", x)
  x <- gsub("C2_Mesenchymal",   "C2", x)
  x <- gsub("C3_Immune",        "C3", x)
  x
}

make_discrete_colors <- function(vals) {
  lvls <- sort(unique(na.omit(as.character(vals))))
  n    <- length(lvls)
  pal  <- if (n <= 8) brewer.pal(max(3, n), "Set2")[seq_len(n)] else colorRampPalette(brewer.pal(8, "Set2"))(n)
  setNames(pal, lvls)
}

make_continuous_col <- function(vals) {
  rng <- range(vals, na.rm = TRUE)
  colorRamp2(c(rng[1], mean(rng), rng[2]), c("#2166AC", "white", "#D6604D"))
}

# 输出目录 (修正了路径拼接方式，防止 sprintf 报错)
out_dir_single <- file.path(FIG_DIR, paste0("Step2_Ext_Validation_Top", TOP_N_GENES, "/Single"))
out_dir_merged <- file.path(FIG_DIR, paste0("Step2_Ext_Validation_Top", TOP_N_GENES, "/Merged"))
dir.create(out_dir_single, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_merged, showWarnings = FALSE, recursive = TRUE)

all_scores    <- list()
all_predicted <- list()
all_study     <- list()
all_clin      <- list()

# ==============================================================================
# 🎯 模块一：分开验证每个数据集
# ==============================================================================
message("\n================ 开始单队列独立验证 ================")

for (ds in DATASETS) {
  if (!file.exists(ds$expr_path)) {
    warning("⚠️ 路径不存在跳过: ", ds$expr_path)
    next
  }
  message("\n📊 正在处理单队列: ", ds$id)
  
  # 1. 表达矩阵与打分
  expr_mat <- load_expr(ds$expr_path)
  res      <- classify_samples(expr_mat, gene_sets)
  clin_df  <- load_clin(ds$clin_path, ds$sample_col)
  
  # 2. 收集供合并使用
  all_scores[[ds$id]]    <- res$scaled_scores
  all_predicted[[ds$id]] <- res$predicted
  all_study[[ds$id]]     <- rep(ds$label, ncol(res$scaled_scores))
  all_clin[[ds$id]]      <- clin_df
  
  # 3. 排序 (强制 C1 -> C2 -> C3)
  pred_short <- short_subtype(res$predicted)
  sort_idx <- order(factor(pred_short, levels = c("C1", "C2", "C3")))
  
  plot_mat   <- res$scaled_scores[, sort_idx]
  plot_pred  <- pred_short[sort_idx]
  plot_names <- colnames(plot_mat)
  
  plot_mat[plot_mat >  2] <-  2
  plot_mat[plot_mat < -2] <- -2
  rownames(plot_mat) <- c("C1", "C2", "C3")
  
  # 4. 构建单队列临床注释
  anno_df <- data.frame(Predicted = plot_pred, row.names = plot_names, stringsAsFactors = FALSE)
  color_list <- list(Predicted = GLOBAL_SUBTYPE_COL)
  anno_labels_vec <- c(Predicted = "Predicted")
  
  if (!is.null(clin_df)) {
    for (vc in ds$vars) {
      if (vc$key %in% colnames(clin_df)) {
        vals <- as.character(clin_df[plot_names, vc$key])
        if (vc$type == "continuous") {
          num_vals <- suppressWarnings(as.numeric(vals))
          if (all(is.na(num_vals))) next
          anno_df[[vc$key]] <- num_vals
          color_list[[vc$key]] <- make_continuous_col(num_vals)
        } else {
          if (all(is.na(vals))) next
          if (vc$key %in% names(GLOBAL_CLIN_COL)) {
            anno_df[[vc$key]] <- factor(vals, levels = names(GLOBAL_CLIN_COL[[vc$key]]))
            color_list[[vc$key]] <- GLOBAL_CLIN_COL[[vc$key]]
          } else {
            anno_df[[vc$key]] <- factor(vals)
            color_list[[vc$key]] <- make_discrete_colors(vals)
          }
        }
        anno_labels_vec[[vc$key]] <- vc$label
      }
    }
  }
  
  ha_single <- HeatmapAnnotation(
    df      = anno_df,
    col     = color_list,
    na_col  = "grey88",
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
    annotation_label     = anno_labels_vec,
    simple_anno_size     = unit(4, "mm"),
    gap     = unit(1, "mm"),
    border  = TRUE,
    gp      = gpar(col = "white", lwd = 0.5)
  )
  
  # 5. 绘图
  subtype_splits <- factor(plot_pred, levels = c("C1", "C2", "C3"))
  
  ht_single <- Heatmap(
    plot_mat,
    name = "Signature\nScore",
    col  = GLOBAL_SCORE_COL,
    column_split    = subtype_splits,
    column_gap      = unit(2, "mm"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    top_annotation  = ha_single,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    row_names_side  = "left",
    show_column_names = FALSE,
    border  = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    width  = unit(10, "cm"),
    height = unit(3,  "cm")
  )
  
  out_pdf <- file.path(out_dir_single, paste0("Score_Heatmap_", ds$id, ".pdf"))
  pdf_h <- 3 + length(anno_labels_vec) * 0.4
  pdf(out_pdf, width = 10, height = pdf_h)
  draw(ht_single, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  message("  ✅ ", ds$id, " 单队列出图完成。")
}

# ==============================================================================
# 🌍 模块二：
# ==============================================================================
message("\n================ 开始多队列合并验证 ================")

combined_scores    <- do.call(cbind, all_scores)
combined_predicted <- short_subtype(unlist(all_predicted))
combined_study     <- unlist(all_study)

study_order <- sapply(DATASETS, `[[`, "label")
sort_order  <- order(
  factor(combined_predicted, levels = c("C1", "C2", "C3")),
  factor(combined_study, levels = study_order)
)

combined_scores    <- combined_scores[, sort_order]
combined_predicted <- combined_predicted[sort_order]
combined_study     <- combined_study[sort_order]
all_sample_names   <- colnames(combined_scores)

combined_scores[combined_scores >  2] <-  2
combined_scores[combined_scores < -2] <- -2
rownames(combined_scores) <- c("C1", "C2", "C3")

message("✅ 数据集合并完毕，总样本数: ", ncol(combined_scores))

anno_rows <- list()
var_meta <- list()
for (ds in DATASETS) {
  for (vc in ds$vars) {
    if (is.null(var_meta[[vc$key]])) {
      var_meta[[vc$key]] <- list(label = vc$label, type = vc$type)
    }
  }
}

for (vkey in names(var_meta)) {
  vtype  <- var_meta[[vkey]]$type
  vlabel <- var_meta[[vkey]]$label
  
  vals <- rep(NA_character_, length(all_sample_names))
  names(vals) <- all_sample_names
  
  for (ds in DATASETS) {
    clin_df <- all_clin[[ds$id]]
    if (!is.null(clin_df) && vkey %in% colnames(clin_df)) {
      matched <- intersect(names(vals), rownames(clin_df))
      vals[matched] <- as.character(clin_df[matched, vkey])
    }
  }
  
  if (vtype == "continuous") {
    num_vals <- suppressWarnings(as.numeric(vals))
    if (all(is.na(num_vals))) next
    anno_rows[[vkey]] <- list(vals = num_vals, label = vlabel, type = "continuous", colors = make_continuous_col(num_vals))
  } else {
    if (all(is.na(vals))) next
    anno_rows[[vkey]] <- list(vals = vals, label = vlabel, type = "discrete", colors = make_discrete_colors(vals))
  }
}

# 动态提取剩余两个数据集的颜色标签
study_col <- setNames(
  c("#6E9BF8", "#66C2A5")[1:length(DATASETS)],
  sapply(DATASETS, `[[`, "label")
)

anno_df    <- data.frame(Study = combined_study, Predicted = combined_predicted, row.names = all_sample_names, stringsAsFactors = FALSE)
color_list <- list(Study = study_col, Predicted = GLOBAL_SUBTYPE_COL)
anno_labels_vec <- c(Study = "Study", Predicted = "Predicted")

for (rk in names(anno_rows)) {
  ar <- anno_rows[[rk]]
  if (ar$type == "discrete") {
    anno_df[[rk]] <- factor(ar$vals, levels = names(ar$colors))
  } else {
    anno_df[[rk]] <- ar$vals
  }
  color_list[[rk]]      <- ar$colors
  anno_labels_vec[[rk]] <- ar$label
}

# 覆盖为全局临床颜色
if ("Sex" %in% names(color_list))  color_list[["Sex"]]  <- GLOBAL_CLIN_COL[["Sex"]] 
if ("Size" %in% names(color_list)) color_list[["Size"]] <- GLOBAL_CLIN_COL[["Size"]] 
if ("NF2" %in% names(color_list))  color_list[["NF2"]]  <- GLOBAL_CLIN_COL[["NF2"]] 

ha_top <- HeatmapAnnotation(
  df      = anno_df,
  col     = color_list,
  na_col  = "grey88",
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
  annotation_label     = anno_labels_vec,
  simple_anno_size     = unit(4, "mm"),
  gap     = unit(1, "mm"),
  border  = TRUE,
  gp      = gpar(col = "white", lwd = 0.5),
  annotation_legend_param = setNames(
    lapply(names(anno_labels_vec), function(k) {
      list(title = anno_labels_vec[[k]], title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7), ncol = 1)
    }),
    names(anno_labels_vec)
  )
)

subtype_splits <- factor(combined_predicted, levels = c("C1", "C2", "C3"))

ht <- Heatmap(
  combined_scores,
  name = "Signature\nScore",
  col  = GLOBAL_SCORE_COL,
  column_split    = subtype_splits,
  column_gap      = unit(2, "mm"),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  top_annotation  = ha_top,
  row_names_side  = "left",
  row_names_gp    = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  border            = TRUE,
  rect_gp           = gpar(col = "white", lwd = 0.5),  
  heatmap_legend_param = list(
    title         = "Signature\nScore",
    title_gp      = gpar(fontsize = 10, fontface = "bold"),
    labels_gp     = gpar(fontsize = 9),
    legend_height = unit(3, "cm"),
    at            = c(-2, -1, 0, 1, 2),
    border        = TRUE
  ),
  width  = unit(16, "cm"),
  height = unit(5,  "cm")
)

# 保存合并 PDF
out_file <- file.path(out_dir_merged, sprintf("Combined_Validation_Top%d.pdf", TOP_N_GENES))
n_anno  <- length(anno_labels_vec)
pdf_h   <- 3.5 + n_anno * 0.40   

pdf(out_file, width = 14, height = pdf_h)
draw(ht,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 35, 5, 5), "mm"))
dev.off()

message("✅ 联合验证热图已保存: ", out_file, " (", round(file.size(out_file)/1024), " KB)")
message("🎉 全部运行完毕！双队列的单项与合并结果已成功生成。")





# ==============================================================================
# Step 2.5-Ext:  (Hallmark Pathway Enrichment)
# ==============================================================================

suppressPackageStartupMessages({
  library(msigdbr)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
})

message("\n==============================")

# ==============================================================================
# 1. 准备 MSigDB Hallmark 基因集
# ==============================================================================

m_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)

# 
names(hallmark_list) <- gsub("HALLMARK_", "", names(hallmark_list))
message("✅ Hallmark 基因集获取成功 (共 ", length(hallmark_list), " 个通路)")

# ==============================================================================
# 2. 对每个外部数据集进行通路 ssGSEA 打分
# ==============================================================================
all_hw_scores <- list()

for (ds in DATASETS) {
  if (!file.exists(ds$expr_path)) next
  message("📊 正在对 ", ds$id, " 进行通路打分...")
  
  expr_mat <- load_expr(ds$expr_path)
  
  if (packageVersion("GSVA") >= "1.50.0") {
    hw_raw <- gsva(ssgseaParam(expr_mat, hallmark_list), verbose = FALSE)
  } else {
    hw_raw <- gsva(expr_mat, hallmark_list, method = "ssgsea", kcdf = "Gaussian")
  }

  hw_scaled <- t(scale(t(hw_raw)))
  all_hw_scores[[ds$id]] <- hw_scaled
  
  message("  ✅ ", ds$id, " 通路打分完成")
}

combined_hw_scores <- do.call(cbind, all_hw_scores)

combined_hw_scores <- combined_hw_scores[, sort_order]

message("✅ 通路得分矩阵合并并排序完毕，维度: ", 
        nrow(combined_hw_scores), " 通路 × ", ncol(combined_hw_scores), " 样本")

# ==============================================================================
# 4.  limma 
# ==============================================================================
message("🔍 正在计算各亚型特异性高表达通路...")


design <- model.matrix(~ 0 + factor(combined_predicted) + factor(combined_study))
colnames(design) <- c("C1", "C2", "C3", "Study_Covariate")

fit <- lmFit(combined_hw_scores, design)


contrast.matrix <- makeContrasts(
  C1_vs_Rest = C1 - (C2 + C3)/2,
  C2_vs_Rest = C2 - (C1 + C3)/2,
  C3_vs_Rest = C3 - (C1 + C2)/2,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


top_pathways <- c()
pathway_anno <- character() 

for (sub in c("C1", "C2", "C3")) {
  coef_name <- paste0(sub, "_vs_Rest")
  res <- topTable(fit2, coef = coef_name, number = Inf) %>%
    filter(adj.P.Val < 0.05 & logFC > 0) %>%
    arrange(desc(logFC)) %>%
    head(7) 
  sig_paths <- rownames(res)
  top_pathways <- c(top_pathways, sig_paths)
  pathway_anno <- c(pathway_anno, rep(sub, length(sig_paths)))
}


top_pathways <- unique(top_pathways)
plot_hw_mat <- combined_hw_scores[top_pathways, ]


plot_hw_mat[plot_hw_mat > 2]  <- 2
plot_hw_mat[plot_hw_mat < -2] <- -2


pathway_col <- colorRamp2(c(-2, 0, 2), c("#7B3294", "white", "#008837"))

pathway_split <- factor(pathway_anno, levels = c("C1", "C2", "C3"))


ht_pathway <- Heatmap(
  plot_hw_mat,
  name = "Pathway\nEnrichment",
  col  = pathway_col,
  

  column_split    = subtype_splits,
  column_gap      = unit(2, "mm"),
  top_annotation  = ha_top, 
  show_column_names = FALSE,
  cluster_columns = FALSE,
  
  # 行（通路）分割
  row_split       = pathway_split,
  row_gap         = unit(2, "mm"),
  row_title_gp    = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows    = FALSE, 
  
  row_names_side  = "left",
  row_names_gp    = gpar(fontsize = 10, fontface = "italic"), #
  

  border          = TRUE,
  rect_gp         = gpar(col = "white", lwd = 0.5),
  
  heatmap_legend_param = list(
    title         = "GSVA\nZ-score",
    title_gp      = gpar(fontsize = 10, fontface = "bold"),
    labels_gp     = gpar(fontsize = 9),
    legend_height = unit(3, "cm"),
    border        = TRUE
  )
)

# 
out_dir_merged <- file.path(FIG_DIR, paste0("Step2_Ext_Validation_Top", TOP_N_GENES, "/Merged"))
out_file_pw <- file.path(out_dir_merged, sprintf("Combined_Pathway_Enrichment_Top%d.pdf", TOP_N_GENES))


pdf_h_pw <- 3.5 + n_anno * 0.40 + length(top_pathways) * 0.15

pdf(out_file_pw, width = 14, height = pdf_h_pw)
draw(ht_pathway,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 35, 5, 5), "mm"))
dev.off()

