# ==============================================================================
# 06_SingleCell
# ==============================================================================
source("00_single cell config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
  library(reshape2)
  library(ggrepel)
  library(tibble)
})

set.seed(PARAMS$global_seed)
SC_FIG_DIR <- file.path(FIG_DIR, "Step06_SingleCell")
if (!dir.exists(SC_FIG_DIR)) dir.create(SC_FIG_DIR, recursive = TRUE)

# ─── 1. 加载数据 ──────────────────────────────────────────────────────────────
message("⏳ 加载单细胞数据...")
stopifnot(file.exists(SC_QS_PATH))
sc_new <- qread(SC_QS_PATH)
message("✅ sc_new loaded. Cells=", ncol(sc_new), " Genes=", nrow(sc_new))

# 加载 Bulk DEG
load(file.path(RDATA_DIR, "Step3_DEA_Results.RData"))
stopifnot(exists("deg_c1"), exists("deg_c2"), exists("deg_c3"))
message("✅ Step3 DEG loaded.")

# ─── 2. 修复 UMAP reduction ──────────────────────────────────────────────────
if (!"umap" %in% Reductions(sc_new)) {
  if (all(c("UMAP_1", "UMAP_2") %in% colnames(sc_new@meta.data))) {
    umap_coords <- as.matrix(sc_new@meta.data[, c("UMAP_1", "UMAP_2")])
    sc_new[["umap"]] <- CreateDimReducObject(
      embeddings = umap_coords, key = "UMAP_",
      assay = DefaultAssay(sc_new)
    )
    message("✅ UMAP reduction 已从 meta.data 修复")
  } else {
    stop("❌ 无法找到 UMAP 坐标，请检查数据")
  }
}

# ─── 3. 设置细胞注释 ─────────────────────────────────────────────────────────
if ("final_label" %in% colnames(sc_new@meta.data)) {
  Idents(sc_new) <- "final_label"
  message("✅ Idents → final_label")
} else {
  message("⚠️ 未找到 final_label, 使用默认 Idents")
}

# 确保 CELLTYPE_COLORS 覆盖所有 ident
ct_levels <- levels(Idents(sc_new))
CELLTYPE_COLORS <- complete_palette(ct_levels, CELLTYPE_COLORS)

# ─── 4. 提取特征基因并打分（使用与 Bulk 一致的 Top100 严格筛选）──────────────
message("🚀 提取 Bulk Top100 特征基因...")

gene_sets <- build_bulk_gene_sets(deg_c1, deg_c2, deg_c3, n = PARAMS$sc_n_sig_genes)
message("📊 Raw gene counts: ",
        paste(names(gene_sets), sapply(gene_sets, length), sep = "=", collapse = " / "))

# 大小写匹配到单细胞行名
sc_genes <- rownames(sc_new)
sig_list_final <- lapply(gene_sets, function(g) {
  m_idx <- match(toupper(g), toupper(sc_genes))
  unique(sc_genes[m_idx[!is.na(m_idx)]])
})

# C1 补救：如果匹配基因过少，注入核心增殖基因
if (length(sig_list_final$C1_Proliferative) < 5) {
  message("⚠️ C1 匹配基因不足 5 个，注入增殖基因补救")
  rescue <- c("MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "MCM6", "CCNA2")
  m_idx <- match(toupper(rescue), toupper(sc_genes))
  sig_list_final$C1_Proliferative <- unique(c(
    sig_list_final$C1_Proliferative, sc_genes[m_idx[!is.na(m_idx)]]
  ))
}

message("✅ Matched gene counts: ",
        paste(names(sig_list_final), sapply(sig_list_final, length),
              sep = "=", collapse = " / "))
if (any(sapply(sig_list_final, length) < 10)) {
  warning("⚠️ 部分 signature <10 基因, module score 可能不稳定")
}

# 清理旧列 → 打分
old_cols <- grep("Score_C|Bulk_Score", colnames(sc_new@meta.data), value = TRUE)
if (length(old_cols) > 0) {
  sc_new@meta.data <- sc_new@meta.data[, setdiff(colnames(sc_new@meta.data), old_cols),
                                       drop = FALSE]
}

sc_new <- AddModuleScore(sc_new, features = sig_list_final,
                         name = "Bulk_Score", seed = PARAMS$global_seed)
# AddModuleScore 按 list 顺序生成 Bulk_Score1/2/3
sc_new$Score_C1 <- sc_new$Bulk_Score1
sc_new$Score_C2 <- sc_new$Bulk_Score2
sc_new$Score_C3 <- sc_new$Bulk_Score3

# 归类
score_mat <- sc_new@meta.data[, c("Score_C1", "Score_C2", "Score_C3")]
sc_new$Bulk_Class <- factor(
  c("C1", "C2", "C3")[apply(score_mat, 1, which.max)],
  levels = c("C1", "C2", "C3")
)

# ─── 2.5 Seurat v5 兼容性修复 ────────────────────────────────────────────────
# 更新旧版对象结构
if (inherits(try(UpdateSeuratObject(sc_new), silent = TRUE), "Seurat")) {
  sc_new <- UpdateSeuratObject(sc_new)
  message("✅ SeuratObject 已更新到 v5 格式")
}

# Seurat v5 中 layers 可能是按样本分裂的，需要 JoinLayers 后才能正常计算
DefaultAssay(sc_new) <- "RNA"
if (packageVersion("SeuratObject") >= "5.0.0") {
  # 检查是否有分裂的 layers
  all_layers <- Layers(sc_new[["RNA"]])
  message("📋 当前 RNA assay layers: ", paste(all_layers, collapse = ", "))
  
  # 如果存在多个 counts layer（如 counts.S1, counts.S2...），需要合并
  if (length(grep("^counts", all_layers)) > 1) {
    sc_new[["RNA"]] <- JoinLayers(sc_new[["RNA"]])
    message("✅ JoinLayers 完成，已合并分裂的 layers")
  }
}


if (!"percent.mt" %in% colnames(sc_new@meta.data)) {
  # 方案 A：尝试用 Seurat 原生函数
  tryCatch({
    sc_new[["percent.mt"]] <- PercentageFeatureSet(sc_new, pattern = "^MT-")
    message("✅ percent.mt 通过 PercentageFeatureSet 计算")
  }, error = function(e) {
    # 方案 B：手动计算（完全绕过 GetAssayData slot 问题）
    message("⚠️ PercentageFeatureSet 失败，手动计算 percent.mt...")
    if (packageVersion("SeuratObject") >= "5.0.0") {
      counts_mat <- LayerData(sc_new, assay = "RNA", layer = "counts")
    } else {
      counts_mat <- GetAssayData(sc_new, assay = "RNA", slot = "counts")
    }
    mt_genes <- grep("^MT-", rownames(counts_mat), value = TRUE)
    if (length(mt_genes) == 0) {
      mt_genes <- grep("^mt-", rownames(counts_mat), value = TRUE)
    }
    if (length(mt_genes) > 0) {
      sc_new$percent.mt <<- colSums(counts_mat[mt_genes, , drop = FALSE]) /
        colSums(counts_mat) * 100
      message("✅ percent.mt 手动计算完成，MT 基因数: ", length(mt_genes))
    } else {
      sc_new$percent.mt <<- 0
      warning("⚠️ 未找到 MT- 基因，percent.mt 设为 0")
    }
  })
}
p_qc1 <- FeatureScatter(sc_new, "nCount_RNA", "Score_C1") + ggtitle("nCount vs C1")
p_qc2 <- FeatureScatter(sc_new, "nFeature_RNA", "Score_C1") + ggtitle("nFeature vs C1")
p_qc3 <- FeatureScatter(sc_new, "percent.mt", "Score_C1") + ggtitle("percent.mt vs C1")
p_qc <- p_qc1 | p_qc2 | p_qc3
print(p_qc)
save_pdf("QC_TechBias_ScoreC1.pdf", p_qc, width = 12, height = 4, subdir = "Step06_SingleCell")

# ─── 6. 图 1: UMAP（半矢量图）──────────────────────────────────────────────
message("🎨 绘制 UMAP...")

umap_df <- as.data.frame(Embeddings(sc_new, "umap"))
umap_df$cluster <- as.character(Idents(sc_new))
set.seed(PARAMS$global_seed)
umap_df <- umap_df[sample(nrow(umap_df)), ]

label_df <- umap_df %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.1, alpha = 0.8, shape = 16, stroke = 0) +
  scale_color_manual(values = CELLTYPE_COLORS) +
  geom_text_repel(data = label_df, aes(label = cluster),
                  color = "black", size = 4, fontface = "bold",
                  max.overlaps = 20) +
  theme_void() +
  ggtitle("Single-Cell Landscape") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# 保存矢量版 PDF（文件大，适合编辑）
save_pdf("Fig06A_UMAP_vector.pdf", p_umap, width = 5.5, height = 5.5,
         subdir = "Step06_SingleCell")

# 保存半矢量版 PDF（文件小，适合投稿）
save_semi_vector("Fig06A_UMAP_semi.pdf", p_umap, width = 5.5, height = 5.5,
                 dpi = 300, subdir = "Step06_SingleCell")
───────────────────────────
message("🎨 绘制 FeaturePlot...")
target_names <- c("Score_C1", "Score_C2", "Score_C3")

p_feat_list <- FeaturePlot(sc_new, features = target_names, combine = FALSE,
                           order = TRUE, min.cutoff = "q10", max.cutoff = "q95",
                           raster = FALSE)
p_feat_list <- lapply(p_feat_list, function(x) {
  x + scale_color_gradientn(colors = c("#E0E0E0", "firebrick3")) +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
})
p_feat_combined <- wrap_plots(p_feat_list, ncol = 3)
save_half_vector("Fig06B_FeaturePlot.pdf", p_feat_combined, width = 12, height = 4,
                 subdir = "Step06_SingleCell")





# 确保安装了 ggrastr 包（如果没有，请先运行 install.packages("ggrastr")）
library(ggrastr)

message("🎨 正在生成真正的 高清半矢量 + 无描边 FeaturePlot...")

# 遍历已经画好的图表列表
p_feat_list_raster <- lapply(p_feat_list, function(p) {
  
  # 1. 第一步：依然去描边，保证像素纯净
  p$layers[[1]]$aes_params$stroke <- 0
  p$layers[[1]]$aes_params$shape <- 16
  
  # 2. 第二步 (核心)：使用 ggrastr 将散点图层强行“栅格化” (转为 600 dpi 的位图)
  p <- ggrastr::rasterise(p, layers = "Point", dpi = 600)
  
  return(p)
})

# 重新用 patchwork 拼接
p_feat_combined_raster <- wrap_plots(p_feat_list_raster, ncol = 3)

# 现在你可以放心保存了，出来的绝对是中间像素点、外围矢量的标准半矢量图
# 注意：使用 ggrastr 后，直接用普通 ggsave 或 pdf() 保存即可生效
ggsave(file.path("Step06_SingleCell", "Fig06B_FeaturePlot_TrueRaster.pdf"), 
       plot = p_feat_combined_raster, 
       width = 12, height = 4)

message("✅ 完美半矢量版本已保存至: Fig06B_FeaturePlot_TrueRaster.pdf")
# ─── 8. 图 3: Violin + Boxplot（缩尾处理）───────────────────────────────────
message("🎨 绘制 Violin Plot...")
score_colors <- c("Score_C1" = SUBTYPE_COLORS["C1"],
                  "Score_C2" = SUBTYPE_COLORS["C2"],
                  "Score_C3" = SUBTYPE_COLORS["C3"])

plot_vln <- FetchData(sc_new, vars = c("ident", target_names)) %>%
  reshape2::melt(id.vars = "ident", variable.name = "Bulk_Class", value.name = "Score")

upper_lim <- quantile(plot_vln$Score, 0.99, na.rm = TRUE)
lower_lim <- quantile(plot_vln$Score, 0.01, na.rm = TRUE)

p_vln <- ggplot(plot_vln, aes(x = ident, y = Score, fill = Bulk_Class)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.8, linewidth = 0.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", alpha = 0.5) +
  scale_fill_manual(values = score_colors) +
  facet_grid(Bulk_Class ~ ., scales = "free_y") +
  coord_cartesian(ylim = c(lower_lim, upper_lim)) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = "Module Score", title = "Bulk Signature Enrichment")

save_pdf("Fig06C_Violin.pdf", p_vln, width = 8, height = 8, subdir = "Step06_SingleCell")

# ─── 9. 图 4: Highlight UMAP（Top 5% 高分细胞按 cluster 原色点亮）──────────
message("🎨 绘制 Highlight UMAP...")

plot_highlight <- function(obj, feature_col, title, quantile_thr, colors) {
  df <- FetchData(obj, vars = c("UMAP_1", "UMAP_2", "ident", feature_col))
  colnames(df) <- c("UMAP_1", "UMAP_2", "Cluster", "Score")
  df$Cluster <- factor(df$Cluster, levels = names(colors))
  
  cut_val <- quantile(df$Score, probs = quantile_thr, na.rm = TRUE)
  bg <- df[df$Score <= cut_val, ]
  fg <- df[df$Score >  cut_val, ]
  
  ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2), color = "grey92", size = 0.01) +
    geom_point(data = fg, aes(UMAP_1, UMAP_2, color = Cluster), size = 0.05) +
    scale_color_manual(values = colors, drop = FALSE) +
    theme_void() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none")
}

q_thr <- PARAMS$sc_highlight_quantile
p_h1 <- plot_highlight(sc_new, "Score_C1", "C1 High (Top 5%)", q_thr, CELLTYPE_COLORS)
p_h2 <- plot_highlight(sc_new, "Score_C2", "C2 High (Top 5%)", q_thr, CELLTYPE_COLORS)
p_h3 <- plot_highlight(sc_new, "Score_C3", "C3 High (Top 5%)", q_thr, CELLTYPE_COLORS)

# 生成图例
p_leg <- DimPlot(sc_new, reduction = "umap", cols = CELLTYPE_COLORS) +
  theme(legend.position = "right")
leg_part <- cowplot::get_legend(p_leg)

p_highlight <- plot_grid(
  plot_grid(p_h1, p_h2, p_h3, nrow = 1),
  leg_part, ncol = 2, rel_widths = c(3, 0.6)
)
save_half_vector("Fig06D_Highlight.pdf", p_highlight, width = 15, height = 5,
                 subdir = "Step06_SingleCell")

# 单独保存
save_half_vector("Fig06D_C1_Highlight.pdf", p_h1, width = 5, height = 5,
                 subdir = "Step06_SingleCell")
save_half_vector("Fig06D_C2_Highlight.pdf", p_h2, width = 5, height = 5,
                 subdir = "Step06_SingleCell")
save_half_vector("Fig06D_C3_Highlight.pdf", p_h3, width = 5, height = 5,
                 subdir = "Step06_SingleCell")

# ─── 10. 图 5: Bubble Plot（气泡图）────────────────────────────────────────
message("🎨 绘制 Bubble Plot...")
meta_bub <- sc_new@meta.data
meta_bub$ident_chr <- as.character(Idents(sc_new))

dot_data <- meta_bub %>%
  dplyr::select(ident = ident_chr, Score_C1, Score_C2, Score_C3) %>%
  pivot_longer(cols = starts_with("Score_C"), names_to = "Bulk_Group", values_to = "Score") %>%
  group_by(ident, Bulk_Group) %>%
  summarise(
    Avg_Score = mean(Score, na.rm = TRUE),
    Pct_Positive = mean(Score > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# 颜色上限防极端值
max_thr <- quantile(dot_data$Avg_Score[dot_data$Avg_Score > 0], 0.90, na.rm = TRUE)
dot_data$Avg_Capped <- pmin(pmax(dot_data$Avg_Score, 0), max_thr)
dot_data$Bulk_Group <- factor(dot_data$Bulk_Group,
                              levels = c("Score_C1", "Score_C2", "Score_C3"))

p_bubble <- ggplot(dot_data, aes(x = Bulk_Group, y = ident)) +
  geom_point(aes(size = Pct_Positive, fill = Avg_Capped),
             shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_gradient(low = "grey92", high = "#D73027",
                      limits = c(0, max_thr), name = "Avg Score") +
  scale_size_continuous(range = c(1, 8), name = "% Score > 0") +
  scale_x_discrete(labels = c("C1", "C2", "C3")) +
  theme_publication() +
  theme(
    panel.grid.major = element_line(color = "grey92", linetype = "dashed"),
    axis.title = element_blank()
  )
save_pdf("Fig06E_Bubble.pdf", p_bubble, width = 6, height = 8, subdir = "Step06_SingleCell")

# ─── 11. 保存打分后的对象 ─────────────────────────────────────────────────
saveRDS(sc_new, file.path(RDATA_DIR, "Step06_sc_new_scored.rds"))
message("✅ 模块 06 完成! 输出目录: ", SC_FIG_DIR)

rm(umap_df, label_df, plot_vln, meta_bub, dot_data, bg, fg); gc()










# ==============================================================================
# 07_C3high_Myeloid_Analysis.R
# 功能: 提取C3-high myeloid细胞 → 深入分析
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
setwd("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/")
set.seed(123)
OUT_DIR <- file.path( "Figures/Step07_C3Myeloid")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
getwd()
# ─── 1. 加载数据并提取C3-high myeloid ────────────────────────────────────
message("⏳ 加载打分后的单细胞数据...")
sc_new <- readRDS(file.path( "RData/Step06_sc_new_scored.rds"))

# 识别myeloid细胞（根据你的celltype命名调整）
myeloid_pattern <- "myeloid|Myeloid|Macrophage|Monocyte|DC|Dendritic"
myeloid_cells <- grep(myeloid_pattern, levels(Idents(sc_new)), value = TRUE, ignore.case = TRUE)

if (length(myeloid_cells) == 0) {
  stop("❌ 未找到myeloid细胞类型，请检查细胞注释")
}
message("✅ 识别到myeloid类型: ", paste(myeloid_cells, collapse = ", "))

# 提取myeloid子集
sc_myeloid <- subset(sc_new, idents = myeloid_cells)
message("📊 Myeloid细胞数: ", ncol(sc_myeloid))

















# ─── 定义 C3-high（GMM 自动阈值）────────────────────────────────────────
library(mclust)
set.seed(123)

# 
gmm_fit     <- Mclust(sc_myeloid$Score_C3, G = 2)
gmm_means   <- sort(gmm_fit$parameters$mean)
c3_threshold <- mean(gmm_means)   # 两高斯分量均值的中点作为分界阈值

message("✅ GMM 自动阈值: ", round(c3_threshold, 4),
        "  (Component means: ", round(gmm_means[1], 4),
        " & ", round(gmm_means[2], 4), ")")

# 基于 GMM 阈值标注 C3_status
sc_myeloid$C3_status <- ifelse(sc_myeloid$Score_C3 > c3_threshold,
                               "C3-high", "C3-low")

sc_c3high <- subset(sc_myeloid, subset = C3_status == "C3-high")

message("✅ C3-high myeloid 细胞数: ", ncol(sc_c3high),
        "  (占比: ", round(ncol(sc_c3high) / ncol(sc_myeloid) * 100, 1), "%)")







# # ─── 2. 重新聚类（精细分型）──────────────────────────────────────────────
# message("🔬 对C3-high myeloid进行重新聚类...")
sc_c3high <- FindNeighbors(sc_c3high, dims = 1:30)
sc_c3high <- FindClusters(sc_c3high, resolution = 0.5)
sc_c3high <- RunUMAP(sc_c3high, dims = 1:30)
# message("🔬 对C3-high myeloid进行重新聚类...")

# 
sc_c3high <- FindVariableFeatures(sc_c3high, selection.method = "vst", nfeatures = 2000)

# 2. 重新缩放数据 (Scale)
sc_c3high <- ScaleData(sc_c3high)

# 3. 运行 PCA
sc_c3high <- RunPCA(sc_c3high, features = VariableFeatures(object = sc_c3high), verbose = FALSE)

# 4. 然后再运行你之前报错的代码及后续聚类
sc_c3high <- FindNeighbors(sc_c3high, dims = 1:30)
sc_c3high <- FindClusters(sc_c3high, resolution = 0.5) # 这里的 resolution 可以根据你想要的亚群颗粒度调整

# 5. 最后重新运行 UMAP 以便后续可视化
sc_c3high <- RunUMAP(sc_c3high, dims = 1:30)
p_umap_recluster <- DimPlot(sc_c3high, label = TRUE, pt.size = 0.5) +
  ggtitle("C3-high Myeloid Re-clustering") +
  theme_void()
# save_pdf("Fig07A_C3high_Reclustering.pdf", p_umap_recluster, 
#          width = 6, height = 5, subdir = "Step07_C3Myeloid")
# 1. 确保输出的子文件夹存在（这里假设你想保存在当前工作目录下的 Step07_C3Myeloid 文件夹）
out_dir <- "Step07_C3Myeloid"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. 用 ggsave 保存
ggsave(filename = file.path("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/Step07_C3Myeloid/Fig07A_C3high_Reclustering.pdf"), 
       plot = p_umap_recluster, 
       width = 5, height = 5)






# 确保加载必要的包
library(ggplot2)
library(mclust)
library(dplyr)

message("🎨 正在生成无歧义的“一刀切”硬分类 GMM 分布图...")

# 1. 提取打分数据
c3_scores <- sc_myeloid$Score_C3

# 2. 运行 GMM 获取分界点
set.seed(123)
gmm_fit <- Mclust(c3_scores, G = 2)
means <- sort(gmm_fit$parameters$mean)
gmm_cutoff <- mean(means)

# 3. 构建总拟合曲线数据
x_range <- seq(min(c3_scores), max(c3_scores), length.out = 500)
props <- gmm_fit$parameters$pro
vars  <- gmm_fit$parameters$variance$sigmasq
if (length(vars) == 1) vars <- rep(vars, 2)
if(gmm_fit$parameters$mean[1] > gmm_fit$parameters$mean[2]) {
  vars <- rev(vars)
  props <- rev(props)
}

gmm_curve <- data.frame(
  x = x_range,
  y_total = props[1] * dnorm(x_range, means[1], sqrt(vars[1])) +
    props[2] * dnorm(x_range, means[2], sqrt(vars[2]))
)

# 计算 Y 轴最大值，用于美化排版
max_y <- max(gmm_curve$y_total)

# 4. 绘制终极版“一刀切”图表
p_gmm_hard <- ggplot(data.frame(Score = c3_scores), aes(x = Score)) +
  # 原始数据的直方图
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "#E0E0E0", color = "white") +
  
  # GMM 总拟合曲线 (黑色实线)
  geom_line(data = gmm_curve, aes(x = x, y = y_total), color = "#2C3E50", linewidth = 1.2) +
  
  # 唯一且强势的分界线 (Cutoff)
  geom_vline(xintercept = gmm_cutoff, color = "#D73027", linewidth = 1.2, linetype = "dashed") +
  
  # 在绿线左侧标注 C3-low (使用冷色调表示非目标态)
  annotate("text", x = gmm_cutoff - 0.5, y = max_y * 0.5, 
           label = "C3-low", color = "#3182BD", fontface = "bold", size = 6, hjust = 1) +
  
  # 在绿线右侧标注 C3-high (使用暖色调表示目标高表达态)
  annotate("text", x = gmm_cutoff + 0.5, y = max_y * 0.5, 
           label = "C3-high", color = "#D73027", fontface = "bold", size = 6, hjust = 0) +
  
  # 标注 Cutoff 具体数值
  annotate("text", x = gmm_cutoff, y = max_y * 0.95, 
           label = paste0("Cutoff = ", round(gmm_cutoff, 3)), 
           color = "black", angle = 90, vjust = -0.8, hjust = 1, size = 4.5, fontface = "italic") +
  
  # 调整 X 轴显示范围 (为了让两边看着更对称些，可以根据你的实际数据调整)
  coord_cartesian(xlim = c(min(c3_scores), max(c3_scores)*0.8)) +
  
  # 主题与标签
  theme_classic(base_size = 14) +
  labs(
    title = "Identification of C3-high Myeloid Subpopulation",
    subtitle = "Hard classification based on GMM-derived objective threshold",
    x = "Module Score (Score_C3)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text = element_text(color = "black")
  )

# 5. 保存图表
out_dir <- "Figures/Step07_C3Myeloid"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "Fig07_GMM_Hard_Classification.pdf"), plot = p_gmm_hard, width = 5, height = 5)

message("✅ 图表已生成并保存至: ", file.path(out_dir, "Fig07_GMM_Hard_Classification.pdf"))









# ─── 3. 差异基因分析（C3-high vs C3-low）─────────────────────────────────
message("🧬 差异基因分析...")
library(Seurat)
Idents(sc_myeloid) <- "C3_status"
deg_c3 <- FindMarkers(sc_myeloid, ident.1 = "C3-high", ident.2 = "C3-low",
                      logfc.threshold = 0.25, min.pct = 0.1)
deg_c3$gene <- rownames(deg_c3)
deg_c3_sig <- deg_c3 %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)

write.csv(deg_c3_sig, file.path(OUT_DIR, "DEG_C3high_vs_C3low.csv"))
message("✅ 显著差异基因数: ", nrow(deg_c3_sig))

# # Top基因可视化
# library(Seurat)
# 
# # 如果报错说找不到 ggplot2，顺便也加载一下它（VlnPlot 底层依赖 ggplot2）
# library(ggplot2)
top_genes <- head(deg_c3_sig[order(-deg_c3_sig$avg_log2FC), ], 10)$gene


library(patchwork)

# 设置 combine = FALSE，返回一个包含 4 张图的列表
p_list <- VlnPlot(sc_myeloid, 
                  features = top_genes[1:10], 
                  group.by = "C3_status", 
                  combine = FALSE)

# 手动用 patchwork 拼图
p_violin <- wrap_plots(p_list, ncol = 2)
p_violin



# ─── 4. 通路富集分析 ──────────────────────────────────────────────────────
message("📊 通路富集分析...")
library(clusterProfiler)
library(org.Hs.eg.db)

# gene_list <- deg_c3_sig$avg_log2FC
# names(gene_list) <- deg_c3_sig$gene
# gene_list <- sort(gene_list, decreasing = TRUE)
# 
# # GSEA
# gsea_go <- gseGO(geneList = gene_list,OrgDb = org.Hs.eg.db,
#                  keyType = "SYMBOL",
#                  ont = "BP",
#                  pvalueCutoff = 0.05)
# 
# if (nrow(gsea_go@result) > 0) {

# }
# 1. 准备数据：使用你之前生成的全量结果 deg_c3 (不要用过滤后的 deg_c3_sig)
# 提取基因名和 Fold Change
gsea_data <- deg_c3

# 2. 数据清洗：剔除 NA 和 Inf 值
gsea_data <- gsea_data[!is.na(gsea_data$avg_log2FC) & is.finite(gsea_data$avg_log2FC), ]

# 3. 构建 GSEA 需要的 named vector (命名向量)
gene_list <- gsea_data$avg_log2FC
names(gene_list) <- gsea_data$gene # 如果你的基因名在 rownames，这里换成 rownames(gsea_data)

# 4. 排序：这是 GSEA 最关键的一步！必须从大到小严格降序排列
gene_list <- sort(gene_list, decreasing = TRUE)

# 💡 你可以运行这两行自检一下，输出必须都是 0 才安全：
# print(paste("无穷大数量:", sum(is.infinite(gene_list))))
# print(paste("NA数量:", sum(is.na(gene_list))))

message("🏃‍♂️ 正在重新运行 GSEA 分析...")

# 5. 重新运行 GSEA
library(clusterProfiler)
library(org.Hs.eg.db)

gsea_go <- gseGO(geneList     = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "SYMBOL", # 确认你的基因名是 SYMBOL (如 TP53) 而不是 ENSEMBL
                 ont          = "BP",     # Biological Process
                 pvalueCutoff = 0.05,
                 minGSSize    = 10,       # 过滤掉包含基因太少的通路
                 maxGSSize    = 500)      # 过滤掉包含基因太多的宽泛通路

message("✅ GSEA 运行完毕！")
p_gsea <- dotplot(gsea_go, showCategory = 20) +
  ggtitle("C3-high Myeloid Enriched Pathways")


save_pdf("Fig07C_GSEA_GO.pdf", p_gsea,
         width = 8, height = 6, subdir = "Step07_C3Myeloid")






# 确保加载了必要的包
library(dplyr)
library(ggplot2)
library(stringr)

message("🎨 正在绘制顶级期刊精简版 GSEA 棒棒糖图...")

# 1. 提取 GSEA 结果数据框
gsea_df <- as.data.frame(gsea_go)

# 2. 空间压缩策略：只提取最显著的 Top 4 上调 和 Top 4 下调通路
top_up <- gsea_df %>% 
  filter(NES > 0) %>% 
  arrange(p.adjust) %>% 
  slice_head(n = 5) %>% 
  mutate(Direction = "Up in C3-high")

top_down <- gsea_df %>% 
  filter(NES < 0) %>% 
  arrange(p.adjust) %>% 
  slice_head(n =10) %>% 
  mutate(Direction = "Down in C3-high")

# 合并数据，并计算 -log10(FDR)
plot_df <- bind_rows(top_up, top_down) %>%
  mutate(
    logFDR = -log10(p.adjust),
    # 针对过长的 GO 描述，强制在 35 个字符处换行，极限节省横向空间
    Description = str_wrap(Description, width = 35) 
  ) %>%
  # 严格按照 NES 从小到大排序，让图表呈完美的阶梯状
  arrange(NES) %>%
  mutate(Description = factor(Description, levels = Description))

# # 3. 绘制紧凑型高级棒棒糖图
# p_lollipop_compact <- ggplot(plot_df, aes(x = NES, y = Description)) +
#   
#   # 中心零线：作为上下调的绝对分水岭
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.6) +
#   
#   # 棒棒糖的“棍子”
#   geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), 
#                color = "grey70", linewidth = 1.2) +
#   
#   # 棒棒糖的“糖”：气泡
#   geom_point(aes(size = logFDR, fill = NES), 
#              shape = 21, color = "black", stroke = 0.6) +
#   
#   # 顶级色彩映射
#   scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, 
#                        name = "NES") +
#   
#   # 调整气泡大小范围
#   scale_size_continuous(range = c(3, 8), name = expression(-log[10](FDR))) +
#   
#   # 精简主题设置
#   theme_classic(base_size = 12) +
#   labs(
#     x = "Normalized Enrichment Score (NES)",
#     y = NULL, 
#     title = "Key GO Pathways in C3-high Myeloid"
#   ) +
#   theme(
#     plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
#     axis.text.y = element_text(color = "black", face = "bold", size = 10, lineheight = 0.8),
#     axis.text.x = element_text(color = "black"),
#     
#     # 🌟 修复点 1：强制使用 ggplot2::margin
#     axis.title.x = element_text(face = "bold", margin = ggplot2::margin(t = 10)),
#     
#     legend.position = "right",
#     legend.background = element_blank(),
#     
#     # 🌟 修复点 2：强制使用 ggplot2::margin
#     legend.box.margin = ggplot2::margin(0, 0, 0, -10), 
#     
#     panel.grid.major.y = element_line(color = "grey95", linetype = "dotted") 
#   )
# 
# # 4. 导出 PDF
# out_dir <- "Figures/Step07_C3Myeloid"
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# 
# ggsave(file.path(out_dir, "Fig07_GSEA_Lollipop_Compact.pdf"), 
#        plot = p_lollipop_compact, width = 8, height =5)
# 
# message("✅ 完美！紧凑排版的高级棒棒糖图已保存至: ", file.path(out_dir, "Fig07_GSEA_Lollipop_Compact.pdf"))
# ------------------------------------------------------------------------------
# 🌟 关键新增：计算 NES 的最大绝对值，并额外放大 30% 的空间
# 这能保证 0 刻度完美居中，且棒棒糖有足够长的距离伸展
# ------------------------------------------------------------------------------
max_nes <- max(abs(plot_df$NES), na.rm = TRUE) * 1.3 

# 3. 绘制紧凑型高级棒棒糖图
p_lollipop_compact <- ggplot(plot_df, aes(x = NES, y = Description)) +
  
  # 中心零线：作为上下调的绝对分水岭
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  
  # 棒棒糖的“棍子”
  geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), 
               color = "grey70", linewidth = 1.2) +
  
  # 棒棒糖的“糖”：气泡
  geom_point(aes(size = logFDR, fill = NES), 
             shape = 21, color = "black", stroke = 0.6) +
  
  # 顶级色彩映射
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, 
                       name = "NES") +
  
  # 调整气泡大小范围
  scale_size_continuous(range = c(3, 8), name = expression(-log[10](FDR))) +
  
  # 🌟 关键修复：强制设定 X 轴范围，拉长棒子的视觉比例
  scale_x_continuous(limits = c(-max_nes, max_nes), expand = expansion(mult = 0.02)) +
  
  # 精简主题设置
  theme_classic(base_size = 12) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL, 
    title = "Key GO Pathways in C3-high Myeloid"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text.y = element_text(color = "black", face = "bold", size = 10, lineheight = 0.8),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = ggplot2::margin(t = 10)),
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.margin = ggplot2::margin(0, 0, 0, -10), 
    panel.grid.major.y = element_line(color = "grey95", linetype = "dotted") 
  )

# 4. 导出 PDF
out_dir <- "Figures/Step07_C3Myeloid"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 🌟 关键修复：把宽度 (width) 从 8 增加到 10，给右侧的数据绘图区更多物理空间
# 高度 (height) 根据你选择了 Top 10，适当从 5 增加到 6.5 防止上下拥挤
ggsave(file.path(out_dir, "Fig07_GSEA_Lollipop_Compact.pdf"), 
       plot = p_lollipop_compact, width = 10, height = 6.5)

message("✅ 完美！延长棒子的高级棒棒糖图已保存至: ", file.path(out_dir, "Fig07_GSEA_Lollipop_Compact.pdf"))








# ─── 5. 代谢通路评分 ──────────────────────────────────────────────────────
message("🔋 代谢特征分析...")
# # 定义代谢基因集
# glycolysis_genes <- c("HK2", "PFKP", "ALDOA", "PKM", "LDHA")
# oxphos_genes <- c("COX5A", "COX6C", "ATP5F1A", "NDUFA4", "UQCRH")
# 
# sc_c3high <- AddModuleScore(sc_c3high, 
#                             features = list(glycolysis_genes, oxphos_genes),
#                             name = c("Glycolysis", "OXPHOS"))
# 
# p_metab <- FeaturePlot(sc_c3high, 
#                        features = c("Glycolysis1", "OXPHOS2"),
#                        ncol = 2) &
#   scale_color_gradientn(colors = c("grey90", "red"))
# 
# # save_pdf("Fig07D_Metabolism.pdf", p_metab, 
#           # width = 10, height = 4, subdir = "Step07_C3Myeloid")
#  # 定义 save_pdf 函数，跑一次就行！
# save_pdf <- function(filename, plot, width, height, subdir = "") {
#   # 假设你的基础输出目录是当前文件夹
#   full_dir <- subdir
#   
#   # 如果文件夹不存在则创建
#   if (!dir.exists(full_dir) && full_dir != "") {
#     dir.create(full_dir, recursive = TRUE)
#   }
#   
#   # 拼接完整的文件路径
#   full_path <- file.path(full_dir, filename)
#   
#   # 保存图片
#   ggplot2::ggsave(filename = full_path, plot = plot, width = width, height = height)
#   message("✅ 搞定！图片已保存至: ", full_path)
# }
# 
# # 现在你可以直接运行你原来的代码了：
# save_pdf("Fig07D_Metabolism.pdf", p_metab, 
#          width = 10, height = 4, subdir = "Step07_C3Myeloid")
# ==============================================================================
# ─── 补充：GSEA 进阶高级可视化 (山峰叠加图、基因网络图、富集网络图) ───────
# ==============================================================================

# # 确保加载专门用于富集分析可视化的核心包
# library(enrichplot)
# library(ggplot2)
# 
# message("🎨 正在生成 GSEA 进阶高级可视化图表...")
# 
# # 确保输出文件夹存在
# out_dir <- "Figures/Step07_C3Myeloid"
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# 
# 
# # ──────────────────────────────────────────────────────────────────────────────
# # 🎯 图 A：多通路 GSEA 叠加“山峰图” (对应你上传的图片)
# # ──────────────────────────────────────────────────────────────────────────────
# message("绘制多通路叠加山峰图...")
# 
# # 提取你最想展示的几条核心通路的 ID。
# # 方法 1：自动提取 P 值最小的 Top 5 通路
# top_pathway_ids <- gsea_go@result$ID[1:5]
# 
# # 方法 2：(推荐) 手动指定你关心的具有生物学意义的通路 ID
# # 例如你之前关心的补体或代谢通路，你可以去 gsea_go@result 里查阅对应的 GO:XXXXXXX
# # top_pathway_ids <- c("GO:0006956", "GO:0006958", "GO:0002250") 
# 
# # 使用 gseaplot2 绘制多通路叠加图
# p_gsea_multi <- gseaplot2(gsea_go, 
#                           geneSetID = top_pathway_ids, 
#                           title = "C3-high Myeloid vs C3-low: Key Pathways",
#                           color = scales::hue_pal()(length(top_pathway_ids)), # 自动生成好看的分类配色
#                           pvalue_table = FALSE, # 关掉 P 值表格，让图面更干净（像你上传的图一样）
#                           base_size = 12)
# 
# # 保存 (由于叠在一起较长，建议高度设置大一点)
# ggsave(file.path(out_dir, "Fig07_GSEA_Multi_Mountain.pdf"), plot = p_gsea_multi, width = 8, height = 6)




# # 1. 基因集定义 (保持不变)
# glycolysis_genes <- c("HK2", "PFKP", "ALDOA", "PKM", "LDHA")
# oxphos_genes <- c("COX5A", "COX6C", "ATP5F1A", "NDUFA4", "UQCRH")
# 
# # 2. 模块打分 (保持不变)
# sc_c3high <- AddModuleScore(sc_c3high, 
#                             features = list(glycolysis_genes, oxphos_genes),
#                             name = c("Glycolysis", "OXPHOS"))
# 
# # 3. 绘制 UMAP 图 (🌟 核心修改在这里)
# p_metab <- FeaturePlot(sc_c3high, 
#                        features = c("Glycolysis1", "OXPHOS2"),
#                        cols = c("grey90", "red"),  # 直接在函数内指定颜色，兼容性更好
#                        min.cutoff = "q5",          # 将最低 5% 的值强行压到底色 (灰色)
#                        max.cutoff = "q95",         # 将最高 5% 的极值强行压成最高色 (纯红)
#                        ncol = 2) 
# # 注意：把后面的 & scale_color_gradientn(...) 删掉，用 cols 参数替代
# 
# # # 4. 保存图片 (建议换个新名字对比一下)
# # save_pdf("Fig07D_Metabolism_Adjusted.pdf", p_metab, 
# #          width = 10, height = 4, subdir = "Step07_C3Myeloid")
# 
# # 定义 save_pdf 函数
# save_pdf <- function(filename, plot, width, height, subdir = "") {
#   # 拼接完整的目录路径
#   full_dir <- subdir
#   
#   # 如果文件夹不存在则自动创建
#   if (!dir.exists(full_dir) && full_dir != "") {
#     dir.create(full_dir, recursive = TRUE)
#   }
#   
#   # 拼接完整的文件路径
#   full_path <- file.path(full_dir, filename)
#   
#   # 使用 ggsave 保存图片
#   ggplot2::ggsave(filename = full_path, plot = plot, width = width, height = height)
#   message("✅ 搞定！图片已保存至: ", full_path)
# }
# 
# save_pdf("Fig07D_Metabolism_Adjusted.pdf", p_metab, 
#          width = 10, height = 4, subdir = "Step07_C3Myeloid")
# 
# 
# 








# ─── 6. 细胞通讯分析（需要完整数据集）────────────────────────────────────
message("📡 准备细胞通讯分析数据...")
# 标记C3-high myeloid在原始对象中
sc_new$C3high_myeloid <- colnames(sc_new) %in% colnames(sc_c3high)

# 保存用于CellChat的对象
saveRDS(sc_new, file.path(RDATA_DIR, "Step07_sc_for_cellchat.rds"))
message("💡 提示: 使用CellChat分析C3-high myeloid与其他细胞的互作")

# ─── 7. 保存结果 ──────────────────────────────────────────────────────────
saveRDS(sc_c3high, file.path(RDATA_DIR, "Step07_C3high_myeloid.rds"))
saveRDS(sc_myeloid, file.path(RDATA_DIR, "Step07_all_myeloid.rds"))

message("✅ 分析完成! 输出目录: ", OUT_DIR)


# ─── 5. 改进版：代谢通路评分（C3-high vs C3-low 严格对比） ────────────────
message("🔋 代谢特征分析 (全量 Myeloid 组间对比)...")

# 1. 定义代谢核心基因集
glycolysis_genes <- c("HK2", "PFKP", "ALDOA", "PKM", "LDHA")
oxphos_genes <- c("COX5A", "COX6C", "ATP5F1A", "NDUFA4", "UQCRH")

# 2. 🌟 关键改变：对全体髓系细胞 (sc_myeloid) 统一打分，而非仅对子集打分
sc_myeloid <- AddModuleScore(sc_myeloid, 
                             features = list(glycolysis_genes, oxphos_genes),
                             name = c("Glycolysis", "OXPHOS"))

# 整理列名，让后续代码和图例更清爽 (AddModuleScore 默认会加数字后缀)
sc_myeloid$Glycolysis_Score <- sc_myeloid$Glycolysis1
sc_myeloid$OXPHOS_Score <- sc_myeloid$OXPHOS2

# ──────────────────────────────────────────────────────────────────────────────
# 🎨 方式一：小提琴图 (Violin Plot) —— SCI 审稿人最爱的统计学证据图
# 目的：直接对比高低两组的分数差异，C3-high用红色，C3-low用蓝色
# ──────────────────────────────────────────────────────────────────────────────
library(ggplot2)

p_metab_vln <- VlnPlot(sc_myeloid, 
                       features = c("Glycolysis_Score", "OXPHOS_Score"),
                       group.by = "C3_status", 
                       cols = c("C3-high" = "#DE2D26", "C3-low" = "#3182BD"), # 🔴红蓝对比
                       pt.size = 0,    # 隐藏杂乱的单细胞散点，让图更干净高级
                       ncol = 2) & 
  geom_boxplot(width = 0.15, fill = "white", color = "black") & # 在小提琴内部加个白色的箱型图，展示中位数
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

save_pdf("Fig07D_Metabolism_VlnPlot_Compare.pdf", p_metab_vln, 
         width = 8, height = 5, subdir = "Step07_C3Myeloid")

# ──────────────────────────────────────────────────────────────────────────────
# 🎨 方式二：拆分 UMAP 涂色图 (Split FeaturePlot) —— 最强烈的视觉空间冲击
# 目的：左边画 C3-low，右边画 C3-high，看看高分细胞（红色）是不是全挤在右边
# ──────────────────────────────────────────────────────────────────────────────
p_metab_umap <- FeaturePlot(sc_myeloid, 
                            features = c("Glycolysis_Score", "OXPHOS_Score"),
                            split.by = "C3_status", # 🌟 核心：按照高低状态左右拆分画图
                            cols = c("grey90", "red"),
                            min.cutoff = "q5", 
                            max.cutoff = "q95",
                            keep.scale = "all") # 🌟 必须加这个：保证左右两张图的红颜色标准是完全一样的！

save_pdf("Fig07D_Metabolism_UMAP_Split.pdf", p_metab_umap, 
         width = 10, height = 8, subdir = "Step07_C3Myeloid")

message("✅ 代谢对比分析完成！两种美图已生成！")



# 
# 
# 
# # ==============================================================================
# # 08_CellChat_C3high_Myeloid.R
# # 功能: 基于Step07打标的C3-high myeloid细胞，进行细胞通讯分析
# #       重点分析 C3-high Myeloid ↔ 肿瘤细胞（myeSC / nSMC）的免疫通讯
# # 输入: Step07_sc_for_cellchat.rds
# # 输出: 全局通讯图 + C3high-Tumor专项分析图
# # ==============================================================================
# 
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(CellChat)
#   library(dplyr)
#   library(ggplot2)
#   library(patchwork)
#   library(ComplexHeatmap)
#   library(NMF)
#   library(ggalluvial)
# })
# 
# setwd("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/")
# set.seed(123)
# 
# OUT_DIR   <- "Figures/Step08_CellChat"
# RDATA_DIR <- "RData"
# if (!dir.exists(OUT_DIR))   dir.create(OUT_DIR,   recursive = TRUE)
# if (!dir.exists(RDATA_DIR)) dir.create(RDATA_DIR, recursive = TRUE)
# 
# # ─── 辅助函数 ─────────────────────────────────────────────────────────────────
# save_fig <- function(filename, plot, width = 8, height = 6) {
#   ggsave(file.path(OUT_DIR, filename), plot = plot,
#          width = width, height = height, dpi = 300)
#   message("✅ 已保存: ", filename)
# }
# 
# # ==============================================================================
# # ─── 1. 加载数据 & 构建细胞标签 ──────────────────────────────────────────────
# # ==============================================================================
# message("⏳ 加载数据...")
# sc <- readRDS(file.path(RDATA_DIR, "Step07_sc_for_cellchat.rds"))
# 
# # ── 1.1 构建新的细胞标签列：将 C3-high myeloid 单独标注 ──────────────────────
# # 原始 ident 中 Myeloid 细胞全部合并，这里细分为 C3-high 和 C3-low 两组
# sc$cellchat_label <- as.character(Idents(sc))
# 
# # 利用Step07打的标记（C3high_myeloid 为 TRUE/FALSE）
# if ("C3high_myeloid" %in% colnames(sc@meta.data)) {
#   sc$cellchat_label[sc$C3high_myeloid == TRUE]  <- "C3high_Myeloid"
#   sc$cellchat_label[sc$C3high_myeloid == FALSE &
#                       sc$cellchat_label == "Myeloid"] <- "C3low_Myeloid"
# } else {
#   # 备用方案：从 C3_status 恢复（如果 Step07 重新运行后有这一列）
#   if ("C3_status" %in% colnames(sc@meta.data)) {
#     sc$cellchat_label <- ifelse(
#       sc$cellchat_label == "Myeloid" & sc$C3_status == "C3-high",
#       "C3high_Myeloid",
#       ifelse(sc$cellchat_label == "Myeloid" & sc$C3_status == "C3-low",
#              "C3low_Myeloid", sc$cellchat_label)
#     )
#   } else {
#     stop("❌ 未找到 C3high_myeloid 或 C3_status 列，请先运行 Step07")
#   }
# }
# 
# # 确认肿瘤细胞类型名称（根据你的数据实际命名调整）
# tumor_types <- c("myeSC", "nSMC")    # ← 如果命名不同请修改
# all_types   <- unique(sc$cellchat_label)
# message("📋 当前所有细胞类型: ", paste(sort(all_types), collapse = " | "))
# message("🎯 肿瘤细胞类型: ",
#         paste(intersect(tumor_types, all_types), collapse = " | "))
# 
# # ── 1.2 定义颜色方案 ──────────────────────────────────────────────────────────
# color_palette <- c(
#   "C3high_Myeloid" = "#D73027",   # 深红 - 分析主角
#   "C3low_Myeloid"  = "#FC8D59",   # 橙红
#   "myeSC"          = "#4575B4",   # 深蓝 - 肿瘤
#   "nSMC"           = "#74ADD1",   # 蓝 - 肿瘤
#   "nmSC"           = "#91BFDB",   # 浅蓝
#   "Endothelial"    = "#E0F3F8",   # 最浅蓝（但用深色描边）
#   "Fibroblast"     = "#FEE090",   # 黄
#   "PC_VSMC"        = "#ABD9E9",
#   "Cycling"        = "#A6D96A",
#   "T_cell"         = "#1A9641",
#   "B_cell"         = "#FDAE61"
# )
# # 为未定义类型补充灰色
# missing_types <- setdiff(all_types, names(color_palette))
# if (length(missing_types) > 0) {
#   extra_colors <- setNames(
#     colorRampPalette(c("#CCCCCC", "#888888"))(length(missing_types)),
#     missing_types
#   )
#   color_palette <- c(color_palette, extra_colors)
# }
# 
# # ==============================================================================
# # ─── 2. 创建 CellChat 对象 ────────────────────────────────────────────────────
# # ==============================================================================
# message("🔨 创建 CellChat 对象...")
# 
# # 提取表达矩阵
# if (packageVersion("SeuratObject") >= "5.0.0") {
#   data_input <- LayerData(sc[["RNA"]], layer = "data")
# } else {
#   data_input <- GetAssayData(sc, assay = "RNA", slot = "data")
# }
# 
# # 细胞 meta
# meta_df <- data.frame(
#   labels   = sc$cellchat_label,
#   row.names = colnames(sc)
# )
# 
# cellchat <- createCellChat(
#   object   = data_input,
#   meta     = meta_df,
#   group.by = "labels"
# )
# 
# # ── 2.1 使用人类 CellChatDB ──────────────────────────────────────────────────
# CellChatDB <- CellChatDB.human
# # 使用全库（包含 Secreted Signaling、Cell-Cell Contact、ECM-Receptor）
# # 如只想用分泌型信号：CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# cellchat@DB <- CellChatDB
# 
# # ── 2.2 预处理 ────────────────────────────────────────────────────────────────
# message("⚙️ 预处理表达数据...")
# cellchat <- subsetData(cellchat)
# options(future.globals.maxSize = 10 * 1024^3)  # 设为 10 GiB，按实际内存调整
# future::plan("multisession", workers = 4)  # 并行加速，按服务器CPU核数调整
# 
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# 
# # ==============================================================================
# # ─── 3. 计算通讯概率 ──────────────────────────────────────────────────────────
# # ==============================================================================
# message("📡 计算细胞通讯概率...")
# cellchat <- computeCommunProb(cellchat, type = "triMean")
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)
# 
# # 保存中间结果
# saveRDS(cellchat, file.path(RDATA_DIR, "Step08_cellchat_object.rds"))
# message("✅ CellChat 对象已保存")
# 
# # ==============================================================================
# # ─── 4. 全局通讯概览 ──────────────────────────────────────────────────────────
# # ==============================================================================
# message("🎨 绘制全局通讯概览图...")
# 
# groupSize <- as.numeric(table(cellchat@idents))
# 
# # ── Fig 4A: 互作数量 + 强度汇总圆圈图 ────────────────────────────────────────
# pdf(file.path(OUT_DIR, "Fig08A_GlobalNet_Count_Strength.pdf"), width = 12, height = 5)
# par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
# netVisual_circle(cellchat@net$count,
#                  vertex.weight = groupSize,
#                  weight.scale  = TRUE,
#                  label.edge    = FALSE,
#                  color.use     = color_palette[levels(cellchat@idents)],
#                  title.name    = "Number of interactions")
# netVisual_circle(cellchat@net$weight,
#                  vertex.weight = groupSize,
#                  weight.scale  = TRUE,
#                  label.edge    = FALSE,
#                  color.use     = color_palette[levels(cellchat@idents)],
#                  title.name    = "Interaction strength")
# dev.off()
# message("✅ Fig08A 已保存")
# 
# # ── Fig 4B: Heatmap - 互作数量 ────────────────────────────────────────────────
# p_heatmap_count <- netVisual_heatmap(cellchat, measure = "count",
#                                      color.use = color_palette[levels(cellchat@idents)])
# p_heatmap_weight <- netVisual_heatmap(cellchat, measure = "weight",
#                                       color.use = color_palette[levels(cellchat@idents)])
# 
# pdf(file.path(OUT_DIR, "Fig08B_Heatmap_Count_Strength.pdf"), width = 12, height = 5)
# print(p_heatmap_count + p_heatmap_weight)
# dev.off()
# message("✅ Fig08B 已保存")
# 
# # ==============================================================================
# # ─── 5. C3-high Myeloid 专项通讯分析 ─────────────────────────────────────────
# # ==============================================================================
# message("🔬 C3-high Myeloid 专项分析...")
# 
# # ── 5.1 提取 C3high_Myeloid 的所有输出/输入通讯 ──────────────────────────────
# c3high_idx <- which(levels(cellchat@idents) == "C3high_Myeloid")
# 
# # C3-high 作为信号发送方
# pdf(file.path(OUT_DIR, "Fig08C_C3high_as_Sender.pdf"), width = 7, height = 6)
# par(mar = c(1, 1, 2, 1))
# netVisual_circle(cellchat@net$weight,
#                  sources.use   = c3high_idx,
#                  vertex.weight = groupSize,
#                  weight.scale  = TRUE,
#                  color.use     = color_palette[levels(cellchat@idents)],
#                  title.name    = "C3-high Myeloid as Sender")
# dev.off()
# 
# # C3-high 作为信号接收方
# pdf(file.path(OUT_DIR, "Fig08D_C3high_as_Receiver.pdf"), width = 7, height = 6)
# par(mar = c(1, 1, 2, 1))
# netVisual_circle(cellchat@net$weight,
#                  targets.use   = c3high_idx,
#                  vertex.weight = groupSize,
#                  weight.scale  = TRUE,
#                  color.use     = color_palette[levels(cellchat@idents)],
#                  title.name    = "C3-high Myeloid as Receiver")
# dev.off()
# message("✅ Fig08C/D 已保存")
# 
# # ── 5.2 C3-high Myeloid ↔ 肿瘤细胞（myeSC + nSMC）专项 ─────────────────────
# tumor_idx <- which(levels(cellchat@idents) %in% tumor_types)
# 
# if (length(tumor_idx) > 0) {
#   # 肿瘤相关细胞组合索引（包含 C3high）
#   focus_idx <- c(c3high_idx, tumor_idx)
#   
#   # 聚焦圆圈图：只看这三类细胞之间的通讯
#   pdf(file.path(OUT_DIR, "Fig08E_C3high_Tumor_Circle.pdf"), width = 7, height = 6)
#   par(mar = c(1, 1, 2, 1))
#   netVisual_circle(cellchat@net$weight,
#                    sources.use   = focus_idx,
#                    targets.use   = focus_idx,
#                    vertex.weight = groupSize,
#                    weight.scale  = TRUE,
#                    remove.isolate = FALSE,
#                    color.use     = color_palette[levels(cellchat@idents)],
#                    title.name    = "C3-high Myeloid ↔ Tumor Cells")
#   dev.off()
#   message("✅ Fig08E 已保存")
# } else {
#   warning("⚠️ 未找到肿瘤细胞类型 ", paste(tumor_types, collapse = "/"),
#           "，请检查 cellchat_label 中的实际名称")
# }
# 
# # ==============================================================================
# # ─── 6. 信号通路层面分析 ──────────────────────────────────────────────────────
# # ==============================================================================
# message("🛤️ 信号通路分析...")
# 
# # ── 6.1 通路汇总 Heatmap（发送 vs 接收）──────────────────────────────────────
# p_outgoing <- netAnalysis_signalingRole_heatmap(
#   cellchat, pattern = "outgoing",
#   color.use = color_palette[levels(cellchat@idents)],
#   width = 8, height = 8
# )
# p_incoming <- netAnalysis_signalingRole_heatmap(
#   cellchat, pattern = "incoming",
#   color.use = color_palette[levels(cellchat@idents)],
#   width = 8, height = 8
# )
# 
# pdf(file.path(OUT_DIR, "Fig08F_Pathway_Heatmap_Out_In.pdf"), width = 16, height = 8)
# print(p_outgoing + p_incoming)
# dev.off()
# message("✅ Fig08F 已保存")
# 
# # ── 6.2 气泡图：C3-high Myeloid → 肿瘤细胞的配体-受体对 ────────────────────
# if (length(tumor_idx) > 0) {
#   
#   # C3-high → 肿瘤
#   p_bubble_send <- netVisual_bubble(
#     cellchat,
#     sources.use = c3high_idx,
#     targets.use = tumor_idx,
#     remove.isolate = FALSE,
#     angle.x = 45
#   ) + ggtitle("C3-high Myeloid → Tumor (myeSC / nSMC)") +
#     theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#   
#   # 肿瘤 → C3-high
#   p_bubble_recv <- netVisual_bubble(
#     cellchat,
#     sources.use = tumor_idx,
#     targets.use = c3high_idx,
#     remove.isolate = FALSE,
#     angle.x = 45
#   ) + ggtitle("Tumor (myeSC / nSMC) → C3-high Myeloid") +
#     theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#   
#   save_fig("Fig08G_Bubble_C3high_to_Tumor.pdf",  p_bubble_send, width = 10, height = 8)
#   save_fig("Fig08H_Bubble_Tumor_to_C3high.pdf",  p_bubble_recv, width = 10, height = 8)
# }
# 
# # ── 6.3 重要通路的 Chord 图（弦图）──────────────────────────────────────────
# # 自动筛选 C3-high Myeloid 参与的前10条通路
# all_pathways <- cellchat@netP$pathways
# 
# # 找与 C3-high 高度相关的通路（出现在其 sender 或 receiver 信号中）
# pathway_contrib <- sapply(all_pathways, function(pw) {
#   mat <- cellchat@netP$prob[, , pw, drop = FALSE]
#   # C3-high 行 or 列的贡献
#   sum(mat[c3high_idx, , ]) + sum(mat[, c3high_idx, ])
# })
# top_pathways <- names(sort(pathway_contrib, decreasing = TRUE))[1:min(10, length(all_pathways))]
# message("🏆 C3-high Myeloid 参与的 Top 通路: ", paste(top_pathways, collapse = ", "))
# 
# # 为每个 Top 通路画弦图
# pdf(file.path(OUT_DIR, "Fig08I_Chord_TopPathways.pdf"), width = 8, height = 8)
# for (pw in top_pathways) {
#   tryCatch({
#     par(mar = c(1, 1, 3, 1))
#     netVisual_aggregate(
#       cellchat, signaling = pw,
#       layout       = "chord",
#       color.use    = color_palette[levels(cellchat@idents)],
#       title.name   = paste0(pw, " signaling pathway")
#     )
#   }, error = function(e) {
#     message("⚠️ 通路 ", pw, " 弦图绘制失败: ", e$message)
#   })
# }
# dev.off()
# message("✅ Fig08I 已保存")
# 
# # ==============================================================================
# # ─── 7. 深度聚焦：C3-high Myeloid & 肿瘤细胞专项通路 ─────────────────────────
# # ==============================================================================
# message("🔭 深度聚焦 C3-high ↔ 肿瘤细胞通路...")
# 
# # ── 7.1 提取 C3-high ↔ 肿瘤细胞间所有 L-R 对 ────────────────────────────────
# df_lr_c3_tumor <- subsetCommunication(
#   cellchat,
#   sources.use = c("C3high_Myeloid"),
#   targets.use = tumor_types
# )
# 
# df_lr_tumor_c3 <- subsetCommunication(
#   cellchat,
#   sources.use = tumor_types,
#   targets.use = c("C3high_Myeloid")
# )
# 
# df_lr_all <- bind_rows(
#   df_lr_c3_tumor %>% mutate(direction = "C3high→Tumor"),
#   df_lr_tumor_c3 %>% mutate(direction = "Tumor→C3high")
# )
# 
# write.csv(df_lr_all, file.path(OUT_DIR, "LR_C3high_Tumor_interactions.csv"),
#           row.names = FALSE)
# message("✅ L-R 互作列表已保存 (", nrow(df_lr_all), " 对)")
# 
# # ── 7.2 Top L-R 对气泡图（按 prob 强度排序，展示最强互作）────────────────────
# top_lr <- df_lr_all %>%
#   arrange(desc(prob)) %>%
#   slice_head(n = 30) %>%
#   mutate(pair = paste0(ligand, " → ", receptor))
# 
# p_top_lr <- ggplot(top_lr,
#                    aes(x = direction, y = reorder(pair, prob),
#                        size = prob, color = pathway_name)) +
#   geom_point(alpha = 0.85) +
#   scale_size_continuous(range = c(2, 8), name = "Interaction prob") +
#   scale_color_brewer(palette = "Set3", name = "Pathway") +
#   theme_classic(base_size = 12) +
#   labs(
#     title    = "Top L-R Interactions: C3-high Myeloid ↔ Tumor",
#     subtitle = "Size = interaction probability | Color = pathway",
#     x = NULL, y = NULL
#   ) +
#   theme(
#     plot.title    = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
#     axis.text.y   = element_text(size = 8),
#     legend.position = "right"
#   )
# save_fig("Fig08J_TopLR_C3high_Tumor.pdf", p_top_lr, width = 10, height = 10)
# 
# # ── 7.3 按通路汇总的 barplot（C3-high → 肿瘤 vs 肿瘤 → C3-high）─────────────
# pathway_summary <- df_lr_all %>%
#   group_by(pathway_name, direction) %>%
#   summarise(total_prob = sum(prob, na.rm = TRUE),
#             n_pairs    = n(), .groups = "drop")
# 
# p_pathway_bar <- ggplot(pathway_summary,
#                         aes(x = reorder(pathway_name, total_prob),
#                             y = total_prob, fill = direction)) +
#   geom_col(position = "dodge", width = 0.7) +
#   scale_fill_manual(values = c("C3high→Tumor" = "#D73027",
#                                "Tumor→C3high" = "#4575B4")) +
#   coord_flip() +
#   theme_classic(base_size = 12) +
#   labs(
#     title = "Signaling Pathways: C3-high Myeloid ↔ Tumor",
#     x = "Pathway", y = "Total interaction probability",
#     fill = "Direction"
#   ) +
#   theme(plot.title = element_text(face = "bold", hjust = 0.5))
# 
# save_fig("Fig08K_Pathway_Bar_C3high_Tumor.pdf", p_pathway_bar, width = 9, height = 6)
# 
# # ── 7.4 重要免疫抑制通路的细节可视化 ─────────────────────────────────────────
# # 自动检测是否存在经典免疫检查点/免疫抑制通路
# key_immune_pathways <- c(
#   "MHC-II", "MHC-I",      # 抗原呈递
#   "CD86",   "CD80",       # 共刺激
#   "PD-L1",  "PD-L2",      # 免疫检查点
#   "TGFb",   "TGFB",       # 免疫抑制
#   "CXCL",   "CCL",        # 趋化因子
#   "CSF",    "IL1",        # 细胞因子
#   "COMPLEMENT",            # 补体（与C3相关！）
#   "SPP1",   "VEGF",       # TME重塑
#   "ANGPTL", "FN1"
# )
# available_immune_pw <- intersect(key_immune_pathways, all_pathways)
# message("🎯 检测到的关键免疫通路: ", paste(available_immune_pw, collapse = ", "))
# 
# if (length(available_immune_pw) > 0) {
#   pdf(file.path(OUT_DIR, "Fig08L_KeyImmune_Pathways_Chord.pdf"), width = 8, height = 8)
#   for (pw in available_immune_pw) {
#     tryCatch({
#       par(mar = c(1, 1, 3, 1))
#       netVisual_aggregate(
#         cellchat, signaling = pw,
#         layout     = "chord",
#         color.use  = color_palette[levels(cellchat@idents)],
#         title.name = paste0(pw, " (key immune pathway)")
#       )
#     }, error = function(e) {
#       message("⚠️ ", pw, " 绘制失败: ", e$message)
#     })
#   }
#   dev.off()
#   message("✅ Fig08L 已保存")
# }
# 
# # ── 7.5 补体通路专项（与 C3 高度相关）────────────────────────────────────────
# complement_pws <- grep("COMPLEMENT|C1Q|C3|C5|COLLAGEN",
#                        all_pathways, value = TRUE, ignore.case = TRUE)
# if (length(complement_pws) > 0) {
#   message("💎 补体相关通路: ", paste(complement_pws, collapse = ", "))
#   
#   pdf(file.path(OUT_DIR, "Fig08M_Complement_Pathway.pdf"), width = 8, height = 8)
#   for (pw in complement_pws) {
#     tryCatch({
#       par(mar = c(1, 1, 3, 1))
#       netVisual_aggregate(
#         cellchat, signaling = pw,
#         layout     = "chord",
#         color.use  = color_palette[levels(cellchat@idents)],
#         title.name = paste0(pw, " signaling (complement-related)")
#       )
#     }, error = function(e) message("⚠️ ", pw, ": ", e$message))
#   }
#   dev.off()
#   
#   # 补体通路配体-受体明细
#   df_complement <- subsetCommunication(cellchat, signaling = complement_pws)
#   write.csv(df_complement, file.path(OUT_DIR, "Complement_LR_details.csv"),
#             row.names = FALSE)
#   message("✅ Fig08M + 补体通路 L-R 明细已保存")
# }
# 
# # ==============================================================================
# # ─── 8. 信号通路模式分析（NMF 分解）─────────────────────────────────────────
# # ==============================================================================
# message("🧩 通讯模式 NMF 分析...")
# 
# # 确认最优 pattern 数（通常 3~5）
# # selectK(cellchat, pattern = "outgoing")  # 可先解注释看 CophenCor 图
# 
# cellchat <- identifyCommunicationPatterns(
#   cellchat, pattern = "outgoing", k = 4,    # k 根据 selectK 结果调整
#   width = 8, height = 6,
#   color.use = color_palette[levels(cellchat@idents)]
# )
# cellchat <- identifyCommunicationPatterns(
#   cellchat, pattern = "incoming", k = 4,
#   width = 8, height = 6,
#   color.use = color_palette[levels(cellchat@idents)]
# )
# 
# # River plot（展示细胞-通路-模式关系）
# p_river_out <- netAnalysis_river(cellchat, pattern = "outgoing",
#                                  color.use = color_palette[levels(cellchat@idents)])
# p_river_in  <- netAnalysis_river(cellchat, pattern = "incoming",
#                                  color.use = color_palette[levels(cellchat@idents)])
# 
# save_fig("Fig08N_Pattern_River_Outgoing.pdf", p_river_out, width = 10, height = 6)
# save_fig("Fig08O_Pattern_River_Incoming.pdf", p_river_in,  width = 10, height = 6)
# 
# # ==============================================================================
# # ─── 9. 2D 信号功能分析（Scatter）──────────────────────────────────────────
# # ==============================================================================
# message("📊 信号功能 2D 分析...")
# 
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# 
# # Scatter: 出入强度（找到 C3-high Myeloid 的位置）
# p_scatter <- netAnalysis_signalingRole_scatter(
#   cellchat,
#   color.use = color_palette[levels(cellchat@idents)]
# ) +
#   geom_point(data = . %>% filter(labels == "C3high_Myeloid"),
#              aes(x = x, y = y), color = "red", size = 5, shape = 18) +
#   ggtitle("Signaling Role: Sender vs Receiver\n(Diamond = C3-high Myeloid)")
# 
# save_fig("Fig08P_SignalingRole_Scatter.pdf", p_scatter, width = 7, height = 6)
# message("✅ Fig08P 已保存")
# 
# # ==============================================================================
# # ─── 10. 保存最终结果 ─────────────────────────────────────────────────────────
# # ==============================================================================
# saveRDS(cellchat, file.path(RDATA_DIR, "Step08_cellchat_final.rds"))
# 
# # 汇总通讯强度表格
# df_all_comm <- subsetCommunication(cellchat)
# write.csv(df_all_comm, file.path(OUT_DIR, "All_LR_interactions.csv"), row.names = FALSE)
# 
# message("\n", paste(rep("=", 60), collapse = ""))
# message("✅ 所有分析完成！输出目录: ", OUT_DIR)
# message("📁 主要输出文件:")
# message("   Fig08A  全局互作圆圈图（数量 + 强度）")
# message("   Fig08B  全局互作 Heatmap")
# message("   Fig08C  C3-high 作为信号发送方")
# message("   Fig08D  C3-high 作为信号接收方")
# message("   Fig08E  C3-high ↔ 肿瘤细胞聚焦圆圈图")
# message("   Fig08F  通路 Heatmap（发送 vs 接收）")
# message("   Fig08G/H C3-high ↔ 肿瘤 L-R 气泡图")
# message("   Fig08I  Top 通路弦图")
# message("   Fig08J  Top L-R 互作点图")
# message("   Fig08K  通路强度 barplot")
# message("   Fig08L  关键免疫通路弦图")
# message("   Fig08M  补体通路专项分析")
# message("   Fig08N/O 通讯模式 River plot")
# message("   Fig08P  信号功能 2D Scatter")
# message("   LR_C3high_Tumor_interactions.csv")
# message(paste(rep("=", 60), collapse = ""))



# ==============================================================================
# 08_CellChat_C3high_Myeloid.R (🚀 提速优化版)
# 功能: 基于Step07打标的C3-high myeloid细胞，进行细胞通讯分析
#       重点分析 C3-high Myeloid ↔ 肿瘤细胞（myeSC / nSMC）的免疫通讯
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(NMF)
  library(ggalluvial)
})

setwd("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/")
set.seed(123)

OUT_DIR   <- "Figures/Step08_CellChat"
RDATA_DIR <- "RData"
if (!dir.exists(OUT_DIR))   dir.create(OUT_DIR,   recursive = TRUE)
if (!dir.exists(RDATA_DIR)) dir.create(RDATA_DIR, recursive = TRUE)

# ─── 辅助函数 ─────────────────────────────────────────────────────────────────
save_fig <- function(filename, plot, width = 8, height = 6) {
  ggsave(file.path(OUT_DIR, filename), plot = plot,
         width = width, height = height, dpi = 300)
  message("✅ 已保存: ", filename)
}

# ==============================================================================
# ─── 1. 加载数据 & 构建细胞标签 ──────────────────────────────────────────────
# ==============================================================================
message("⏳ 加载数据...")
sc <- readRDS(file.path(RDATA_DIR, "Step07_sc_for_cellchat.rds"))

sc$cellchat_label <- as.character(Idents(sc))

if ("C3high_myeloid" %in% colnames(sc@meta.data)) {
  sc$cellchat_label[sc$C3high_myeloid == TRUE]  <- "C3high_Myeloid"
  sc$cellchat_label[sc$C3high_myeloid == FALSE &
                      sc$cellchat_label == "Myeloid"] <- "C3low_Myeloid"
} else {
  if ("C3_status" %in% colnames(sc@meta.data)) {
    sc$cellchat_label <- ifelse(
      sc$cellchat_label == "Myeloid" & sc$C3_status == "C3-high",
      "C3high_Myeloid",
      ifelse(sc$cellchat_label == "Myeloid" & sc$C3_status == "C3-low",
             "C3low_Myeloid", sc$cellchat_label)
    )
  } else {
    stop("❌ 未找到 C3high_myeloid 或 C3_status 列，请先运行 Step07")
  }
}

tumor_types <- c("myeSC", "nSMC")
all_types   <- unique(sc$cellchat_label)
message("📋 当前所有细胞类型: ", paste(sort(all_types), collapse = " | "))
message("🎯 肿瘤细胞类型: ",
        paste(intersect(tumor_types, all_types), collapse = " | "))

# color_palette <- c(
#   "C3high_Myeloid" = "#D73027",   
#   "C3low_Myeloid"  = "#FC8D59",   
#   "myeSC"          = "#4575B4",   
#   "nSMC"           = "#74ADD1",   
#   "nmSC"           = "#91BFDB",   
#   "Endothelial"    = "#E0F3F8",   
#   "Fibroblast"     = "#FEE090",   
#   "PC_VSMC"        = "#ABD9E9",
#   "Cycling"        = "#A6D96A",
#   "T_cell"         = "#1A9641",
#   "B_cell"         = "#FDAE61"
# )

color_palette <- c(
  "C3high_Myeloid" = "#D73027",   
  "C3low_Myeloid"  = "#FC8D59",   
  "myeSC"       = "#8ECFC9",
  "nmSC"        = "#E5A84B",
  "Fibroblast"  = "#87CEEB",
  "PC_VSMC"     = "#82B0D2",
  "Endothelial" = "#4169E1",
  # "Myeloid"     = "#F7BDBE",
  "Mast"        = "#F6CAE5",
  "TC"          = "#96C37D",
  "NKC"         = "#F3D266",
  "BC"          = "#DDA0DD",
  "Cycling"     = "#A4312A",
  "Mucosa"      = "#FF7F50"
)





missing_types <- setdiff(all_types, names(color_palette))
if (length(missing_types) > 0) {
  extra_colors <- setNames(
    colorRampPalette(c("#CCCCCC", "#888888"))(length(missing_types)),
    missing_types
  )
  color_palette <- c(color_palette, extra_colors)
}

# ==============================================================================
# ─── 2. 创建 CellChat 对象 ────────────────────────────────────────────────────
# ==============================================================================
message("🔨 创建 CellChat 对象...")

if (packageVersion("SeuratObject") >= "5.0.0") {
  data_input <- LayerData(sc[["RNA"]], layer = "data")
} else {
  data_input <- GetAssayData(sc, assay = "RNA", slot = "data")
}

meta_df <- data.frame(
  labels   = sc$cellchat_label,
  row.names = colnames(sc)
)

cellchat <- createCellChat(
  object   = data_input,
  meta     = meta_df,
  group.by = "labels"
)

# 使用人类全库，保证信号通路的完整性（不删减，靠算法提速）
cellchat@DB <- CellChatDB.human

# ── 2.2 预处理 (🚀 加速优化点 1：调整 Linux 并行策略) ──────────────────────
message("⚙️ 预处理表达数据...")
cellchat <- subsetData(cellchat)

# 【修改说明】: 内存配额提升至 20G，使用 Linux 更高效的 multicore 替换 multisession，进程数设为 8
options(future.globals.maxSize = 20 * 1024^3) 
future::plan("multicore", workers = 8)  # 如果你的服务器有几十个核心，这里可以改大比如 16 或 24

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ==============================================================================
# ─── 3. 计算通讯概率 (🚀 加速优化点 2：更改核心计算算法) ────────────────────
# ==============================================================================
message("📡 计算细胞通讯概率...")

# 【修改说明】: type 从 triMean 改为 官方推荐的针对大数据集的 truncatedMean (截断前10%极值)，大幅提速且可靠
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path(RDATA_DIR, "Step08_cellchat_object.rds"))
message("✅ CellChat 对象已保存")

# ==============================================================================
# ─── 4. 全局通讯概览 ──────────────────────────────────────────────────────────
# ==============================================================================
message("🎨 绘制全局通讯概览图...")

groupSize <- as.numeric(table(cellchat@idents))

pdf(file.path(OUT_DIR, "Fig08A_GlobalNet_Count_Strength.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, color.use = color_palette[levels(cellchat@idents)],
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, color.use = color_palette[levels(cellchat@idents)],
                 title.name = "Interaction strength")
dev.off()
message("✅ Fig08A 已保存")

p_heatmap_count <- netVisual_heatmap(cellchat, measure = "count", color.use = color_palette[levels(cellchat@idents)])
p_heatmap_weight <- netVisual_heatmap(cellchat, measure = "weight", color.use = color_palette[levels(cellchat@idents)])

pdf(file.path(OUT_DIR, "Fig08B_Heatmap_Count_Strength.pdf"), width = 12, height = 5)
print(p_heatmap_count + p_heatmap_weight)
dev.off()
message("✅ Fig08B 已保存")

# ==============================================================================
# ─── 5. C3-high Myeloid 专项通讯分析 ─────────────────────────────────────────
# ==============================================================================
message("🔬 C3-high Myeloid 专项分析...")

c3high_idx <- which(levels(cellchat@idents) == "C3high_Myeloid")

pdf(file.path(OUT_DIR, "Fig08C_C3high_as_Sender.pdf"), width = 7, height = 6)
par(mar = c(1, 1, 2, 1))
netVisual_circle(cellchat@net$weight, sources.use = c3high_idx, vertex.weight = groupSize,
                 weight.scale = TRUE, color.use = color_palette[levels(cellchat@idents)],
                 title.name = "C3-high Myeloid as Sender")
dev.off()

pdf(file.path(OUT_DIR, "Fig08D_C3high_as_Receiver.pdf"), width = 7, height = 6)
par(mar = c(1, 1, 2, 1))
netVisual_circle(cellchat@net$weight, targets.use = c3high_idx, vertex.weight = groupSize,
                 weight.scale = TRUE, color.use = color_palette[levels(cellchat@idents)],
                 title.name = "C3-high Myeloid as Receiver")
dev.off()

tumor_idx <- which(levels(cellchat@idents) %in% tumor_types)

if (length(tumor_idx) > 0) {
  focus_idx <- c(c3high_idx, tumor_idx)
  
  pdf(file.path(OUT_DIR, "Fig08E_C3high_Tumor_Circle.pdf"), width = 7, height = 6)
  par(mar = c(1, 1, 2, 1))
  netVisual_circle(cellchat@net$weight, sources.use = focus_idx, targets.use = focus_idx,
                   vertex.weight = groupSize, weight.scale = TRUE, remove.isolate = FALSE,
                   color.use = color_palette[levels(cellchat@idents)],
                   title.name = "C3-high Myeloid ↔ Tumor Cells")
  dev.off()
} else {
  warning("⚠️ 未找到肿瘤细胞类型，请检查 cellchat_label 中的实际名称")
}

# ==============================================================================
# ─── 6. 信号通路层面分析 ──────────────────────────────────────────────────────
# ==============================================================================
message("🛤️ 信号通路分析...")
# 补算网络中心度打分
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
p_outgoing <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = color_palette[levels(cellchat@idents)], width = 16, height =25)
p_incoming <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = color_palette[levels(cellchat@idents)], width = 16, height = 25)

pdf(file.path(OUT_DIR, "Fig08F_Pathway_Heatmap_Out_In.pdf"), width = 16, height = 40)
print(p_outgoing + p_incoming)
dev.off()

if (length(tumor_idx) > 0) {
  p_bubble_send <- netVisual_bubble(cellchat, sources.use = c3high_idx, targets.use = tumor_idx, remove.isolate = FALSE, angle.x = 45) + 
    ggtitle("C3-high Myeloid → Tumor (myeSC / nSMC)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p_bubble_recv <- netVisual_bubble(cellchat, sources.use = tumor_idx, targets.use = c3high_idx, remove.isolate = FALSE, angle.x = 45) + 
    ggtitle("Tumor (myeSC / nSMC) → C3-high Myeloid") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  save_fig("Fig08G_Bubble_C3high_to_Tumor.pdf",  p_bubble_send, width = 10, height = 8)
  save_fig("Fig08H_Bubble_Tumor_to_C3high.pdf",  p_bubble_recv, width = 10, height = 8)
}

all_pathways <- cellchat@netP$pathways

pathway_contrib <- sapply(all_pathways, function(pw) {
  mat <- cellchat@netP$prob[, , pw, drop = FALSE]
  sum(mat[c3high_idx, , ]) + sum(mat[, c3high_idx, ])
})
top_pathways <- names(sort(pathway_contrib, decreasing = TRUE))[1:min(10, length(all_pathways))]

pdf(file.path(OUT_DIR, "Fig08I_Chord_TopPathways.pdf"), width = 8, height = 8)
for (pw in top_pathways) {
  tryCatch({
    par(mar = c(1, 1, 3, 1))
    netVisual_aggregate(cellchat, signaling = pw, layout = "chord", color.use = color_palette[levels(cellchat@idents)], title.name = paste0(pw, " signaling pathway"))
  }, error = function(e) {
    message("⚠️ 通路 ", pw, " 弦图绘制失败: ", e$message)
  })
}
dev.off()

# ==============================================================================
# ─── 7. 深度聚焦：C3-high Myeloid & 肿瘤细胞专项通路 ─────────────────────────
# ==============================================================================
message("🔭 深度聚焦 C3-high ↔ 肿瘤细胞通路...")

df_lr_c3_tumor <- subsetCommunication(cellchat, sources.use = c("C3high_Myeloid"), targets.use = tumor_types)
df_lr_tumor_c3 <- subsetCommunication(cellchat, sources.use = tumor_types, targets.use = c("C3high_Myeloid"))

df_lr_all <- bind_rows(
  df_lr_c3_tumor %>% mutate(direction = "C3high→Tumor"),
  df_lr_tumor_c3 %>% mutate(direction = "Tumor→C3high")
)

write.csv(df_lr_all, file.path(OUT_DIR, "LR_C3high_Tumor_interactions.csv"), row.names = FALSE)

top_lr <- df_lr_all %>%
  arrange(desc(prob)) %>%
  slice_head(n = 50) %>%
  mutate(pair = paste0(ligand, " → ", receptor))

p_top_lr <- ggplot(top_lr, aes(x = direction, y = reorder(pair, prob), size = prob, color = pathway_name)) +
  geom_point(alpha = 0.85) + scale_size_continuous(range = c(2, 8), name = "Interaction prob") +
  scale_color_brewer(palette = "Set3", name = "Pathway") + theme_classic(base_size = 12) +
  labs(title = "Top L-R Interactions: C3-high Myeloid ↔ Tumor", subtitle = "Size = interaction probability | Color = pathway", x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "grey40"), axis.text.y = element_text(size = 8))
save_fig("Fig08J_TopLR_C3high_Tumor.pdf", p_top_lr, width = 10, height = 10)

pathway_summary <- df_lr_all %>%
  group_by(pathway_name, direction) %>%
  summarise(total_prob = sum(prob, na.rm = TRUE), n_pairs = n(), .groups = "drop")

p_pathway_bar <- ggplot(pathway_summary, aes(x = reorder(pathway_name, total_prob), y = total_prob, fill = direction)) +
  geom_col(position = "dodge", width = 0.7) + scale_fill_manual(values = c("C3high→Tumor" = "#D73027", "Tumor→C3high" = "#4575B4")) +
  coord_flip() + theme_classic(base_size = 12) + labs(title = "Signaling Pathways: C3-high Myeloid ↔ Tumor", x = "Pathway", y = "Total interaction probability", fill = "Direction") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_fig("Fig08K_Pathway_Bar_C3high_Tumor.pdf", p_pathway_bar, width = 9, height = 6)

key_immune_pathways <- c("MHC-II", "MHC-I", "CD86", "CD80", "PD-L1", "PD-L2", "TGFb", "TGFB", "CXCL", "CCL", "CSF", "IL1", "COMPLEMENT", "SPP1", "VEGF", "ANGPTL", "FN1")
available_immune_pw <- intersect(key_immune_pathways, all_pathways)

if (length(available_immune_pw) > 0) {
  pdf(file.path(OUT_DIR, "Fig08L_KeyImmune_Pathways_Chord.pdf"), width = 8, height = 8)
  for (pw in available_immune_pw) {
    tryCatch({
      par(mar = c(1, 1, 3, 1))
      netVisual_aggregate(cellchat, signaling = pw, layout = "chord", color.use = color_palette[levels(cellchat@idents)], title.name = paste0(pw, " (key immune pathway)"))
    }, error = function(e) { message("⚠️ ", pw, " 绘制失败: ", e$message) })
  }
  dev.off()
}

complement_pws <- grep("COMPLEMENT|C1Q|C3|C5|COLLAGEN", all_pathways, value = TRUE, ignore.case = TRUE)
if (length(complement_pws) > 0) {
  pdf(file.path(OUT_DIR, "Fig08M_Complement_Pathway.pdf"), width = 8, height = 8)
  for (pw in complement_pws) {
    tryCatch({
      par(mar = c(1, 1, 3, 1))
      netVisual_aggregate(cellchat, signaling = pw, layout = "chord", color.use = color_palette[levels(cellchat@idents)], title.name = paste0(pw, " signaling (complement-related)"))
    }, error = function(e) message("⚠️ ", pw, ": ", e$message))
  }
  dev.off()
  
  df_complement <- subsetCommunication(cellchat, signaling = complement_pws)
  write.csv(df_complement, file.path(OUT_DIR, "Complement_LR_details.csv"), row.names = FALSE)
}

# ==============================================================================
# ─── 8. 信号通路模式分析（NMF 分解）─────────────────────────────────────────
# ==============================================================================
message("🧩 通讯模式 NMF 分析...")

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 4, width = 8, height = 6, color.use = color_palette[levels(cellchat@idents)])
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 4, width = 8, height = 6, color.use = color_palette[levels(cellchat@idents)])

p_river_out <- netAnalysis_river(cellchat, pattern = "outgoing", color.use = color_palette[levels(cellchat@idents)])
p_river_in  <- netAnalysis_river(cellchat, pattern = "incoming", color.use = color_palette[levels(cellchat@idents)])

save_fig("Fig08N_Pattern_River_Outgoing.pdf", p_river_out, width = 10, height = 6)
save_fig("Fig08O_Pattern_River_Incoming.pdf", p_river_in,  width = 10, height = 6)

# ==============================================================================
# ─── 9. 2D 信号功能分析（Scatter）──────────────────────────────────────────
# ==============================================================================
message("📊 信号功能 2D 分析...")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

p_scatter <- netAnalysis_signalingRole_scatter(cellchat, color.use = color_palette[levels(cellchat@idents)]) +
  geom_point(data = . %>% filter(labels == "C3high_Myeloid"), aes(x = x, y = y), color = "red", size = 5, shape = 18) +
  ggtitle("Signaling Role: Sender vs Receiver\n(Diamond = C3-high Myeloid)")

save_fig("Fig08P_SignalingRole_Scatter.pdf", p_scatter, width = 7, height = 6)




# ==============================================================================
# ─── 10. 保存最终结果 ─────────────────────────────────────────────────────────
# ==============================================================================
saveRDS(cellchat, file.path(RDATA_DIR, "Step08_cellchat_final.rds"))

df_all_comm <- subsetCommunication(cellchat)
write.csv(df_all_comm, file.path(OUT_DIR, "All_LR_interactions.csv"), row.names = FALSE)

message("\n", paste(rep("=", 60), collapse = ""))
message("✅ 所有分析完成！")











# ==============================================================================
# 08B_CellChat_Replot_HQ.R (高清重绘王牌图专用)
# 功能: 直接读取已计算好的 CellChat RDS，应用顶级 SCI 配色，重新生成核心主图
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
})

# 设置工作目录和输出路径
setwd("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/")
OUT_DIR <- "Figures/Step08_CellChat/Step08_CellChat_HQ"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ─── 1. 加载已经算好的 CellChat 对象 ──────────────────────────────────────────
message("⏳ 正在加载 CellChat 结果数据 (这可能需要几秒钟)...")
cellchat <- readRDS("RData/Step08_cellchat_final.rds")

# ─── 2. 定义顶级 SCI 配色方案 (Nature Publishing Group 风格) ────────────────
# 核心原则：主角用醒目暖色(红/橙)，肿瘤用深浅冷色(蓝)，基质细胞用大地色系，其他辅助色
# color_palette <- c(
#   "C3high_Myeloid" = "#E64B35",  # 经典高亮红 (主角，极度醒目)
#   "C3low_Myeloid"  = "#F39B7F",  # 浅砖红/珊瑚橘
#   "myeSC"          = "#4DBBD5",  # 亮孔雀蓝 (主要肿瘤)
#   "nSMC"           = "#3C5488",  # 沉稳深藏青 (次要肿瘤)
#   "nmSC"           = "#8491B4",  # 高级灰蓝
#   "Endothelial"    = "#00A087",  # 祖母绿 (血管)
#   "Fibroblast"     = "#7E6148",  # 深褐色/咖啡色 (成纤维基质)
#   "PC_VSMC"        = "#B09C85",  # 浅卡其色
#   "T_cell"         = "#4DAF4A",  # 鲜艳的绿色 (如果有)
#   "B_cell"         = "#984EA3",  # 紫色 (如果有)
#   "Cycling"        = "#CCCCCC"   # 灰色 (增殖细胞)
# )
color_palette <- c(
  "C3high_Myeloid" = "#F09997",   
  "C3low_Myeloid"  = "#FC8D59",   
  "myeSC"       = "#8ECFC9",
  "nmSC"        = "#E5A84B",
  "Fibroblast"  = "#87CEEB",
  "PC_VSMC"     = "#82B0D2",
  "Endothelial" = "#4169E1",
  # "Myeloid"     = "#F7BDBE",
  "Mast"        = "#F6CAE5",
  "TC"          = "#96C37D",
  "NKC"         = "#F3D266",
  "BC"          = "#DDA0DD",
  "Cycling"     = "#A4312A",
  "Mucosa"      = "#FF7F50"
)


# 智能补全可能遗漏的细胞颜色（用中性灰兜底）
missing_types <- setdiff(levels(cellchat@idents), names(color_palette))
if (length(missing_types) > 0) {
  extra_colors <- setNames(rep("#D3D3D3", length(missing_types)), missing_types)
  color_palette <- c(color_palette, extra_colors)
}
# 确保颜色顺序与数据 levels 一致
color_use <- color_palette[levels(cellchat@idents)]

# 提前准备好细胞索引
c3high_idx <- which(levels(cellchat@idents) == "C3high_Myeloid")
tumor_types <- c("myeSC", "nSMC")
tumor_idx <- which(levels(cellchat@idents) %in% tumor_types)
focus_idx <- c(c3high_idx, tumor_idx)


# ─── 3. 重新输出 Panel A: 2D 信号地位散点图 ────────────────────────────────
message("🎨 正在重绘 Panel A: 2D 散点图...")
p_scatter <- netAnalysis_signalingRole_scatter(cellchat, color.use = color_use) +
  geom_point(data = . %>% filter(labels == "C3high_Myeloid"), 
             aes(x = x, y = y), color = "#E64B35", size = 6, shape = 18) + # 超大红色菱形高亮
  ggtitle("Signaling Role: Sender vs Receiver") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(OUT_DIR, "Fig_PanelA_Scatter.pdf"), p_scatter, width = 7, height = 6)


# ─── 4. 重新输出 Panel B: 无自环的放射状和弦图 ──────────────────────────────
message("🎨 正在重绘 Panel B: 放射状双向奔赴图...")
if (length(tumor_idx) > 0) {
  cellchat_temp <- cellchat
  diag(cellchat_temp@net$weight) <- 0 # 挖去自环
  diag(cellchat_temp@net$count) <- 0
  
  pdf(file.path(OUT_DIR, "Fig_PanelB_Radial_Chord.pdf"), width = 7, height = 7)
  par(mar = c(1, 1, 1, 1))
  netVisual_chord_cell(cellchat_temp,
                       net = cellchat_temp@net$weight,
                       sources.use = focus_idx,
                       targets.use = focus_idx,
                       color.use = color_use,
                       title.name = "Interactions: C3-high Myeloid <-> Tumor")
  dev.off()
}


# ─── 5. 重新输出 Panel C: L-R 核心机制气泡图 ────────────────────────────────
message("🎨 正在重绘 Panel C: Top 配体受体气泡图...")
if (length(tumor_idx) > 0) {
  # C3-high 发给 Tumor 的信号
  p_bubble_send <- netVisual_bubble(cellchat, sources.use = c3high_idx, targets.use = tumor_idx, 
                                    remove.isolate = FALSE, angle.x = 45) + 
    ggtitle("C3-high Myeloid -> Tumor") + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "#E64B35"))
  
  # Tumor 发给 C3-high 的信号
  p_bubble_recv <- netVisual_bubble(cellchat, sources.use = tumor_idx, targets.use = c3high_idx, 
                                    remove.isolate = FALSE, angle.x = 45) + 
    ggtitle("Tumor -> C3-high Myeloid") + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "#3C5488"))
  
  # 拼接到一张长图里
  library(patchwork)
  p_bubble_combine <- p_bubble_send + p_bubble_recv
  ggsave(file.path(OUT_DIR, "Fig_PanelC_Bubble_Combined.pdf"), p_bubble_combine, width = 14, height = 7)
}


# # ─── 6. 重新输出 Panel D-F: SPP1, MIF, COMPLEMENT 三大核心通路特写 ────────
# message("🎨 正在重绘 三大核心机制通路特写...")
# key_pathways <- c("SPP1", "MIF", "COMPLEMENT","APP")
# 
# for (pw in key_pathways) {
#   # 检查该通路是否存在于当前分析中
#   if (pw %in% cellchat@netP$pathways) {
#     file_name <- paste0("Fig_PanelD_", pw, "_Chord.pdf")
#     pdf(file.path(OUT_DIR, file_name), width = 7, height = 7)
#     par(mar = c(1, 1, 2, 1))
#     
#     # 绘制特定通路的弦图，带上新颜色
#     netVisual_aggregate(cellchat, signaling = pw, layout = "chord", 
#                         color.use = color_use, 
#                         title.name = paste0(pw, " signaling pathway"))
#     dev.off()
#   } else {
#     message("⚠️ 提示: 当前对象中未检测到 ", pw, " 通路，已跳过。")
#   }
# }
# 
# message("\n✅ 搞定！所有超高颜值的主图都已经存到 `Figures/Step08_Ce
#         
#         llChat_HQ` 文件夹下啦！")
# ─── 6. 重新输出 Panel D-F: SPP1, MIF, COMPLEMENT, APP 核心通路特写 ────────
message("🎨 正在重绘 三大核心机制通路特写 (固定细胞位置版)...")
key_pathways <- c("SPP1", "MIF", "COMPLEMENT", "APP")

# 🌟 关键锁 1：提取当前对象中所有细胞类型的固定顺序
fixed_order <- levels(cellchat@idents)

for (pw in key_pathways) {
  # 检查该通路是否存在于当前分析中
  if (pw %in% cellchat@netP$pathways) {
    file_name <- paste0("Fig_PanelD_", pw, "_Chord.pdf")
    pdf(file.path(OUT_DIR, file_name), width = 7, height = 7)
    par(mar = c(1, 1, 2, 1))
    
    # 绘制特定通路的弦图，带上新颜色并固定位置
    netVisual_aggregate(cellchat, 
                        signaling = pw, 
                        layout = "chord", 
                        color.use = color_use, 
                        title.name = paste0(pw, " signaling pathway"),
                        reduce = FALSE,         # 🌟 关键锁 2：保留所有细胞群，哪怕没有互作
                        order = fixed_order)    # 🌟 关键锁 3：严格按照固定顺序排列
    dev.off()
  } else {
    message("⚠️ 提示: 当前对象中未检测到 ", pw, " 通路，已跳过。")
  }
}

message("\n✅ 搞定！位置固定的超高颜值主图都已经存到 `Figures/Step08_CellChat_HQ` 文件夹下啦！")




# ==============================================================================
# ─── 接上文：补充宏观全景图与单向网络图 (使用相同的顶级 SCI 配色) ─────────────
# ==============================================================================

# 计算每个细胞群的大小，用于后续控制圆圈图中基座的宽度
groupSize <- as.numeric(table(cellchat@idents))

# ─── 7. 补充输出：全局网络互作圆圈图 (Count & Weight) ─────────────────────────
message("🎨 正在重绘 总体全局互作圆圈图 (Count & Weight)...")

pdf(file.path(OUT_DIR, "Fig_Extra_GlobalNet_Circle.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
# 左图：互作数量 (Count)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, color.use = color_use,
                 title.name = "Number of interactions")
# 右图：互作强度 (Weight/Strength)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, color.use = color_use,
                 title.name = "Interaction strength")
dev.off()


# ─── 8. 补充输出：全局网络互作热图 (Count & Weight) ───────────────────────────
message("🎨 正在重绘 总体全局互作热图 (Count & Weight)...")

# 使用同样的配色生成热图
p_heatmap_count <- netVisual_heatmap(cellchat, measure = "count", color.use = color_use)
p_heatmap_weight <- netVisual_heatmap(cellchat, measure = "weight", color.use = color_use)

pdf(file.path(OUT_DIR, "Fig_Extra_GlobalNet_Heatmap.pdf"), width = 12, height = 5)
print(p_heatmap_count + p_heatmap_weight)
dev.off()


# ─── 9. 补充输出：C3-high 单向发送与接收圆圈图 ────────────────────────────────
message("🎨 正在重绘 C3-high 单向作为 Sender / Receiver 的圆圈图...")

# 为了方便排版，把这两张图拼在一个 PDF 的左右两边
pdf(file.path(OUT_DIR, "Fig_Extra_C3high_Sender_Receiver.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))

# 左图：只看 C3-high 发送给别人的信号 (Sender)
netVisual_circle(cellchat@net$weight, sources.use = c3high_idx, vertex.weight = groupSize,
                 weight.scale = TRUE, color.use = color_use,
                 title.name = "C3-high Myeloid as Sender")

# 右图：只看 C3-high 接收别人的信号 (Receiver)
netVisual_circle(cellchat@net$weight, targets.use = c3high_idx, vertex.weight = groupSize,
                 weight.scale = TRUE, color.use = color_use,
                 title.name = "C3-high Myeloid as Receiver")
dev.off()

message("✅ 完美！全局视图和 C3-high 单向图均已使用高级配色重新保存至 `Figures/Step08_CellChat_HQ`。")








# ==============================================================================
# ─── 补充代码：修改 Top 50 气泡图 & 放宽 nSMC 显示阈值 ──────────────────────
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# 确保 tumor_idx 和 c3high_idx 存在
tumor_types <- c("myeSC", "nSMC")
tumor_idx <- which(levels(cellchat@idents) %in% tumor_types)
c3high_idx <- which(levels(cellchat@idents) == "C3high_Myeloid")

# 重新提取上下游通讯数据
df_lr_c3_tumor <- subsetCommunication(cellchat, sources.use = c("C3high_Myeloid"), targets.use = tumor_types)
df_lr_tumor_c3 <- subsetCommunication(cellchat, sources.use = tumor_types, targets.use = c("C3high_Myeloid"))

df_lr_all <- bind_rows(
  df_lr_c3_tumor %>% mutate(direction = "C3high -> Tumor"),
  df_lr_tumor_c3 %>% mutate(direction = "Tumor -> C3high")
)

# ──────────────────────────────────────────────────────────────────────────────
# 🎯 动作 1：把你原来的代码 n=30 改成 n=50，生成真正的综合排名图
# ──────────────────────────────────────────────────────────────────────────────
message("🎨 正在绘制真正的 Top 50 综合气泡图...")

top_lr_50 <- df_lr_all %>%
  arrange(desc(prob)) %>%
  slice_head(n = 30) %>%      # 🌟 这里正式把 n 调成了 50！
  mutate(pair = paste0(ligand, " -> ", receptor))

p_top50_lr <- ggplot(top_lr_50, aes(x = direction, y = reorder(pair, prob), size = prob, color = pathway_name)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(2, 8), name = "Interaction prob") +
  theme_bw(base_size = 12) +
  labs(title = "Top 50 L-R Interactions: C3-high Myeloid <-> Tumor",
       subtitle = "Ranked by probability", x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        axis.text.y = element_text(size = 9))

ggsave(file.path(OUT_DIR, "Fig_Extra_Top30_CustomBubble.pdf"), p_top50_lr, width = 6, height = 8)

# 
# # ──────────────────────────────────────────────────────────────────────────────
# # 🎯 动作 2：放宽 CellChat 自带图的阈值，强行逼出 nSMC
# # ──────────────────────────────────────────────────────────────────────────────
# message("🎨 正在放宽阈值重绘 CellChat 原生气泡图...")
# 
# if (length(tumor_idx) > 0) {
#   # 肿瘤 -> C3-high (放宽 p 值)
#   p_bubble_recv_relaxed <- netVisual_bubble(
#     cellchat, 
#     sources.use = tumor_idx, 
#     targets.use = c3high_idx, 
#     remove.isolate = FALSE, 
#     angle.x = 45,
#     thresh = 0.1   # 🌟 核心：默认是 0.05，我们放宽到 0.1，允许更多微弱互作显现
#   ) + 
#     ggtitle("Tumor -> C3-high Myeloid (Relaxed P<0.1)") + 
#     theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "#3C5488"))
#   
#   ggsave(file.path(OUT_DIR, "Fig_Extra_CellChat_Bubble_nSMC_Relaxed.pdf"), p_bubble_recv_relaxed, width = 8, height = 10)
# }
# 
# message("✅ 搞定！两个版本的 Top50 / 放宽版气泡图已保存。")


# 专门过滤：只看 C3high 和 Tumor 之间的总体通讯强度热图
focus_cells <- c("C3high_Myeloid", "myeSC", "nSMC")

p_custom_heatmap <- netVisual_heatmap(cellchat, 
                                      measure = "weight", 
                                      color.use = color_use,
                                      sources.use = focus_cells,   # 🌟 强行只显示这三个细胞发送的
                                      targets.use = focus_cells,   # 🌟 强行只显示这三个细胞接收的
                                      title.name = "Interaction strength (Focus)")

pdf(file.path(OUT_DIR, "Fig_Extra_Focus_Heatmap.pdf"), width = 6, height = 5)
print(p_custom_heatmap)
dev.off()










# ==============================================================================
# ─── 终极方案：手动提取矩阵绘制 定制化信号通路角色热图 (修正拼写版) ─────────
# ==============================================================================
library(ComplexHeatmap)
library(circlize)

message("🎨 正在提取底层矩阵并绘制高度定制版热图...")

# 1. 🌟 注意这里：已经把 nSMC 修正为了 nmSC ！！！, "Endothelial", "TC" "C3low_Myeloid",
focus_idents <- c("C3high_Myeloid", "myeSC", "nmSC")

# 2. 从 CellChat 的 3D 概率矩阵中，计算细胞在各个通路的发送(out)和接收(in)总和
out_mat <- t(apply(cellchat@netP$prob, c(1, 3), sum)) 
in_mat  <- t(apply(cellchat@netP$prob, c(2, 3), sum))

# 🌟 安全锁：只提取确实存在于矩阵中的细胞，彻底杜绝“下标出界”报错
focus_idents <- intersect(focus_idents, colnames(out_mat))

# 3. 计算 Top 40 核心通路
pathway_contrib <- rowSums(out_mat)
top_pathways <- names(sort(pathway_contrib, decreasing = TRUE))[1:min(10, length(pathway_contrib))]

# 提取子矩阵 (加上 drop = FALSE 防止单列变成向量)
out_mat_sub <- out_mat[top_pathways, focus_idents, drop = FALSE]
in_mat_sub  <- in_mat[top_pathways, focus_idents, drop = FALSE]

# 4. 构建顶部条形图
focus_colors <- color_use[focus_idents]

top_anno_out <- HeatmapAnnotation(
  Strength = anno_barplot(colSums(out_mat_sub), gp = gpar(fill = focus_colors, col = focus_colors)),
  show_annotation_name = FALSE
)
top_anno_in <- HeatmapAnnotation(
  Strength = anno_barplot(colSums(in_mat_sub), gp = gpar(fill = focus_colors, col = focus_colors)),
  show_annotation_name = FALSE
)

# 5. 定义渐变颜色 (经典的白->深绿色)
max_val <- max(c(out_mat_sub, in_mat_sub))
color_func <- colorRamp2(c(0, max_val), c("white", "#1B7837"))

# 6. 开始画图！
pdf(file.path(OUT_DIR, "Fig_Extra_SignalingRole_CustomMatrix.pdf"), width = 10, height = 5)

ht1 <- Heatmap(out_mat_sub,
               name = "Outgoing",
               col = color_func,
               top_annotation = top_anno_out,
               cluster_rows = FALSE, cluster_columns = FALSE,
               column_title = "Outgoing signaling (Focus)",
               row_names_side = "left",
               column_names_rot = 45,
               rect_gp = gpar(col = "white", lwd = 0.5)) 

ht2 <- Heatmap(in_mat_sub,
               name = "Incoming",
               col = color_func,
               top_annotation = top_anno_in,
               cluster_rows = FALSE, cluster_columns = FALSE,
               column_title = "Incoming signaling (Focus)",
               column_names_rot = 45,
               rect_gp = gpar(col = "white", lwd = 0.5))

# 把左右两张图拼在一起
draw(ht1 + ht2)
dev.off()

message("✅ 完美！绝对自由的定制版热图已保存。")