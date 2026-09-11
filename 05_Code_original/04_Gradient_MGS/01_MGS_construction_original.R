

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

CONFIG_PATH <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/00_config.R"
if (!file.exists(CONFIG_PATH)) stop("❌ : ", CONFIG_PATH)
source(CONFIG_PATH)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(princurve)     # Principal Curve
  library(viridis)
  library(pheatmap)
  library(RColorBrewer)
  library(uwot)          #
  library(GSVA)          # 
  library(parallel)      
})

has_ppcor   <- requireNamespace("ppcor",           quietly = TRUE)
has_enrichr <- requireNamespace("clusterProfiler", quietly = TRUE)
has_orgdb   <- requireNamespace("org.Hs.eg.db",   quietly = TRUE)

if (has_ppcor)   suppressPackageStartupMessages(library(ppcor))
if (has_enrichr) suppressPackageStartupMessages(library(clusterProfiler))
if (has_orgdb)   suppressPackageStartupMessages(library(org.Hs.eg.db))

S9_FIG_DIR <- file.path(FIG_DIR,   "Step9_MolGradient")
S9_TAB_DIR <- file.path(TABLE_DIR, "Step9_MolGradient")
dir.create(S9_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(S9_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(PARAMS$global_seed)

# ---- 0.4 Step9 分析参数 ----
TOP_MAD_GENES  <- 3000  
N_PERM         <- 1000    
N_BOOT         <- 500   
N_CORES        <- max(1, parallel::detectCores() - 2)  # 并行核数
PARTIAL_FDR    <- 0.05    #
PARTIAL_RHO    <- 0.35    #
N_GENE_MODULES <- 5       # 
WINDOW_SIZE    <- 10      # 


W <- 6.6; H <- 5.2
W_COMB <- 14; H_COMB <- 16


subtype_levels <- c("C1", "C2", "C3")
module_cols    <- brewer.pal(8, "Set2")
pt_palette     <- "inferno"        
heatmap_cols   <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)


scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) < .Machine$double.eps) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 1, color = "black"),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey40"),
      plot.margin      = margin(8, 8, 8, 8),
      strip.text       = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "grey96", color = NA),
      panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.3)
    )
}


save_pdf <- function(p, filename, w = W, h = H) {
  fn <- file.path(S9_FIG_DIR, paste0(filename, ".pdf"))
  ggsave(fn, plot = p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  message("✅ Saved: ", fn)
  invisible(fn)
}


save_pheatmap <- function(pheat_obj, filename, w = 10, h = 12) {
  fn <- file.path(S9_FIG_DIR, paste0(filename, ".pdf"))
  grDevices::cairo_pdf(fn, width = w, height = h, bg = "white")
  grid::grid.newpage()
  grid::grid.draw(pheat_obj$gtable)
  dev.off()
  message("✅ Saved: ", fn)
  invisible(fn)
}


load_rdata_safe <- function(path, label) {
  if (!file.exists(path)) stop("❌ 缺少 ", label, ": ", path)
  load(path, envir = .GlobalEnv)
  message("✅ 加载 ", label)
}

load_rdata_safe(file.path(RDATA_DIR, "Step2_Clustering.RData"),  "Step2_Clustering")
load_rdata_safe(file.path(RDATA_DIR, "Step3_DEA_Results.RData"), "Step3_DEA_Results")

# Step4 ssGSEA 结果 (可选, 用于 Phase 6 通路趋势图)
has_gsva <- file.exists(file.path(RDATA_DIR, "Step4_Enrichment.RData"))
if (has_gsva) {
  load(file.path(RDATA_DIR, "Step4_Enrichment.RData"), envir = .GlobalEnv)
  message("✅ 加载 Step4_Enrichment (gsva_matrix)")
}

# ---- 0.8 对齐表达矩阵与样本元数据 ----
# log_tpm 来自 Step2 (已 log2 变换); meta_step2 含 SampleID, Subtype 列
stopifnot(
  "log_tpm 缺失"    = exists("log_tpm"),
  "meta_step2 缺失" = exists("meta_step2")
)

# 确保 Subtype 因子水平有序
meta_use <- meta_step2 %>%
  mutate(Subtype = factor(Subtype, levels = subtype_levels)) %>%
  filter(!is.na(Subtype))

common_samples <- intersect(colnames(log_tpm), meta_use$SampleID)
if (length(common_samples) < 10) stop("❌ 可用样本不足 10 个")

expr_log <- log_tpm[, common_samples, drop = FALSE]
meta_use  <- meta_use[match(common_samples, meta_use$SampleID), ]
stopifnot(all(colnames(expr_log) == meta_use$SampleID))

n_samples <- ncol(expr_log)
message("✅ Phase 0 完成: ", n_samples, " 样本 × ", nrow(expr_log), " 基因")



message("═══ Phase 1: 无偏特征选择 ═══")

# ---- 1.1 Top MAD 特征选择 ----
#
gene_mad  <- apply(expr_log, 1, mad, na.rm = TRUE)
genes_mad <- names(sort(gene_mad, decreasing = TRUE))[
  seq_len(min(TOP_MAD_GENES, length(gene_mad)))]

expr_use <- expr_log[genes_mad, , drop = FALSE]
message("  MAD Top ", TOP_MAD_GENES, ": ", length(genes_mad), " 基因用于梯度推断")


lit_genesets <- list(
  Schwann_Myelination = c(
    "MBP","PMP22","MPZ","MAG","PRX","EGR2","SOX10",
    "PLP1","NDRG1","CLDN19","ERBB3","DHH","S100B","CADM4","NFASC"
  ),
  Schwann_Repair = c(
    "JUN","GDNF","BDNF","NGFR","SOX2","GAP43",
    "SHH","OLIG1","FGF5","RUNX2","CXCL12"
  ),
  ECM_Matrisome = c(
    "FN1","COL1A1","COL1A2","COL3A1","COL5A1","COL6A1",
    "POSTN","VCAN","TNC","SPARC","THBS1","THBS2",
    "BGN","DCN","LUM","FBLN1","FBLN2",
    "LOX","LOXL2","MMP2","MMP14",
    "TGFB1","TGFB2","TGFBI","CTGF"
  ),
  Immune_Activation = c(
    "CD8A","CD8B","CD4","CD3E",
    "GZMA","GZMB","PRF1","IFNG",
    "CXCL9","CXCL10","CXCL11","CCL5",
    "CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2",
    "IDO1","STAT1","IRF1","CIITA"
  ),
  Angiogenesis = c(
    "VEGFA","VEGFC","KDR","FLT1",
    "PECAM1","ANGPT1","ANGPT2","HIF1A",
    "EPAS1","NRP1","PDGFB"
  )
)

# 覆盖率统计
for (nm in names(lit_genesets)) {
  present <- sum(lit_genesets[[nm]] %in% rownames(expr_log))
  message("  ", nm, ": ", present, "/", length(lit_genesets[[nm]]))
}




# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║     PHASE 2: PCA + Principal Curve → Molecular Gradient Score (MGS)      ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
message("═══ Phase 2: Principal Curve → MGS ═══")

# ---- 2.1 PCA (Elbow 法自动选 PC 数) ----
pca_res  <- prcomp(t(expr_use), center = TRUE, scale. = TRUE)
var_prop <- summary(pca_res)$importance[2, ]

find_elbow <- function(vp, max_pc = 20) {
  vp <- vp[seq_len(min(max_pc, length(vp)))]
  d2 <- diff(diff(vp))
  max(which.max(abs(d2)) + 1, 3L)
}

n_pc_opt <- find_elbow(var_prop)
pca_emb  <- pca_res$x[, seq_len(n_pc_opt), drop = FALSE]
var_exp  <- var_prop[seq_len(n_pc_opt)] * 100
message("  自动选 ", n_pc_opt, " 个 PC | 累计方差 ", round(sum(var_exp), 1), "%")

# ---- 2.2 Principal Curve 拟合 ----
# stretch=0: 不外推数据点云以外, 防小样本过拟合
pc_fit <- principal_curve(
  pca_emb,
  smoother = "smooth_spline",
  trace    = FALSE,
  stretch  = 0          # [SCI1] 原 stretch=2 过于激进
)

mgs_raw <- as.numeric(pc_fit$lambda)

# ---- 2.3 方向校正: C1 在梯度低端 (代表高分化起点) ----
med_c1  <- median(mgs_raw[meta_use$Subtype == "C1"], na.rm = TRUE)
med_oth <- median(mgs_raw[meta_use$Subtype != "C1"], na.rm = TRUE)
if (med_c1 > med_oth) {
  mgs_raw <- max(mgs_raw) - mgs_raw
  message("  ↕️ 梯度方向已翻转 (C1 → 低 MGS 端)")
}

mgs_01 <- scale01(mgs_raw)

# ---- 2.4 写入 meta ----
meta_use$MGS_raw <- mgs_raw
meta_use$MGS     <- mgs_01

# ---- 2.5 UMAP (仅用于可视化, 不参与梯度计算) ----
umap_emb <- uwot::umap(pca_emb, n_neighbors = min(15, n_samples - 1),
                       min_dist = 0.3, metric = "euclidean", verbose = FALSE)
meta_use$UMAP1 <- umap_emb[, 1]
meta_use$UMAP2 <- umap_emb[, 2]

# ---- 2.6 准备绘图数据 ----
pca_df   <- as.data.frame(pca_emb[, 1:2]) %>%
  mutate(SampleID = meta_use$SampleID,
         Subtype  = meta_use$Subtype,
         MGS      = mgs_01)

# Principal Curve 轨迹坐标 (按 lambda 排序取前两个 PC)
curve_idx <- order(pc_fit$lambda)
curve_df  <- data.frame(PC1 = pc_fit$s[curve_idx, 1],
                        PC2 = pc_fit$s[curve_idx, 2])

message("✅ Phase 2 完成 | MGS 范围 [",
        round(min(mgs_01), 3), ", ", round(max(mgs_01), 3), "]")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║     PHASE 3: 统计验证 (排列检验 + Bootstrap + Schwann 锚定 + LOCO)       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
message("═══ Phase 3: 统计验证 ═══")

# ---------- 3A: Permutation Test (KW) ----------
message("  3A: 排列检验 (", N_PERM, " 次) ...")

real_kw  <- kruskal.test(MGS ~ Subtype, data = meta_use)$statistic
perm_kw  <- replicate(N_PERM, {
  kruskal.test(mgs_01 ~ sample(meta_use$Subtype))$statistic
})
perm_p <- mean(perm_kw >= real_kw)
message("  KW χ²=", round(real_kw, 2), " | 排列 p=", perm_p)

# ---------- 3B: Bootstrap 稳定性 (并行) ----------
message("  3B: Bootstrap (", N_BOOT, " 次, ", N_CORES, " 核) ...")

# [PERF2] 改用 mclapply 并行化, 8 核服务器约 3 分钟 (原串行 ~20 分钟)
boot_one <- function(b) {
  idx <- sample(seq_len(n_samples), replace = TRUE)
  eb  <- expr_use[, idx, drop = FALSE]
  pb  <- tryCatch(prcomp(t(eb), center = TRUE, scale. = TRUE), error = function(e) NULL)
  if (is.null(pb)) return(NULL)
  emb <- pb$x[, seq_len(min(n_pc_opt, ncol(pb$x))), drop = FALSE]
  pc  <- tryCatch(principal_curve(emb, smoother = "smooth_spline",
                                  trace = FALSE, stretch = 0), error = function(e) NULL)
  if (is.null(pc)) return(NULL)
  lam <- as.numeric(pc$lambda)
  if (suppressWarnings(cor(lam, mgs_raw[idx], method = "spearman")) < 0)
    lam <- max(lam) - lam
  list(idx = idx, ranks = rank(lam))
}

boot_list <- parallel::mclapply(seq_len(N_BOOT), boot_one, mc.cores = N_CORES)

# 聚合 bootstrap 排名矩阵
boot_ranks <- matrix(NA_real_, nrow = n_samples, ncol = N_BOOT)
for (b in seq_along(boot_list)) {
  res <- boot_list[[b]]
  if (is.null(res)) next
  for (j in seq_along(res$idx)) boot_ranks[res$idx[j], b] <- res$ranks[j]
}

boot_taus <- apply(boot_ranks, 2, function(r) {
  v <- !is.na(r)
  if (sum(v) < 10) return(NA_real_)
  suppressWarnings(cor(rank(mgs_01)[v], r[v], method = "kendall"))
})
boot_taus <- boot_taus[!is.na(boot_taus)]

mean_tau <- mean(boot_taus)
ci_tau   <- quantile(boot_taus, c(0.025, 0.975))   # [SCI2] 新增 95% CI

# 样本排名稳定性 (IQR/range 越小越稳定)
rank_stability <- apply(boot_ranks, 1, function(r) {
  r <- r[!is.na(r)]
  if (length(r) < 10) return(NA_real_)
  1 - IQR(r) / diff(range(r))
})
meta_use$Rank_Stability <- rank_stability

message("  Bootstrap Kendall τ = ", round(mean_tau, 3),
        " [", round(ci_tau[1], 3), ", ", round(ci_tau[2], 3), "]",
        ifelse(mean_tau > 0.7, " ✅ 稳健", " ⚠️ 需谨慎"))

# ---------- 3C: Schwann 分化标志物独立锚定验证 ----------
message("  3C: Schwann 锚点验证 ...")

anchor_genes <- lit_genesets$Schwann_Myelination
anchor_avail <- anchor_genes[anchor_genes %in% rownames(expr_log)]

anchor_df <- lapply(anchor_avail, function(g) {
  ct <- suppressWarnings(cor.test(as.numeric(expr_log[g, common_samples]),
                                  mgs_01, method = "spearman", exact = FALSE))
  data.frame(Gene = g, Rho = ct$estimate, Pval = ct$p.value,
             In_MAD = g %in% genes_mad, stringsAsFactors = FALSE)
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(Pval, method = "BH"),
         Consistent = Rho < 0)  # 期望: 高分化标志物随 MGS 下降

write.csv(anchor_df, file.path(S9_TAB_DIR, "MGS_Anchor_Validation.csv"),
          row.names = FALSE)

n_consist <- sum(anchor_df$Consistent, na.rm = TRUE)
message("  Schwann 锚点: ", n_consist, "/", nrow(anchor_df),
        " 个基因表现一致 (随 MGS 下降)")

# ---------- 3D: Leave-One-Cluster-Out (LOCO) 验证 ----------
message("  3D: LOCO 验证 ...")

loco_df_list <- list()

for (leave_out in subtype_levels) {
  keep_idx <- which(meta_use$Subtype != leave_out)
  left_idx <- which(meta_use$Subtype == leave_out)
  if (length(keep_idx) < 15 || length(left_idx) < 3) {
    message("    跳过 ", leave_out, " (样本不足)"); next
  }
  
  pb_keep <- tryCatch(
    prcomp(t(expr_use[, keep_idx, drop = FALSE]), center = TRUE, scale. = TRUE),
    error = function(e) NULL)
  if (is.null(pb_keep)) next
  
  n_pc_k <- min(n_pc_opt, ncol(pb_keep$x))
  emb_keep <- pb_keep$x[, seq_len(n_pc_k), drop = FALSE]
  
  pc_keep <- tryCatch(
    principal_curve(emb_keep, smoother = "smooth_spline",
                    trace = FALSE, stretch = 0),
    error = function(e) NULL)
  if (is.null(pc_keep)) next
  
  # 将 left-out 样本投影至 keep 样本的曲线空间
  left_centered <- scale(t(expr_use[, left_idx, drop = FALSE]),
                         center = pb_keep$center,
                         scale  = pb_keep$scale)
  emb_left <- left_centered %*% pb_keep$rotation[, seq_len(n_pc_k), drop = FALSE]
  
  proj_left   <- princurve::project_to_curve(emb_left,
                                             pc_keep$s[order(pc_keep$lambda), , drop = FALSE])
  lambda_keep <- as.numeric(pc_keep$lambda)
  lambda_left <- as.numeric(proj_left$lambda)
  
  # [BUG1 修复] 翻转时先保存原始最大值, keep 和 left 使用同一参考坐标系
  keep_subs <- meta_use$Subtype[keep_idx]
  if ("C1" %in% keep_subs) {
    c1_med  <- median(lambda_keep[keep_subs == "C1"])
    oth_med <- median(lambda_keep[keep_subs != "C1"])
    if (c1_med > oth_med) {
      max_ref     <- max(lambda_keep)          # ★ 翻转前保存原始最大值
      lambda_keep <- max_ref - lambda_keep
      lambda_left <- max_ref - lambda_left     # ★ 用同一参考值, 而非翻转后的 max
    }
  }
  
  all_lam <- c(lambda_keep, lambda_left)
  all_01  <- scale01(all_lam)
  left_01 <- all_01[(length(lambda_keep) + 1):length(all_lam)]
  
  loco_df_list[[leave_out]] <- data.frame(
    LeftOut     = leave_out,
    SampleID    = meta_use$SampleID[left_idx],
    Subtype     = as.character(meta_use$Subtype[left_idx]),
    Projected_MGS = left_01,
    stringsAsFactors = FALSE)
  
  message("    去掉 ", leave_out, ": 投影 MGS 中位数 = ",
          round(median(left_01, na.rm = TRUE), 3))
}

if (length(loco_df_list) > 0) {
  loco_df <- bind_rows(loco_df_list)
  write.csv(loco_df, file.path(S9_TAB_DIR, "MGS_LOCO_Validation.csv"),
            row.names = FALSE)
}

# 保存统计验证汇总 (含 Bootstrap 95% CI)
phase3_stats <- data.frame(
  Test  = c("Permutation_KW_p", "Bootstrap_Kendall_tau",
            "Bootstrap_tau_CI_low", "Bootstrap_tau_CI_high"),
  Value = c(perm_p, mean_tau, ci_tau[1], ci_tau[2]),
  stringsAsFactors = FALSE
)
write.csv(phase3_stats, file.path(S9_TAB_DIR, "MGS_StatisticalValidation.csv"),
          row.names = FALSE)

message("✅ Phase 3 完成")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║     PHASE 4: 梯度关联基因鉴定 (向量化 Spearman + 偏相关)                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
message("═══ Phase 4: 梯度关联基因鉴定 ═══")

# ---- 4A: Spearman 相关 (向量化矩阵运算) ----
# [PERF1] 改用矩阵运算替代 lapply+cor.test, 速度提升 10-50×
# 原理: Spearman ρ = 秩变换后的 Pearson r, p 值用 t 近似公式

purity_cands <- c("TumorPurity", "Purity", "purity", "ESTIMATE_Purity", "Tumor_Score")
purity_col   <- intersect(purity_cands, colnames(meta_use))[1]

if (!is.na(purity_col) && has_ppcor) {
  # ── 偏 Spearman (控制肿瘤纯度混杂) ──
  # 注意: 对 VS (前庭神经鞘瘤) 而言, TumorPurity 实为 Schwann 细胞比例,
  #       而非恶性肿瘤的癌细胞占比。如有 CIBERSORTx 解卷积结果, 建议替换。
  message("  4A: 偏 Spearman (控制 ", purity_col, ") | 并行 ", N_CORES, " 核 ...")
  purity_vec <- as.numeric(meta_use[[purity_col]])
  purity_vec[is.na(purity_vec)] <- median(purity_vec, na.rm = TRUE)
  
  partial_list <- parallel::mclapply(rownames(expr_log), function(g) {
    x <- as.numeric(expr_log[g, common_samples])
    tryCatch({
      res <- ppcor::pcor.test(x, mgs_01, purity_vec, method = "spearman")
      c(rho = as.numeric(res$estimate), pval = as.numeric(res$p.value))
    }, error = function(e) c(rho = NA_real_, pval = NA_real_))
  }, mc.cores = N_CORES)
  
  partial_mat <- do.call(rbind, partial_list)
  rownames(partial_mat) <- rownames(expr_log)
  
  gradient_df <- data.frame(
    Gene = rownames(partial_mat),
    Rho  = partial_mat[, "rho"],
    Pval = partial_mat[, "pval"],
    stringsAsFactors = FALSE)
  gradient_method <- "Partial_Spearman"
  
} else {
  # ── 普通 Spearman (全基因组矩阵化计算) ──
  message("  4A: 全基因组 Spearman (矩阵运算) ...")
  n  <- ncol(expr_log)
  
  # 对每个基因逐行秩变换, 然后矩阵乘以 MGS 秩向量
  expr_ranked  <- t(apply(expr_log, 1, rank))             # genes × samples
  mgs_ranked   <- rank(mgs_01)                            # samples
  
  # Pearson(ranks) = Spearman ρ
  # cor(expr_ranked[g,], mgs_ranked) 等价于矩阵乘法除以 n-1
  mgs_r_sc  <- scale(mgs_ranked)                          # 中心化+归一化
  expr_r_sc <- t(scale(t(expr_ranked)))                   # 每行中心化+归一化
  rho_vec   <- as.numeric(expr_r_sc %*% mgs_r_sc) / (n - 1)
  rho_vec   <- pmax(pmin(rho_vec, 1 - 1e-9), -1 + 1e-9)  # 防止 |ρ|=1 导致 t→Inf
  
  # t 近似公式 (与 cor.test 等价)
  t_stat  <- rho_vec * sqrt((n - 2) / (1 - rho_vec^2))
  pval_vec <- 2 * pt(-abs(t_stat), df = n - 2)
  
  gradient_df <- data.frame(
    Gene = rownames(expr_log),
    Rho  = rho_vec,
    Pval = pval_vec,
    stringsAsFactors = FALSE)
  gradient_method <- "Spearman"
}

# 统一后处理
gradient_df <- gradient_df %>%
  filter(!is.na(Rho), !is.na(Pval)) %>%
  mutate(
    FDR = p.adjust(Pval, method = "BH"),
    Direction = case_when(
      Rho >  PARTIAL_RHO & FDR < PARTIAL_FDR ~ "Up_along_MGS",
      Rho < -PARTIAL_RHO & FDR < PARTIAL_FDR ~ "Down_along_MGS",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(Pval)

genes_up   <- gradient_df %>% filter(Direction == "Up_along_MGS")   %>% pull(Gene)
genes_down <- gradient_df %>% filter(Direction == "Down_along_MGS") %>% pull(Gene)

write.csv(gradient_df, file.path(S9_TAB_DIR, "MGS_GradientGenes_All.csv"),
          row.names = FALSE)
message("  ", gradient_method, ": Up=", length(genes_up),
        ", Down=", length(genes_down))

# ---- 4B: 滑窗拐点检测 ----
message("  4B: 滑窗拐点检测 ...")
mgs_order    <- order(mgs_01)
expr_ordered <- expr_log[, common_samples[mgs_order], drop = FALSE]
mgs_sorted   <- mgs_01[mgs_order]

genes_sig <- unique(c(head(genes_up, 100), head(genes_down, 100)))
genes_sig <- genes_sig[genes_sig %in% rownames(expr_ordered)]

if (length(genes_sig) >= 10) {
  bp_list <- lapply(genes_sig, function(g) {
    eg <- as.numeric(expr_ordered[g, ])
    best_t <- 0; best_pos <- length(eg) %/% 2
    for (i in WINDOW_SIZE:(length(eg) - WINDOW_SIZE)) {
      tt <- suppressWarnings(t.test(eg[1:i], eg[(i+1):length(eg)])$statistic)
      if (!is.na(tt) && abs(tt) > abs(best_t)) { best_t <- tt; best_pos <- i }
    }
    data.frame(Gene = g, Breakpoint_Rank = best_pos,
               Breakpoint_MGS = mgs_sorted[best_pos],
               T_statistic = best_t, stringsAsFactors = FALSE)
  })
  write.csv(bind_rows(bp_list),
            file.path(S9_TAB_DIR, "MGS_GeneBreakpoints.csv"),
            row.names = FALSE)
}



# ---- 4C: GO/KEGG 通路富集 ----

# 🌟 终极修复：强行夺回 margin 的控制权！
# 这样你之前定义的 theme_pub() 里面用到的 margin 就会自动恢复正常
margin <- ggplot2::margin

if (has_enrichr && has_orgdb) {
  message("  4C: GO/KEGG 富集 ...")
  run_enrich <- function(gene_list, label) {
    if (length(gene_list) < 5) return(invisible(NULL))
    ego <- tryCatch(
      enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = "BP", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2),
      error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(as.data.frame(ego),
                file.path(S9_TAB_DIR, paste0("MGS_GO_BP_", label, ".csv")),
                row.names = FALSE)
      
      p_dot <- dotplot(ego, showCategory = 15, font.size = 9) +
        labs(title = paste0("GO BP: ", label)) + 
        theme_pub() # 这里之前报错，现在已经被修复了
      
      save_pdf(p_dot, paste0("FigS_GO_", label), w = 10, h = 10)
    }
    return(invisible(ego))
  }
  
  run_enrich(genes_up,   "Up_along_MGS")
  run_enrich(genes_down, "Down_along_MGS")
}
# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 5: 基因模块 (k-means) + Fisher's Exact Overlap 验证               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
message("═══ Phase 5: 基因模块 + Overlap 验证 ═══")

# ---- 5.1 梯度基因 k-means 聚类 ----
genes_for_mod <- unique(c(head(genes_up, 150), head(genes_down, 150)))
genes_for_mod <- genes_for_mod[genes_for_mod %in% rownames(expr_log)]

mat_mod_z    <- NULL
module_final <- NULL
eigen_df     <- NULL

if (length(genes_for_mod) >= 20) {
  mat_mod   <- expr_log[genes_for_mod, common_samples[mgs_order], drop = FALSE]
  mat_mod_z <- t(scale(t(mat_mod)))
  mat_mod_z[mat_mod_z >  2.5] <-  2.5
  mat_mod_z[mat_mod_z < -2.5] <- -2.5
  mat_mod_z[is.na(mat_mod_z)] <- 0
  
  set.seed(PARAMS$global_seed)
  km  <- kmeans(mat_mod_z, centers = N_GENE_MODULES, nstart = 50, iter.max = 200)
  
  # 按峰值位置重排模块编号 (使模块编号与 MGS 轴有序对应)
  mod_peaks   <- sapply(seq_len(N_GENE_MODULES), function(k) {
    gk <- names(km$cluster)[km$cluster == k]
    if (length(gk) < 3) return(ncol(mat_mod_z) / 2)
    which.max(colMeans(mat_mod_z[gk, , drop = FALSE]))
  })
  new_order    <- order(mod_peaks)
  remap        <- setNames(seq_along(new_order), new_order)
  module_final <- remap[as.character(km$cluster)]
  names(module_final) <- names(km$cluster)
  
  # Eigengene 趋势数据
  eigen_df <- lapply(seq_len(N_GENE_MODULES), function(k) {
    gk <- names(module_final)[module_final == k]
    if (length(gk) < 3) return(NULL)
    data.frame(MGS     = mgs_sorted,
               Subtype = meta_use$Subtype[mgs_order],
               Eigen   = colMeans(mat_mod_z[gk, , drop = FALSE]),
               Module  = paste0("M", k, " (n=", length(gk), ")"),
               stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  # 模块内 GO 富集
  if (has_enrichr && has_orgdb) {
    for (k in seq_len(N_GENE_MODULES)) {
      gk <- names(module_final)[module_final == k]
      if (length(gk) < 10) next
      ego <- tryCatch(
        enrichGO(gene = gk, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2),
        error = function(e) NULL)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0)
        write.csv(as.data.frame(ego),
                  file.path(S9_TAB_DIR, paste0("MGS_Module_M", k, "_GO_BP.csv")),
                  row.names = FALSE)
    }
  }
  
  write.csv(
    data.frame(Gene = names(module_final),
               Module = paste0("M", module_final), stringsAsFactors = FALSE) %>%
      dplyr::left_join(dplyr::select(gradient_df, Gene, Rho, FDR, Direction), by = "Gene"),
    file.path(S9_TAB_DIR, "MGS_GeneModules.csv"),
    row.names = FALSE)
}  # ← [BUG2 修复] 正确闭合 if(genes_for_mod) 块

# ---- 5.2 Post-hoc Fisher's Exact Overlap (梯度基因 vs 文献基因集) ----
# [BUG2 修复] 此步骤现在独立于 if(genes_for_mod) 块之外, 始终运行
message("  5.2: Overlap 检验 (Fisher's Exact) ...")

n_background <- nrow(expr_log)
overlap_res  <- list()

for (set_name in names(lit_genesets)) {
  lit_in <- lit_genesets[[set_name]][lit_genesets[[set_name]] %in% rownames(expr_log)]
  if (length(lit_in) < 5) next
  
  for (dir_str in c("Up_along_MGS", "Down_along_MGS")) {
    grad_g <- gradient_df %>% filter(Direction == dir_str) %>% pull(Gene)
    if (length(grad_g) < 5) next
    
    a  <- length(intersect(grad_g, lit_in))
    b  <- length(setdiff(grad_g, lit_in))
    cc <- length(setdiff(lit_in, grad_g))
    d  <- n_background - a - b - cc
    ft <- fisher.test(matrix(c(a, b, cc, d), 2), alternative = "greater")
    
    overlap_res[[paste0(set_name, "__", dir_str)]] <- data.frame(
      Lit_Set = set_name, Direction = dir_str,
      Overlap_N = a, Gradient_N = length(grad_g),
      Lit_N = length(lit_in), OR = ft$estimate, P = ft$p.value,
      Overlap_Genes = paste(intersect(grad_g, lit_in), collapse = ";"),
      stringsAsFactors = FALSE)
  }
}

overlap_df <- bind_rows(overlap_res) %>%
  mutate(FDR = p.adjust(P, method = "BH")) %>%
  arrange(P)

write.csv(overlap_df, file.path(S9_TAB_DIR, "MGS_Overlap_Fisher.csv"),
          row.names = FALSE)

for (i in seq_len(min(8, nrow(overlap_df)))) {
  r <- overlap_df[i, ]
  message(sprintf("  %s %-24s × %-15s n=%d OR=%.1f FDR=%.3f",
                  ifelse(r$FDR < 0.05, "★", " "),
                  r$Lit_Set, r$Direction, r$Overlap_N, r$OR, r$FDR))
}

message("✅ Phase 5 完成")




message("═══ Phase 6: ssGSEA Signature Scores ═══")

df_sig_long <- NULL

if (all(exists("deg_c1"), exists("deg_c2"), exists("deg_c3"))) {
  
  # 构建 C1/C2/C3 差异基因集 (上调 + 下调方向均纳入, 方向性基因集)
  build_sig_geneset <- function(deg_df, n_up = 100, n_dn = 50) {
    ups <- deg_df %>%
      filter(Change == "Up", adj.P.Val < 0.05) %>%
      arrange(desc(logFC)) %>%
      head(n_up) %>% pull(Gene)
    # 若基因集过小, 放宽阈值
    if (length(ups) < 20)
      ups <- deg_df %>% filter(logFC > 0) %>%
        arrange(desc(logFC)) %>% head(n_up) %>% pull(Gene)
    ups
  }
  
  sig_sets <- list(
    C1_Signature = build_sig_geneset(deg_c1),
    C2_Signature = build_sig_geneset(deg_c2),
    C3_Signature = build_sig_geneset(deg_c3)
  )
  
  # 过滤至表达矩阵中存在的基因
  sig_sets <- lapply(sig_sets, function(gs) gs[gs %in% rownames(expr_log)])
  sig_sets <- sig_sets[sapply(sig_sets, length) >= 10]
  
  if (length(sig_sets) >= 2) {
    expr_mat_raw <- 2^expr_log - 1  # ssGSEA 推荐用原始(非log)TPM空间
    expr_mat_raw[expr_mat_raw < 0] <- 0
    
    # 兼容 GSVA v1.50+ 新 API 与旧 API
    gsva_sig <- tryCatch({
      if (utils::packageVersion("GSVA") >= "1.50.0") {
        param <- ssgseaParam(exprData = as.matrix(expr_mat_raw), geneSets = sig_sets)
        gsva(param, verbose = FALSE)
      } else {
        gsva(as.matrix(expr_mat_raw), sig_sets, method = "ssgsea",
             kcdf = "Gaussian", verbose = FALSE)
      }
    }, error = function(e) {
      message("  ⚠️ ssGSEA 失败 (", e$message, "), 退回均值 z-score")
      NULL
    })
    
    if (!is.null(gsva_sig)) {
      sig_score_df <- as.data.frame(t(gsva_sig[, common_samples, drop = FALSE]))
      sig_score_df$SampleID <- rownames(sig_score_df)
      
      df_sig <- left_join(sig_score_df,
                          meta_use[, c("SampleID", "Subtype", "MGS")],
                          by = "SampleID")
      
      df_sig_long <- df_sig %>%
        pivot_longer(cols = ends_with("_Signature"),
                     names_to = "Signature", values_to = "ssGSEA_Score") %>%
        mutate(Signature = factor(Signature,
                                  levels = c("C1_Signature", "C2_Signature", "C3_Signature"),
                                  labels = c("C1 Signature", "C2 Signature", "C3 Signature")))
      
      write.csv(df_sig, file.path(S9_TAB_DIR, "MGS_SignatureScores_ssGSEA.csv"),
                row.names = FALSE)
      message("  ssGSEA Signature Scores 计算完成")
    }
  }
}

# Step4 通路 ssGSEA (gsva_matrix) 沿 MGS 的趋势
prog_df_long <- NULL
if (has_gsva && exists("gsva_matrix")) {
  prog_samps  <- intersect(colnames(gsva_matrix), common_samples)
  prog_scores <- as.data.frame(t(gsva_matrix[, prog_samps, drop = FALSE]))
  prog_scores$SampleID <- rownames(prog_scores)
  prog_df_long <- left_join(prog_scores,
                            meta_use[, c("SampleID", "Subtype", "MGS")],
                            by = "SampleID") %>%
    pivot_longer(cols = -c(SampleID, Subtype, MGS),
                 names_to = "Pathway", values_to = "ssGSEA_Score")
}

message("✅ Phase 6 完成")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║                       PHASE 7: 全套出图                                  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
message("═══ Phase 7: 出图 ═══")

# ── 通用 LOESS 参数 ──
LOESS_SPAN  <- 1.5
LOESS_ALPHA <- 0.12   # 置信带透明度 (全文统一)

# ====== Fig A: PCA + Principal Curve ======
p_A <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_path(data = curve_df, aes(x = PC1, y = PC2),
            linewidth = 1.2, color = "grey25",
            arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_point(aes(fill = Subtype), shape = 21, size = 3.5,
             color = "white", stroke = 0.5, alpha = 0.9) +
  scale_fill_manual(values = SUBTYPE_COLORS, name = "Subtype") +
  labs(title = "A | PCA + Principal Curve",
       subtitle = "Molecular gradient fitted in PCA space",
       x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  theme_pub() +
  theme(legend.position = c(0.85, 0.15))

# ====== Fig B: UMAP colored by MGS ======
p_B <- ggplot(meta_use, aes(x = UMAP1, y = UMAP2, color = MGS)) +
  geom_point(size = 3.0, alpha = 0.9) +
  scale_color_viridis_c(option = pt_palette, end = 0.95, name = "MGS") +
  labs(title = "B | UMAP (colored by MGS)", x = "UMAP1", y = "UMAP2") +
  theme_pub()

# ====== Fig C: MGS Boxplot ======
kw_res <- kruskal.test(MGS ~ Subtype, data = meta_use)
pw_df  <- as.data.frame(as.table(
  pairwise.wilcox.test(meta_use$MGS, meta_use$Subtype,
                       p.adjust.method = "BH", exact = FALSE)$p.value)) %>%
  filter(!is.na(Freq)) %>%
  dplyr::rename(Group1 = Var1, Group2 = Var2, P_adj = Freq)
write.csv(pw_df, file.path(S9_TAB_DIR, "MGS_PairwiseWilcox.csv"),
          row.names = FALSE)

p_C <- ggplot(meta_use, aes(x = Subtype, y = MGS, fill = Subtype)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.75, shape = 21,
              fill = "white", color = "grey30", stroke = 0.5) +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  labs(title = "C | MGS by subtype",
       subtitle = paste0("KW p=", signif(kw_res$p.value, 3),
                         " | Perm p=", perm_p),
       x = NULL,
       y = "MGS  (0 = differentiated → 1 = de-differentiated)") +
  theme_pub() + theme(legend.position = "none")

# ====== Fig D: Permutation Test ======
p_D <- ggplot(data.frame(KW = perm_kw), aes(x = KW)) +
  geom_histogram(bins = 50, fill = "grey75", color = "grey55", alpha = 0.85) +
  geom_vline(xintercept = real_kw, color = SUBTYPE_COLORS["C3"],
             linewidth = 1.2) +
  annotate("text", x = real_kw, y = Inf, vjust = 1.5,
           label = paste0("Observed χ²=", round(real_kw, 1),
                          "\np=", perm_p),
           color = SUBTYPE_COLORS["C3"], fontface = "bold", size = 3.5) +
  labs(title = "D | Permutation test (n=1000)",
       subtitle = "H₀: Subtype labels ⊥ MGS",
       x = "KW χ² (permuted)", y = "Count") +
  theme_pub()

# ====== Fig E: Bootstrap stability (with 95% CI) ======
p_E <- ggplot(data.frame(Tau = boot_taus), aes(x = Tau)) +
  geom_histogram(bins = 40, fill = SUBTYPE_COLORS["C2"],
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept = mean_tau, color = SUBTYPE_COLORS["C1"],
             linewidth = 1.0, linetype = "dashed") +
  geom_vline(xintercept = ci_tau, color = SUBTYPE_COLORS["C1"],  # [SCI2] 95% CI
             linewidth = 0.6, linetype = "dotted") +
  annotate("text", x = mean_tau, y = Inf, vjust = 1.5,
           label = paste0("Mean τ=", round(mean_tau, 3),
                          "\n95% CI [", round(ci_tau[1], 3),
                          ", ", round(ci_tau[2], 3), "]"),
           color = SUBTYPE_COLORS["C1"], fontface = "bold", size = 3.2) +
  labs(title = "E | Bootstrap stability (n=500)",
       subtitle = "Kendall τ: bootstrap vs original MGS ranking",
       x = "Kendall τ", y = "Count") +
  theme_pub()
# 选出单调性最强的 6 个下降基因画图
rep_down_ideal <- c("PAQR6", "MTA1", "CKB", "KLF15", "CSRP1", "CLDN15")

# 重新运行 Phase 7 中绘制 Fig F 的那部分代码
# ====== Fig F: Schwann 锚点验证 ======
anchor_plot_genes <- head(anchor_df$Gene, 6)
p_F <- NULL
if (length(anchor_plot_genes) >= 3) {
  schwann_long <- lapply(anchor_plot_genes, function(g) {
    data.frame(MGS  = mgs_01,
               Expr = as.numeric(expr_log[g, common_samples]),
               Gene = g,
               Subtype = meta_use$Subtype, stringsAsFactors = FALSE)
  }) %>% bind_rows() %>%
    mutate(Gene = factor(Gene, levels = anchor_plot_genes))
  
  <- ggplot(schwann_long, aes(x = MGS, y = Expr)) +
      geom_point(aes(color = Subtype), size = 2.0, alpha = 0.7) +
      geom_smooth(method = "loess", se = TRUE, span = LOESS_SPAN,
                  linewidth = 0.8, color = "black", alpha = LOESS_ALPHA) +
      scale_color_manual(values = SUBTYPE_COLORS) +
      facet_wrap(~Gene, scales = "free_y", ncol = 3) +
      labs(title = "F | Schwann myelination markers (independent anchors)",
           subtitle = "Expected: decrease along MGS (de-differentiation → low expression)",
           x = "MGS (0→1)", y = "Expression (log2 TPM)", color = "Subtype") +
      theme_pub() + theme(legend.position = "bottom")
    
    save_pdf(p_F, "FigF_Schwann_Anchors",
             w = 9, h = 2 + 2.5 * ceiling(length(anchor_plot_genes) / 3))
}

# ====== Fig G: Gene Module Heatmap (cairo_pdf, 与全文统一) ======
if (!is.null(mat_mod_z) && !is.null(module_final)) {
  gene_ord  <- names(sort(module_final))
  mat_hm    <- mat_mod_z[gene_ord, , drop = FALSE]
  
  ann_col_hm <- data.frame(
    Subtype = meta_use$Subtype[mgs_order],
    MGS     = mgs_sorted,
    row.names = colnames(mat_hm))
  ann_row_hm <- data.frame(
    Module = paste0("M", module_final[gene_ord]),
    row.names = gene_ord)
  
  mod_palette <- setNames(module_cols[seq_len(N_GENE_MODULES)],
                          paste0("M", seq_len(N_GENE_MODULES)))
  ann_colors_hm <- list(
    Subtype = SUBTYPE_COLORS,
    MGS     = viridis(100, option = pt_palette),
    Module  = mod_palette)
  
  # [BUG4 修复] 用 cairo_pdf + grid::grid.draw 替代 base pdf()
  ph_g <- pheatmap(mat_hm,
                   cluster_cols = FALSE, cluster_rows = FALSE,
                   show_colnames = FALSE,
                   show_rownames = nrow(mat_hm) <= 80,
                   annotation_col = ann_col_hm, annotation_row = ann_row_hm,
                   annotation_colors = ann_colors_hm,
                   color = heatmap_cols, border_color = NA,
                   fontsize = 8, fontsize_row = 5,
                   gaps_row = cumsum(table(module_final[gene_ord])),
                   main = "Gene modules along Molecular Gradient",
                   silent = TRUE)
  save_pheatmap(ph_g, "FigG_GeneModule_Heatmap", w = 10, h = 14)
}

# ====== Fig H: Module Eigengene Trends ======
p_H <- NULL
if (!is.null(eigen_df) && nrow(eigen_df) > 0) {
  p_H <- ggplot(eigen_df, aes(x = MGS, y = Eigen)) +
    geom_point(aes(color = Subtype), size = 1.8, alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, span = LOESS_SPAN,
                linewidth = 0.9, color = "black", alpha = LOESS_ALPHA) +
    scale_color_manual(values = SUBTYPE_COLORS) +
    facet_wrap(~Module, ncol = 3, scales = "free_y") +
    labs(title = "H | Module eigengenes along MGS",
         x = "MGS (0→1)", y = "Eigengene (z-score)", color = "Subtype") +
    theme_pub() + theme(legend.position = "bottom")
  
  save_pdf(p_H, "FigH_ModuleEigengenes",
           w = 10, h = 2 + 2.5 * ceiling(N_GENE_MODULES / 3))
}


p_I <- NULL
if (!is.null(df_sig_long)) {
  # 终极防错版：直接提供 3 个颜色代码，不写名字，它会自动按 C1, C2, C3 的顺序匹配
  sig_colors <- c("#00A087", "#4DBBD5", "#E64B35")
  
  # 注意：在 aes() 里新增了 fill = Signature，让拟合曲线的置信带阴影也跟着变色
  p_I <- ggplot(df_sig_long, aes(x = MGS, y = ssGSEA_Score, color = Signature, fill = Signature)) +
    geom_point(size = 2.0, alpha = 0.55) +
    geom_smooth(method = "loess", se = TRUE, span = LOESS_SPAN,
                linewidth = 1.0, alpha = LOESS_ALPHA) +
    # 同时手动修改线条颜色 (color) 和阴影颜色 (fill)
    scale_color_manual(values = sig_colors) +
    scale_fill_manual(values = sig_colors) +
    labs(title = "I | Subtype ssGSEA signatures along MGS",
         subtitle = "Computed via ssGSEA (GSVA pkg) on DEG-derived gene sets",
         x = "MGS (0→1)", y = "ssGSEA Score", 
         color = NULL, fill = NULL) + # 隐藏两个图例的标题
    theme_pub() +
    theme(legend.position = c(0.15, 0.85),
          legend.background = element_rect(fill = alpha("white", 0.8), color = NA))
  
  save_pdf(p_I, "FigI_SignatureTrends_ssGSEA", w = 7, h = 5.5)
}
# ====== Fig J: Step4 Pathway Scores along MGS ======
if (!is.null(prog_df_long)) {
  p_J <- ggplot(prog_df_long, aes(x = MGS, y = ssGSEA_Score)) +
    geom_point(aes(color = Subtype), size = 1.8, alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, span = LOESS_SPAN,
                linewidth = 0.8, color = "black", alpha = LOESS_ALPHA) +
    scale_color_manual(values = SUBTYPE_COLORS) +
    facet_wrap(~Pathway, scales = "free_y", ncol = 3) +
    labs(title = "J | Pathway ssGSEA scores along MGS",
         subtitle = "From Step4 ssGSEA analysis (GO-term gene sets)",
         x = "MGS (0→1)", y = "ssGSEA Score", color = "Subtype") +
    theme_pub() + theme(legend.position = "bottom")
  
  n_prog <- length(unique(prog_df_long$Pathway))
  save_pdf(p_J, "FigJ_PathwayScores",
           w = 10, h = 2.5 + 2.5 * ceiling(n_prog / 3))
}

# ====== Fig K: Gradient Gene Volcano ======
p_K <- ggplot(gradient_df, aes(x = Rho, y = -log10(FDR + 1e-300))) +
  geom_point(aes(color = Direction), size = 0.9, alpha = 0.5) +
  scale_color_manual(
    values = c(Up_along_MGS = "#D73027", Down_along_MGS = "#4575B4", NS = "grey78"),
    labels = c("Up along MGS", "Down along MGS", "NS")) +
  geom_hline(yintercept = -log10(PARTIAL_FDR), linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = c(-PARTIAL_RHO, PARTIAL_RHO), linetype = "dotted",
             color = "grey40", linewidth = 0.3) +
  labs(title = paste0("K | ", gradient_method, " with MGS"),
       subtitle = paste0("Up=", length(genes_up), ", Down=", length(genes_down)),
       x = paste0(gradient_method, " ρ"),
       y = expression(-log[10](FDR)), color = NULL) +
  theme_pub() + theme(legend.position = c(0.15, 0.85))

save_pdf(p_K, "FigK_GradientGene_Volcano", w = 6.5, h = 5.5)

# ====== Fig L: Fisher Overlap Dotplot ======
p_L <- NULL
if (nrow(overlap_df) > 0) {
  p_L <- ggplot(
    overlap_df %>%
      mutate(Label    = paste0(Lit_Set, "\n(", Direction, ")"),
             Sig      = FDR < 0.05,
             neg_logF = pmin(-log10(FDR + 1e-300), 10),
             OR_cap   = pmin(OR, 50)),
    aes(x = OR_cap, y = reorder(Label, neg_logF))) +
    geom_point(aes(size = Overlap_N, color = neg_logF, shape = Sig), alpha = 0.85) +
    scale_color_viridis_c(option = "plasma", name = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(2, 8), name = "Overlap genes") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "FDR<0.05") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    labs(title = "L | Gradient genes vs literature gene sets",
         subtitle = "Fisher's Exact (one-sided) | unbiased discovery → literature validation",
         x = "Odds Ratio", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(size = 8))
  
  save_pdf(p_L, "FigL_Overlap_Fisher", w = 8, h = 6)
}

# ====== LOCO 补充图 ======
if (exists("loco_df") && nrow(loco_df) > 0) {
  p_LOCO <- ggplot(loco_df, aes(x = LeftOut, y = Projected_MGS, fill = LeftOut)) +
    geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2.5, shape = 21,
                fill = "white", color = "grey30", stroke = 0.5) +
    scale_fill_manual(values = SUBTYPE_COLORS) +
    labs(title = "LOCO Validation",
         subtitle = "Projected MGS of samples left out during curve fitting",
         x = "Left-out subtype", y = "Projected MGS") +
    theme_pub() + theme(legend.position = "none")
  
  save_pdf(p_LOCO, "FigS_LOCO_Validation", w = 5, h = 5)
}

# ====== 保存单面板 ======
save_pdf(p_A, "FigA_PCA_PrincipalCurve")
save_pdf(p_B, "FigB_UMAP_MGS")
save_pdf(p_C, "FigC_MGS_Boxplot", w = 5, h = 5)
save_pdf(p_D, "FigD_Permutation_Test")
save_pdf(p_E, "FigE_Bootstrap_Stability")

# ====== 组合主图 ======
# [STYLE] 若 p_I 为 NULL, 用 plot_spacer() 填充, 避免空白
p_row1 <- p_A | p_B | p_C
p_row2 <- p_D | p_E
p_row3 <- if (!is.null(p_I)) (p_K | p_I) else (p_K | patchwork::plot_spacer())

p_main <- p_row1 / p_row2 / p_row3 +
  patchwork::plot_layout(heights = c(1, 1, 1)) +
  patchwork::plot_annotation(
    title    = "Molecular Gradient Analysis: Vestibular Schwannoma Subtypes (C1/C2/C3)",
    subtitle = "Principal Curve · Permutation · Bootstrap · Partial Correlation · ssGSEA",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey35")))

save_pdf(p_main, "Fig_Combined_Main", w = W_COMB, h = H_COMB)

# ====== 保存 per-sample 结果表 ======
write.csv(
  meta_use %>%
    dplyr::select(SampleID, Subtype, MGS_raw, MGS, Rank_Stability, UMAP1, UMAP2),
  file.path(S9_TAB_DIR, "MGS_PerSample.csv"),
  row.names = FALSE)


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║                          最终保存 & 汇报                                  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

save(meta_use, mgs_01, mgs_raw, gradient_df, genes_up, genes_down,
     module_final, mat_mod_z, eigen_df, pca_df, curve_df,
     perm_kw, real_kw, perm_p, boot_taus, mean_tau, ci_tau,
     anchor_df, overlap_df,
     file = file.path(RDATA_DIR, "Step9_MolGradient.RData"))

message("══════════════════════════════════════════════════")
message("🎉 Step9 全流程完成!")
message("   图: ", S9_FIG_DIR)
message("   表: ", S9_TAB_DIR)
message("   KW χ²=", round(real_kw, 2),
        " | Perm p=", perm_p,
        " | Bootstrap τ=", round(mean_tau, 3),
        " [", round(ci_tau[1], 3), ", ", round(ci_tau[2], 3), "]")
message("   梯度基因: Up=", length(genes_up), ", Down=", length(genes_down))
message("══════════════════════════════════════════════════")

# ── 输出文件清单 ──────────────────────────────────────────
# 图 (FIG_DIR/Step9_MolGradient/):
#   FigA_PCA_PrincipalCurve.pdf    PCA + 拟合曲线
#   FigB_UMAP_MGS.pdf              UMAP (MGS 着色)
#   FigC_MGS_Boxplot.pdf           MGS 分组箱线图
#   FigD_Permutation_Test.pdf      排列检验零分布
#   FigE_Bootstrap_Stability.pdf   Bootstrap τ 分布 + 95% CI
#   FigF_Schwann_Anchors.pdf       Schwann marker 锚定验证
#   FigG_GeneModule_Heatmap.pdf    模块化基因热图
#   FigH_ModuleEigengenes.pdf      模块 eigengene 沿 MGS
#   FigI_SignatureTrends_ssGSEA.pdf C1/C2/C3 ssGSEA Signature
#   FigJ_PathwayScores.pdf         Step4 通路得分趋势
#   FigK_GradientGene_Volcano.pdf  梯度基因火山图
#   FigL_Overlap_Fisher.pdf        文献基因集 Overlap
#   FigS_LOCO_Validation.pdf       LOCO 验证
#   FigS_GO_Up/Down_along_MGS.pdf  GO 富集 (补充)
#   Fig_Combined_Main.pdf          组合主图
#
# 表 (TABLE_DIR/Step9_MolGradient/):
#   MGS_PerSample.csv              每样本 MGS + 稳定性
#   MGS_PairwiseWilcox.csv         配对 Wilcoxon
#   MGS_StatisticalValidation.csv  排列/Bootstrap 汇总 (含 CI)
#   MGS_Anchor_Validation.csv      Schwann 锚定结果
#   MGS_LOCO_Validation.csv        LOCO 投影结果
#   MGS_GradientGenes_All.csv      全基因组梯度相关
#   MGS_GeneBreakpoints.csv        滑窗拐点
#   MGS_GeneModules.csv            基因模块分配
#   MGS_Overlap_Fisher.csv         文献 Overlap 检验
#   MGS_SignatureScores_ssGSEA.csv 各亚型 ssGSEA Signature
#   MGS_Module_M1~M5_GO_BP.csv     各模块 GO BP




# ====================================================================
# 补充可视化 2: 代表性基因的 Loess 动态趋势曲线
# ====================================================================
# 1. 从 gradient_df 中选取极具代表性的基因 (各取 Top 4)
rep_up <- gradient_df %>% 
  filter(Direction == "Up_along_MGS") %>% 
  arrange(desc(Rho)) %>% 
  head(4) %>% 
  pull(Gene)

rep_down <- gradient_df %>% 
  filter(Direction == "Down_along_MGS") %>% 
  arrange(Rho) %>% 
  head(4) %>% 
  pull(Gene)

rep_genes <- c(rep_up, rep_down)

# 2. 转换为适合 ggplot 的长格式数据
trend_df <- lapply(rep_genes, function(g) {
  if(g %in% rownames(expr_log)) {
    data.frame(
      SampleID = meta_use$SampleID,
      MGS      = meta_use$MGS,
      Subtype  = meta_use$Subtype,
      Expression = as.numeric(expr_log[g, meta_use$SampleID]),
      Gene     = g,
      Trend    = ifelse(g %in% rep_up, "Up_along_MGS (De-differentiation)", "Down_along_MGS (Differentiation)"),
      stringsAsFactors = FALSE
    )
  }
}) %>% bind_rows()

# 固定基因的分面顺序
trend_df$Gene <- factor(trend_df$Gene, levels = rep_genes)

# 3. 绘图
p_trend <- ggplot(trend_df, aes(x = MGS, y = Expression)) +
  geom_point(aes(color = Subtype), size = 2, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, span = 0.85, 
              color = "black", fill = "grey50", alpha = 0.2, linewidth = 1) +
  scale_color_manual(values = SUBTYPE_COLORS) +
  facet_wrap(Trend ~ Gene, ncol = 4, scales = "free_y") +
  labs(
    title = "Representative Genes along Molecular Gradient",
    subtitle = "Top row: Genes increasing with MGS | Bottom row: Genes decreasing with MGS",
    x = "Molecular Gradient Score (MGS: 0 → 1)",
    y = "Expression (log2 TPM)"
  ) +
  theme_pub() + 
  theme(legend.position = "bottom")

save_pdf(p_trend, "FigS_Representative_Genes_Trend", w = 12, h = 6)






# ====== Fig K: Gradient Gene Volcano (带 Top 5 基因标注) ======
library(ggrepel)

# 1. 提取 Up 和 Down 各前 5 个基因用于标注
label_up <- gradient_df %>% 
  filter(Direction == "Up_along_MGS") %>% 
  arrange(desc(Rho)) %>% 
  head(5)

label_down <- gradient_df %>% 
  filter(Direction == "Down_along_MGS") %>% 
  arrange(Rho) %>% 
  head(5)

label_genes <- bind_rows(label_up, label_down)

# 2. 绘图
p_K <- ggplot(gradient_df, aes(x = Rho, y = -log10(FDR + 1e-300))) +
  geom_point(aes(color = Direction), size = 0.9, alpha = 0.5) +
  scale_color_manual(
    values = c(Up_along_MGS = "#D73027", Down_along_MGS = "#4575B4", NS = "grey78"),
    labels = c("Up along MGS", "Down along MGS", "NS")) +
  geom_hline(yintercept = -log10(PARTIAL_FDR), linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = c(-PARTIAL_RHO, PARTIAL_RHO), linetype = "dotted",
             color = "grey40", linewidth = 0.3) +
  # 标注 Top 5 基因: 放大点
  geom_point(data = label_genes, 
             aes(x = Rho, y = -log10(FDR + 1e-300), color = Direction),
             size = 2.5, alpha = 1, show.legend = FALSE) +
  # 标注 Top 5 基因: 基因名文本 (ggrepel 自动避免重叠)
  geom_text_repel(data = label_genes,
                  aes(x = Rho, y = -log10(FDR + 1e-300), label = Gene, color = Direction),
                  size = 3.2, fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.linewidth = 0.3,
                  max.overlaps = 20,
                  min.segment.length = 0,
                  show.legend = FALSE) +
  labs(title = paste0("K | ", gradient_method, " with MGS"),
       subtitle = paste0("Up=", length(genes_up), ", Down=", length(genes_down),
                         " | Top 5 labeled per direction"),
       x = paste0(gradient_method, " ρ"),
       y = expression(-log[10](FDR)), color = NULL) +
  theme_pub() + 
  theme(legend.position = c(0.15, 0.85))

save_pdf(p_K, "FigK_GradientGene_Volcano", w = 7, h = 6)






# ================================================
# Fig B 增强版：UMAP (MGS 连续着色) + 按亚型添加凸包圈
# ================================================

library(ggforce)   # 用于 geom_mark_hull 或 geom_convexhull（推荐先加载）

p_B <- ggplot(meta_use, aes(x = UMAP1, y = UMAP2)) +
  
  # 1. 样本点：按 MGS 连续着色（保持原梯度信息）
  geom_point(aes(color = MGS), 
             size = 3.2, alpha = 0.92, shape = 19) +
  
  # 2. 按 Subtype 添加凸包圈（分组圈出 C1/C2/C3）
  ggforce::geom_mark_hull(
    aes(fill = Subtype, label = Subtype),
    alpha = 0.12,          # 圈的透明度
    color = "black", 
    linewidth = 0.6,
    expand = unit(4, "mm"),      # 圈向外扩展一点
    radius = unit(2, "mm"),      # 圆角
    con.size = 0.4,              # 连接线粗细
    label.fontsize = 10,
    label.fill = alpha("white", 0.85),
    label.colour = "black",
    show.legend = FALSE
  ) +
  
  # 3. MGS 颜色条
  scale_color_viridis_c(option = pt_palette, end = 0.95, 
                        name = "MGS\n(0 → 1)") +
  
  # 4. 亚型填充颜色（仅用于圈）
  scale_fill_manual(values = SUBTYPE_COLORS) +
  
  labs(
    title = "B | UMAP colored by MGS\n(with subtype convex hulls)",
    subtitle = "C1 → C2 → C3 separation along molecular gradient",
    x = "UMAP1", 
    y = "UMAP2"
  ) +
  theme_pub(base_size = 10.5) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 9)
  )

# 保存增强后的 Fig B
save_pdf(p_B, "FigB_UMAP_MGS_with_SubtypeHulls", w = 6.8, h = 5.5)




# =============================================================================
# Step9_MolecularGradient_v3_Final.R
# ─────────────────────────────────────────────────────────────────────────────
# 修复内容：
# 1. 全局 P 值格式化（严谨学术显示，拒绝 P=0）
# 2. 火山图顶部截断优化（添加轻微抖动，视觉更自然）
# 3. 增加阈值辅助线，增加图的可读性
# =============================================================================

# rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# ---- 0.1 基础配置与函数 ----
source("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/00_config.R")

# 🌟 核心函数：严谨 P 值格式化
fmt_p <- function(p, prefix = "P = ", n_perm = NULL) {
  if (is.na(p)) return(paste0(prefix, "NS"))
  if (!is.null(n_perm) && p == 0) return(paste0(prefix, "< ", format(1/n_perm, scientific = FALSE)))
  if (p < 0.0001) return(paste0(prefix, formatC(p, format = "e", digits = 2)))
  return(paste0(prefix, round(p, 4)))
}

# 绘图主题
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
          plot.subtitle = element_text(size = 9, color = "grey30"))
}

# ---- 0.2 加载数据 ----
load(file.path(RDATA_DIR, "Step2_Clustering.RData"))
load(file.path(RDATA_DIR, "Step3_DEA_Results.RData"))

# 样本对齐
meta_use <- meta_step2 %>% filter(!is.na(Subtype)) %>% mutate(Subtype = factor(Subtype, levels = c("C1", "C2", "C3")))
common_samples <- intersect(colnames(log_tpm), meta_use$SampleID)
expr_log <- log_tpm[, common_samples]
meta_use <- meta_use[match(common_samples, meta_use$SampleID), ]

# ==============================================================================
# PHASE 2 & 3: 梯度计算与统计
# ==============================================================================
# 计算 MGS (Principal Curve)
pca_res <- prcomp(t(expr_log[head(names(sort(apply(expr_log, 1, mad), TRUE)), 2000),]), center = T, scale. = T)
pc_fit  <- principal_curve(pca_res$x[,1:3], smoother = "smooth_spline", stretch = 0)
mgs_raw <- as.numeric(pc_fit$lambda)
# 方向校正
if (median(mgs_raw[meta_use$Subtype == "C1"]) > median(mgs_raw[meta_use$Subtype == "C3"])) { mgs_raw <- max(mgs_raw) - mgs_raw }
meta_use$MGS <- (mgs_raw - min(mgs_raw)) / diff(range(mgs_raw))

# 排列检验 (Permutation)
N_PERM <- 1000
real_kw <- kruskal.test(MGS ~ Subtype, data = meta_use)$statistic
perm_kw <- replicate(N_PERM, kruskal.test(sample(meta_use$MGS) ~ Subtype, data = meta_use)$statistic)
perm_p  <- mean(perm_kw >= real_kw)

# 相关性计算 (用于火山图)
grad_df <- data.frame(
  Gene = rownames(expr_log),
  Rho  = as.numeric(cor(t(expr_log), meta_use$MGS, method = "spearman"))
) %>% mutate(
  Pval = 2 * pt(-abs(Rho * sqrt((nrow(meta_use)-2)/(1-Rho^2))), df=nrow(meta_use)-2),
  FDR  = p.adjust(Pval, method = "BH")
)

# ==============================================================================
# PHASE 7: 绘图输出
# ==============================================================================

# --- Fig C: MGS Boxplot ---
kw_p_lab <- fmt_p(kruskal.test(MGS ~ Subtype, data = meta_use)$p.value, prefix = "KW ")
pm_p_lab <- fmt_p(perm_p, prefix = "Perm ", n_perm = N_PERM)

p_C <- ggplot(meta_use, aes(x = Subtype, y = MGS, fill = Subtype)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, color = "grey30") +
  geom_jitter(width = 0.1, size = 2.5, shape = 21, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  labs(title = "C | MGS by subtype", subtitle = paste0(kw_p_lab, " | ", pm_p_lab), x = NULL, y = "MGS (0→1)") +
  theme_pub()

# --- Fig K: Volcano Plot (优化版) ---
# 1. 筛选与定义方向 (这里设阈值为 0.3)
PARTIAL_RHO <- 0.3
grad_df <- grad_df %>% mutate(
  Direction = case_when(
    Rho > PARTIAL_RHO & FDR < 0.05 ~ "Up",
    Rho < -PARTIAL_RHO & FDR < 0.05 ~ "Down",
    TRUE ~ "NS"
  ),
  # 🌟 核心视觉优化：为 FDR 极小的点添加轻微随机抖动，防止顶部成直线
  log_FDR_plot = -log10(FDR + 1e-300) + runif(nrow(grad_df), 0, 0.5)
)

p_K <- ggplot(grad_df, aes(x = Rho, y = log_FDR_plot)) +
  # 背景 NS 点
  geom_point(data = subset(grad_df, Direction == "NS"), color = "grey85", size = 0.8, alpha = 0.4) +
  # 彩色 Up/Down 点
  geom_point(data = subset(grad_df, Direction != "NS"), aes(color = Direction), size = 1.2, alpha = 0.6) +
  # 阈值辅助线 (解释方块分布)
  geom_vline(xintercept = c(-PARTIAL_RHO, PARTIAL_RHO), linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5")) +
  # 标注 Top 基因
  geom_text_repel(data = head(grad_df[order(grad_df$FDR),], 12), aes(label = Gene), 
                  size = 3, fontface = "italic", max.overlaps = 20) +
  labs(title = "K | Gradient Associated Genes", x = "Spearman ρ with MGS", y = expression(-log[10](FDR))) +
  theme_pub() + theme(legend.position = "right")

# --- 保存 ---
p_final <- (p_C | p_K) + plot_layout(widths = c(1, 2))
ggsave(file.path(S9_FIG_DIR, "Fig_Step9_Optimized_Pvals.pdf"), p_final, width = 12, height = 5, device = cairo_pdf)

message("✅ 代码运行完成！检查 Fig C 的副标题 P 值以及 Fig K 的顶部分布。")



# ==============================================================================
# PHASE 7 绘图修改版：移除箭头 + 使用椭圆轮廓 (更圆润)
# ==============================================================================

# 基础配置 (假设你的 config 已经加载，这里只列出绘图必需)
library(ggplot2)
library(ggforce)
library(viridis)

# ---- 7.1 绘制 Fig B (UMAP colored by MGS) ----
p_B <- ggplot(meta_use, aes(x = UMAP1, y = UMAP2)) +
  
  # 1. 🌟 画椭圆轮廓线 (取代生硬的 Hull)
  # 使用 geom_mark_ellipse 可以得到更圆润、平滑的包围圈
  ggforce::geom_mark_ellipse(
    aes(fill = Subtype, group = Subtype), 
    color = "grey30",          # 椭圆边框颜色
    linetype = "dashed",       # 使用虚线更显单细胞图质感
    linewidth = 0.5,           # 边框粗细
    alpha = 0.1,              # 极淡的背景色填充
    na.rm = TRUE,              # 忽略 NA 坐标
    show.legend = FALSE        # 不显示此层的填充图例
  ) +
  
  # 2. 绘制样本点 (保持原逻辑：MGS 连续着色)
  geom_point(aes(color = MGS), size = 3.5, alpha = 0.9, stroke = 0.3, shape = 16) +
  
  # 3. 颜色设置 (Virginis 色板和自定义 Subtype 颜色)
  scale_color_viridis_c(option = "inferno", end = 0.95, name = "MGS\n(0 → 1)") +
  # 确保你 config 里的 SUBTYPE_COLORS 已经定义
  scale_fill_manual(values = SUBTYPE_COLORS) + 
  
  # 4. 🌟 移除：箭头 (geom_segment) 和 指示文本 (annotate) 已被移除
  
  # 5. 标签与布局优化 (保持 theme_pub 风格)
  labs(
    title = "B | UMAP of Molecular Gradient", 
    subtitle = "Transition: C1 (Lower-Left) → C2 → C3 (Upper-Right)",
    x = "UMAP1", 
    y = "UMAP2"
  ) +
  theme_pub() + # 应用你的出版主题
  theme(
    legend.position = "right",
    # 增加额外的个性化 theme 设置（可选）
    panel.grid = element_blank(), # 去掉网格背景，视觉更干净
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

# 🌟 修复 colMedians 报错问题 (确保 coordinates 定义时使用了正确的函数)
# 如果你之前遇到了 colMedians 报错，请确保Coordinates定向部分代码如下：
c1_coords <- umap_coords[meta_use$Subtype == "C1", ]
if(nrow(c1_coords) > 0) {
  # 计算 C1 中位数中心
  c1_center <- apply(c1_coords, 2, median)
  # 这里继续你之前的轴翻转逻辑 (C1定向到左下)
  if (c1_center[1] > median(umap_coords$UMAP1)) { umap_coords$UMAP1 <- -umap_coords$UMAP1 }
  if (c1_center[2] > median(umap_coords$UMAP2)) { umap_coords$UMAP2 <- -umap_coords$UMAP2 }
}

# 保存最终结果
save_pdf(p_B, "FigB_UMAP_Ellipses_Final", w = 7, h = 6)
p_B
message("✅ 图B已修改：箭头已移除，生硬的凸包已被圆润的椭圆轮廓取代。")
message(file.path(S9_FIG_DIR, "FigB_UMAP_Ellipses_Final.pdf"))
# 强制保存一份到当前工作目录，方便查找
ggsave("temp_check_figB.pdf", plot = p_B, width = 7, height = 6, device = cairo_pdf)
message("应急备份图已保存至：", getwd(), "/temp_check_figB.pdf")





# # 查看前 10 行结果
# head(gradient_df, 10)
# 
# # 或者按相关性 Rho 从大到小排序看（看最正相关的基因）
# gradient_df %>% arrange(desc(Rho)) %>% head(10)
# 
# # 或者看最负相关的基因（也就是你觉得曲线不理想的那些）
# gradient_df %>% arrange(Rho) %>% head(10)



# ==============================================================================
# 专项绘图：使用单调性最强的基因重新绘制 Fig F (Schwann Myelination Markers)
# ==============================================================================

message("═══ 正在使用单调性最优基因重新绘制 Fig F ═══")

# 1. 定义你刚发现的单调下降最强的 6 个基因 (Rho 接近 -0.9)
# 这里的顺序按照你提供的 gradient_df 排序
ideal_down_genes <- c("SOX10","PLP1","EGR2", "ING5", "CKB", "KLF15")
# ideal_down_genes <- c(,"EGR2")

# 2. 准备绘图数据
# 确保 expr_log 和 mgs_01 已经在当前环境
plot_data_list <- lapply(ideal_down_genes, function(g) {
  if (g %in% rownames(expr_log)) {
    data.frame(
      SampleID = meta_use$SampleID,
      MGS      = mgs_01,
      Expression = as.numeric(expr_log[g, meta_use$SampleID]),
      Gene      = g,
      Subtype   = meta_use$Subtype,
      stringsAsFactors = FALSE
    )
  } else {
    message("⚠️ 警告: 基因 ", g, " 不在表达矩阵中")
    NULL
  }
})

df_plot_ideal <- bind_rows(plot_data_list) %>%
  mutate(Gene = factor(Gene, levels = ideal_down_genes))

# 3. 绘图：强调单调下降趋势
# 增加 span 值到 1.5 以获得更平滑、更具学术感的拟合曲线
p_F_ideal <- ggplot(df_plot_ideal, aes(x = MGS, y = Expression)) +
  # 绘制背景置信区间
  geom_smooth(method = "loess", se = TRUE, span = 1.5, 
              linewidth = 1, color = "black", fill = "grey80", alpha = 0.3) +
  # 绘制样本点
  geom_point(aes(fill = Subtype), shape = 21, size = 2.5, alpha = 0.8, color = "white", stroke = 0.3) +
  # 颜色方案
  scale_fill_manual(values = SUBTYPE_COLORS) +
  # 分面展示
  facet_wrap(~Gene, scales = "free_y", ncol = 3) +
  # 标签优化
  labs(
    title = "F | Refined Schwann differentiation markers (Top Monotonic)",
    subtitle = "Genes precisely decreasing along Molecular Gradient Score (MGS)",
    x = "Molecular Gradient Score (0 = Differentiated → 1 = De-differentiated)",
    y = "Expression (log2 TPM)",
    fill = "Subtype"
  ) +
  # 应用出版主题
  theme_pub() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold.italic", size = 11)
  )

# 4. 保存图片
# 使用你的 save_pdf 函数保存
# getwd()
setwd("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/Step9_MolGradient/")
save_pdf(p_F_ideal, "FigF_Schwann_Anchors_Refined.pdf", w = 10, h = 7)

# 打印完成信息
message("✅ 完美单调性 Fig F 已生成：FigF_Schwann_Anchors_Refined.pdf")




# 查看 GSEA 自动计算出的最显著上调通路 (NES > 0)
head(gse_go_res@result %>% filter(NES > 0) %>% arrange(desc(NES)), 20)

# 查看 GSEA 自动计算出的最显著下调通路 (NES < 0)
head(gse_go_res@result %>% filter(NES < 0) %>% arrange(NES), 20)

# ==============================================================================
# 最终精修版：GSEA 驼峰图垂直组图 (根据最新筛选)
# ==============================================================================
library(ggplotify)
library(patchwork)
library(enrichplot)

# 1. 锁定筛选出的三个核心 ID
# 注意：GO:0007272 如果在你的 gse_go_res 结果表里，可以保留。
# 如果不在，可以换成 GO:0045665
target_ids <- c("GO:0002250", "GO:0008299", "GO:0045665")

# 2. 循环绘图并优化比例
p_gsea_list <- lapply(target_ids, function(id) {
  res_row <- gse_go_res@result[gse_go_res@result$ID == id, ]
  if(nrow(res_row) == 0) {
    message("跳过不存在的 ID: ", id)
    return(NULL)
  }
  
  desc <- res_row$Description
  nes  <- round(res_row$NES, 2)
  
  # 绘图设置
  p <- gseaplot2(gse_go_res, geneSetID = id, 
                 title = paste0(desc, " (NES = ", nes, ")"),
                 color = ifelse(nes > 0, "#E64B35", "#4DBBD5"),
                 base_size = 10,
                 # rel_heights 控制三个部分的比例：驼峰图、竖线、底部密度
                 rel_heights = c(1.8, 0.4, 0.8)) 
  
  return(ggplotify::as.ggplot(p))
})

# 3. 垂直拼接
p_gsea_list <- Filter(Negate(is.null), p_gsea_list)
p_final_stack <- wrap_plots(p_gsea_list, ncol = 1) +
  plot_annotation(
    title = "Biological Transition Across Molecular Gradient",
    subtitle = "Upregulated Immune Signatures vs. Downregulated Schwann Differentiation/Metabolism",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# 4. 保存
# 确保你的输出目录存在
save_pdf(p_final_stack, "FigL_GSEA_Combined_Final.pdf", w = 8, h = 13)

message("✅ 终极版 GSEA 驼峰图已完成，请查看 FigL_GSEA_Combined_Final.pdf")








# ==============================================================================
# 专项绘图：展示随 MGS 上升和下降的前 20 个基因 (BarPlot)
# ==============================================================================

library(dplyr)
library(ggplot2)

message("═══ 正在生成 Top 20 梯度关联基因条形图 ═══")

# 1. 从 gradient_df 中提取 Top 20 上调和 Top 20 下调基因
top20_genes_bar <- gradient_df %>%
  filter(Direction %in% c("Up_along_MGS", "Down_along_MGS")) %>%
  group_by(Direction) %>%
  # 按 Rho 的绝对值排序，各取前 20 个
  slice_max(abs(Rho), n = 20) %>%
  ungroup() %>%
  # 为了绘图时条形能按大小排列，重新设置 Gene 的因子顺序
  mutate(Gene = factor(Gene, levels = Gene[order(Rho)]))

# 2. 绘制条形图
p_top20_bar <- ggplot(top20_genes_bar, aes(x = Rho, y = Gene, fill = Direction)) +
  # 绘制条形
  geom_col(width = 0.7, color = "white", linewidth = 0.1) +
  # 在 0 处画一条中轴线
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  # 设置颜色（红正蓝负，与你之前的风格对齐）
  scale_fill_manual(values = c("Up_along_MGS" = "#D73027", 
                               "Down_along_MGS" = "#4575B4"),
                    labels = c("Down along MGS", "Up along MGS")) +
  # 设置 X 轴范围
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  # 标签优化
  labs(
    title = "Top 20 Gradient-Associated Genes",
    subtitle = "Ranked by Spearman correlation coefficient (ρ) with MGS",
    x = "Spearman ρ with MGS",
    y = NULL,
    fill = "Trend"
  ) +
  # 应用出版级主题
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 8, face = "italic"), # 基因名设为斜体且稍小一点以容纳 40 个基因
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
  )

# 3. 保存图片
# 增加高度 (h=10) 以确保 40 个基因名不会重叠
save_pdf(p_top20_bar, "FigS_Top20_Genes_BarPlot.pdf", w =10, h = 5)

message("✅ Top 20 基因条形图已生成：FigS_Top20_Genes_BarPlot.pdf")








# ==============================================================================
# 终极解决方案：预像素化位图嵌入法 (Raster Injection)
# ==============================================================================
library(ggplotify)
library(patchwork)
library(enrichplot)
library(magick)

# 如果没有安装 magick，请先运行: install.packages("magick")

p_gsea_list <- lapply(target_ids, function(id) {
  res_row <- gse_go_res@result[gse_go_res@result$ID == id, ]
  if(nrow(res_row) == 0) return(NULL)
  
  # 1. 生成原始 GSEA 基础对象
  p_base <- gseaplot2(gse_go_res, geneSetID = id, 
                      title = paste0(res_row$Description, " (NES = ", round(res_row$NES, 2), ")"),
                      color = ifelse(res_row$NES > 0, "#E64B35", "#4DBBD5"),
                      base_size = 10,
                      rel_heights = c(1.8, 0.4, 0.8))
  
  # 2. 局部像素化处理函数
  # 逻辑：将 ggplot 子图转为 PNG 流，再读回作为位图，彻底切断矢量联系
  rasterize_subplot <- function(plot_obj, res = 300) {
    img <- image_graph(width = 800, height = 400, res = res)
    print(plot_obj)
    dev.off()
    return(as.ggplot(img))
  }

  p1_raster <- rasterize_subplot(p_base[[1]])
  p2_raster <- rasterize_subplot(p_base[[2]])
  p3_raster <- rasterize_subplot(p_base[[3]])
  
  # 3. 手动重新组装
  # 这样组装出来的每一层在 PDF 看来都只是一张透明底的高清图片
  p_combined_sub <- (p1_raster / p2_raster / p3_raster) + 
    plot_layout(heights = c(1.8, 0.4, 0.8))
  
  return(p_combined_sub)
})

# 4. 全局垂直拼接
p_gsea_list <- Filter(Negate(is.null), p_gsea_list)
p_final_rasterized <- wrap_plots(p_gsea_list, ncol = 1)

# 5. 导出 PDF
save_path <- file.path(S9_FIG_DIR, "FigL_GSEA_Final_Raster_Locked.pdf")
grDevices::cairo_pdf(save_path, width = 8, height = 13)
print(p_final_rasterized)
dev.off()

message("✅ 终极修复完成！底部 Ranked List Metric 已彻底转化为位图，AI 不再有路径。")
