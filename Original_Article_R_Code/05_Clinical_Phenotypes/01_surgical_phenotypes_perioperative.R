# ==========================================
# 1. 加载必要的R包
# ==========================================
library(tidyverse)
library(ggpubr)
library(readxl) 

# ==========================================
# 2. 设置路径与颜色
# ==========================================
input_file <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/cl_bulk_with_Subtype.xlsx"
output_dir <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/Clinical_Subtype/"

# 确保输出文件夹存在
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- paste0(output_dir, "Significant_Indicators_ViolinPlot_8Panels.pdf")

SUBTYPE_COLORS <- c(
  "C1" = "#00A087",   # Proliferative (绿色系)
  "C2" = "#4DBBD5",   # Mesenchymal (蓝色系)
  "C3" = "#E64B35"    # Immune (红色系)
)

# 原始数据的中文列名（从 Excel 提取的基础数据）
target_vars_cn <- c(
  "自然杀伤细胞", 
  "白蛋白", 
  "Age", 
  "平均血小板体积", 
  "总T淋巴细胞", 
  "估算肾小球滤过率(CKD-EPI)", 
  "白/球比"
)

# 基础数据的英文重命名映射
target_vars_en_mapping <- c(
  "NK Cells", 
  "Albumin", 
  "Age", 
  "MPV", 
  "Total T Cells", 
  "eGFR", 
  "A/G Ratio"
)

# 最终要在图上展示的 8 个指标（包含计算得出的新指标）
final_plot_vars <- c(
  "NK Cells", 
  "Albumin", 
  "Age", 
  "MPV", 
  "Total T Cells", 
  "eGFR", 
  "A/G Ratio",
  "NK/Total T Ratio"  # 新增的计算指标
)

# 需要进行两两比较的分组对
my_comparisons <- list( c("C1", "C2"), c("C1", "C3"), c("C2", "C3") )

# ==========================================
# 3. 读取、清洗并计算新指标
# ==========================================
df <- read_excel(input_file)

if(any(grepl("--|×10\\^", as.character(df[1, ])))) {
  df <- df[-1, ]
}

clean_numeric <- function(x) {
  x <- as.character(x)
  x <- gsub("<|>", "", x)
  x <- trimws(x)
  as.numeric(x)
}

plot_data <- df %>%
  select(Subtype, all_of(target_vars_cn)) %>%
  filter(Subtype %in% names(SUBTYPE_COLORS)) %>% 
  mutate(across(all_of(target_vars_cn), clean_numeric)) %>% 
  # 第一步：把原有的7个指标换成英文
  rename_with(~ target_vars_en_mapping, all_of(target_vars_cn)) %>%
  # 第二步：新增计算列 (NK/Total T Ratio)
  mutate(`NK/Total T Ratio` = `NK Cells` / `Total T Cells`) %>%
  # 第三步：将8个指标全部转为长数据格式
  pivot_longer(
    cols = all_of(final_plot_vars),
    names_to = "Indicator",
    values_to = "Value"
  ) %>%
  drop_na(Value)

# 固定这8个指标的排列顺序，确保它出现在右下角的第8个格子
plot_data$Indicator <- factor(plot_data$Indicator, levels = final_plot_vars)

# ==========================================
# 4. 绘制带有两两比较标注的图表 (2行4列完美排版)
# ==========================================
p <- ggplot(plot_data, aes(x = Subtype, y = Value, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.7, color = "gray20") +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  # ncol = 4 会让这 8 个指标正好排成 2 行 4 列
  facet_wrap(~ Indicator, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(y = "Level / Ratio") +
  
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", 
                     label = "p.signif",     
                     hide.ns = TRUE,         
                     step.increase = 0.1,    
                     vjust = 0.5) +
  
  stat_compare_means(method = "kruskal.test", 
                     label = "p.format", 
                     label.y.npc = "bottom", 
                     size = 3.5, 
                     color = "firebrick")

# ==========================================
# 5. 保存高清图片
# ==========================================
ggsave(
  filename = output_file,
  plot = p,
  width = 12,    
  height = 8.5,    
  dpi = 300      
)

cat("分析与绘图完成！带有8个指标的图片已保存至:", output_file, "\n")# ==============================================================================
# 整合分析：外周血指标与分子亚型（已修复 stat_halfeye 兼容性报错）
# ==============================================================================

library(tidyverse)
library(readxl)
library(broom)
library(ggpubr)
library(rstatix)
library(car)

# ------------------------------ 1. 配置 ------------------------------
OUTPUT_DIR <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/自审/Step_Integrated_Blood_Subtype/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

clinical_file <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/cl_bulk_with_Subtype.xlsx"

PAL <- c("C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35")

# ------------------------------ 2. 读取与清洗 ------------------------------
df_raw <- read_excel(clinical_file)
if (any(grepl("--|×10\\^", as.character(df_raw[1,])))) df_raw <- df_raw[-1, ]

clean_numeric <- function(x) as.numeric(trimws(gsub("<|>", "", as.character(x))))

indicators_cn <- c("自然杀伤细胞", "白蛋白", "平均血小板体积",
                   "总T淋巴细胞", "估算肾小球滤过率(CKD-EPI)", "白/球比")
indicators_en <- c("NK_Cells", "Albumin", "MPV", "Total_T", "eGFR", "AG_Ratio")

df <- df_raw %>%
  filter(Subtype %in% c("C1", "C2", "C3")) %>%
  mutate(
    Subtype = factor(Subtype, levels = c("C1", "C2", "C3")),
    Age = clean_numeric(Age),
    across(all_of(indicators_cn), clean_numeric)
  )

# base R 重命名，彻底避免 rename_with 报错
for (i in seq_along(indicators_cn)) {
  idx <- which(colnames(df) == indicators_cn[i])
  if (length(idx) > 0) colnames(df)[idx] <- indicators_en[i]
}

# 肿瘤大小（用于 full model）
size_col <- grep("SizeGrade|肿瘤大小|Size|最大径", colnames(df_raw),
                 value = TRUE, ignore.case = TRUE)[1]

if (!is.na(size_col) && size_col %in% colnames(df_raw)) {
  # 直接在 df 中匹配（按行顺序，两表行数相同）
  df$Max_Dim <- sapply(df_raw[[size_col]][df_raw$Subtype %in% c("C1","C2","C3")], function(s) {
    nums <- suppressWarnings(as.numeric(strsplit(as.character(s), "\\*|×|x|X")[[1]]))
    if (all(is.na(nums))) NA_real_ else max(nums, na.rm = TRUE)
  })
  df$Tumor_Size <- factor(
    case_when(df$Max_Dim > 20 ~ "Large", df$Max_Dim <= 20 ~ "Small", TRUE ~ NA_character_),
    levels = c("Small", "Large")
  )
}

cat("样本量:", nrow(df), "\n"); print(table(df$Subtype))

# ------------------------------ 3. 核心绘图函数（原生 Raincloud）------------------------------
make_raincloud <- function(data, y_var, title_label, subtitle_label,
                           pairwise_df = NULL, overall_p_label = NULL) {
  
  # 计算用于放置 p 值标注的 y 上限
  y_vals <- data[[y_var]]
  y_max  <- max(y_vals, na.rm = TRUE)
  y_range <- diff(range(y_vals, na.rm = TRUE))
  
  p <- ggplot(data, aes(x = Subtype, y = .data[[y_var]], fill = Subtype, color = Subtype)) +
    # 1. 小提琴（模拟云朵部分）
    geom_violin(trim = FALSE, adjust = 1.3, alpha = 0.35, color = NA, width = 0.8) +
    # 2. 箱线图
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white",
                 color = "black", linewidth = 0.5, alpha = 0.9) +
    # 3. 散点（雨滴）
    geom_jitter(width = 0.09, size = 2.0, alpha = 0.75) +
    scale_fill_manual(values = PAL) +
    scale_color_manual(values = PAL) +
    labs(title = title_label,
         subtitle = subtitle_label,
         x = NULL, y = y_var) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
      axis.text.x   = element_text(face = "bold", color = "black", size = 12)
    )
  
  # 两两比较显著性标注
  if (!is.null(pairwise_df)) {
    sig_pairs <- pairwise_df %>% filter(p.adj < 0.05)
    if (nrow(sig_pairs) > 0) {
      comparisons_list <- lapply(seq_len(nrow(sig_pairs)), function(i) {
        c(as.character(sig_pairs$group1[i]), as.character(sig_pairs$group2[i]))
      })
      p <- p + stat_compare_means(
        comparisons = comparisons_list,
        method = "wilcox.test",
        method.args = list(exact = FALSE),
        label = "p.format",
        step.increase = 0.12,
        tip.length = 0.02,
        size = 3.5
      )
    }
  }
  
  # 整体 p 值标注（右上角）
  if (!is.null(overall_p_label)) {
    p <- p + annotate("text",
                      x = 2.5, y = y_max + y_range * 0.12,
                      label = overall_p_label,
                      hjust = 0.5, size = 4, fontface = "italic", color = "grey30")
  }
  
  return(p)
}

# ------------------------------ 4. 主分析循环 ------------------------------
results_summary <- list()

for (var in indicators_en) {
  cat("\n=== 处理指标:", var, "===\n")
  
  df_var <- df %>% drop_na(all_of(c(var, "Age", "Subtype")))
  if (nrow(df_var) < 20) { cat("样本量不足，跳过\n"); next }
  
  # ── 统计检验 ──────────────────────────────────────────────────────────────
  # Unadjusted
  kw       <- kruskal.test(as.formula(paste(var, "~ Subtype")), data = df_var)
  dunn_res <- dunn_test(as.formula(paste(var, "~ Subtype")),
                        data = df_var, p.adjust.method = "BH")
  
  # Age-adjusted
  fit_age   <- lm(as.formula(paste(var, "~ Subtype + Age")), data = df_var)
  anova_age <- car::Anova(fit_age, type = "II")
  
  # 年龄校正后的两两比较（emmeans 思路：用线性模型残差 Wilcoxon）
  df_var$Resid_Age <- residuals(fit_age) + mean(df_var[[var]], na.rm = TRUE)
  dunn_age <- dunn_test(Resid_Age ~ Subtype, data = df_var, p.adjust.method = "BH")
  
  # Full-adjusted (Age + Tumor_Size)
  has_size <- "Tumor_Size" %in% colnames(df_var) && sum(!is.na(df_var$Tumor_Size)) > 5
  fit_full <- NULL; anova_full <- NULL; dunn_full <- NULL
  
  if (has_size) {
    df_full <- df_var %>% drop_na(Tumor_Size)
    if (nrow(df_full) >= 15) {
      fit_full   <- lm(as.formula(paste(var, "~ Subtype + Age + Tumor_Size")), data = df_full)
      anova_full <- car::Anova(fit_full, type = "II")
      df_full$Resid_Full <- residuals(fit_full) + mean(df_full[[var]], na.rm = TRUE)
      dunn_full <- dunn_test(Resid_Full ~ Subtype, data = df_full, p.adjust.method = "BH")
    }
  }
  
  tidy_age  <- tidy(fit_age, conf.int = TRUE) %>% filter(term == "SubtypeC3")
  tidy_full <- if (!is.null(fit_full)) tidy(fit_full, conf.int = TRUE) %>%
    filter(term == "SubtypeC3") else NULL
  
  results_summary[[var]] <- tibble(
    Indicator          = var,
    N                  = nrow(df_var),
    KW_p               = kw$p.value,
    Age_adj_Overall_p  = anova_age["Subtype", "Pr(>F)"],
    Full_adj_Overall_p = if (!is.null(anova_full)) anova_full["Subtype", "Pr(>F)"] else NA,
    C3_vs_C1_Age_beta  = tidy_age$estimate,
    C3_vs_C1_Age_CI_lo = tidy_age$conf.low,
    C3_vs_C1_Age_CI_hi = tidy_age$conf.high,
    C3_vs_C1_Age_p     = tidy_age$p.value,
    C3_vs_C1_Full_beta = if (!is.null(tidy_full)) tidy_full$estimate else NA,
    C3_vs_C1_Full_p    = if (!is.null(tidy_full)) tidy_full$p.value else NA
  )
  
  # ── 可视化 ──────────────────────────────────────────────────────────────────
  fmt_p <- function(p) if (is.na(p)) "" else if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
  
  p1 <- make_raincloud(
    data           = df_var,
    y_var          = var,
    title_label    = "Unadjusted",
    subtitle_label = paste0("KW ", fmt_p(kw$p.value)),
    pairwise_df    = dunn_res,
    overall_p_label = NULL
  )
  
  p2 <- make_raincloud(
    data           = df_var,
    y_var          = "Resid_Age",
    title_label    = "Age-adjusted",
    subtitle_label = paste0("ANCOVA ", fmt_p(anova_age["Subtype", "Pr(>F)"])),
    pairwise_df    = dunn_age,
    overall_p_label = NULL
  ) + labs(y = paste0(var, " (age-adj residuals)"))
  
  plot_list <- list(p1, p2)
  n_cols <- 2
  
  if (!is.null(fit_full) && !is.null(dunn_full)) {
    p3 <- make_raincloud(
      data           = df_full,
      y_var          = "Resid_Full",
      title_label    = "Age+Size-adjusted",
      subtitle_label = paste0("ANCOVA ", fmt_p(anova_full["Subtype", "Pr(>F)"])),
      pairwise_df    = dunn_full,
      overall_p_label = NULL
    ) + labs(y = paste0(var, " (full-adj residuals)"))
    
    plot_list <- list(p1, p2, p3)
    n_cols <- 3
  }
  
  combined <- ggarrange(
    plotlist     = plot_list,
    ncol         = n_cols,
    common.legend = TRUE,
    legend       = "bottom"
  )
  
  out_w <- ifelse(n_cols == 3, 15, 10)
  ggsave(
    file.path(OUTPUT_DIR, paste0("Raincloud_", var, ".pdf")),
    combined, width = out_w, height = 6.5, device = cairo_pdf
  )
  cat("✅", var, "图表已保存\n")
}

# ------------------------------ 5. 输出汇总表格 ------------------------------
final_table <- bind_rows(results_summary)
write.csv(final_table,
          file.path(OUTPUT_DIR, "Integrated_Blood_Subtype_Summary.csv"),
          row.names = FALSE)

cat("\n🎉 全部完成！输出目录：", OUTPUT_DIR, "\n")
print(final_table)












# ==============================================================================
# 顶刊标准多面板组合图：Age + Tumor Size 双校正后的临床指标分布
# (修复版：使用 emmeans 计算校正后的两两比较 P 值，并强制标注在图上)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggpubr)
  library(car)
  library(cowplot)
  library(emmeans) # 🌟 核心新增：用于多变量模型的严格事后检验 (Post-hoc)
})

OUTPUT_DIR <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/自审/Figure5_Final/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

clinical_file <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/cl_bulk_with_Subtype.xlsx"
PAL <- c("C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35")

# --- 数据读取与清洗 (保持不变) ---
df_raw <- read_excel(clinical_file)
if (any(grepl("--|×10\\^", as.character(df_raw[1,])))) df_raw <- df_raw[-1, ]

clean_numeric <- function(x) as.numeric(trimws(gsub("<|>", "", as.character(x))))
indicators_cn <- c("自然杀伤细胞", "平均血小板体积", "总T淋巴细胞", "白蛋白")
indicators_en <- c("NK_Cells", "MPV", "Total_T", "Albumin")

df <- df_raw %>%
  filter(Subtype %in% c("C1", "C2", "C3")) %>%
  mutate(
    Subtype = factor(Subtype, levels = c("C1", "C2", "C3")),
    Age = clean_numeric(Age),
    across(all_of(indicators_cn), clean_numeric)
  )

for (i in seq_along(indicators_cn)) {
  idx <- which(colnames(df) == indicators_cn[i])
  if (length(idx) > 0) colnames(df)[idx] <- indicators_en[i]
}

size_col <- grep("SizeGrade|肿瘤大小|Size|最大径", colnames(df_raw), value = TRUE, ignore.case = TRUE)[1]
if (!is.na(size_col)) {
  df$Max_Dim <- sapply(df_raw[[size_col]][df_raw$Subtype %in% c("C1","C2","C3")], function(s) {
    nums <- suppressWarnings(as.numeric(strsplit(as.character(s), "\\*|×|x|X")[[1]]))
    if (all(is.na(nums))) NA_real_ else max(nums, na.rm = TRUE)
  })
  df$Tumor_Size <- factor(
    case_when(df$Max_Dim > 20 ~ "Large", df$Max_Dim <= 20 ~ "Small", TRUE ~ NA_character_),
    levels = c("Small", "Large")
  )
}

# ------------------------------ 自动化校正与绘图函数 ------------------------------

make_native_adjusted_panel <- function(data, var_name, y_label, title_desc) {
  
  df_tmp <- data %>% drop_na(all_of(c(var_name, "Age", "Tumor_Size", "Subtype")))
  if (nrow(df_tmp) < 15) return(NULL)
  
  # 1. 拟合多变量模型 & 计算残差 (用于画图)
  formula_str <- paste(var_name, "~ Subtype + Age + Tumor_Size")
  fit <- lm(as.formula(formula_str), data = df_tmp)
  df_tmp$Plot_Value <- residuals(fit) + mean(df_tmp[[var_name]], na.rm = TRUE)
  
  # 2. 全局 P 值 (ANCOVA)
  anova_res <- car::Anova(fit, type = "II")
  global_p <- anova_res["Subtype", "Pr(>F)"]
  p_title <- if(global_p < 0.001) "ANCOVA p < 0.001" else sprintf("ANCOVA p = %.4f", global_p)
  
  # 🌟 3. 核心修复：使用 emmeans 进行严格的多变量事后两两比较 🌟
  emm <- emmeans(fit, ~ Subtype)
  pairs_res <- as.data.frame(pairs(emm, adjust = "tukey")) # 使用 Tukey 校正多重比较
  
  # 提取显著的对比，格式化为 ggplot 需要的数据框
  sig_pairs <- pairs_res %>%
    filter(p.value < 0.05) %>% # 只提取显著的
    mutate(
      # 解析 contrasts 列 (例如 "C1 - C2") 提取出 group1 和 group2
      group1 = str_trim(str_split(contrast, "-")[[1]][1]),
      group2 = str_trim(str_split(contrast, "-")[[1]][2]),
      # 生成星号
      p.signif = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  # 4. Y 轴坐标计算
  y_min <- min(df_tmp$Plot_Value, na.rm = TRUE)
  y_max <- max(df_tmp$Plot_Value, na.rm = TRUE)
  y_range <- y_max - y_min
  y_bottom_text <- y_min - (y_range * 0.15)
  
  # 5. 绘图
  p <- ggplot(df_tmp, aes(x = Subtype, y = Plot_Value, fill = Subtype, color = Subtype)) +
    geom_violin(trim = FALSE, adjust = 1.3, alpha = 0.4, color = NA, width = 0.85) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.5, alpha = 0.95) +
    geom_jitter(width = 0.1, size = 1.2, alpha = 0.7) +
    
    geom_text(data = df_tmp %>% count(Subtype), aes(x = Subtype, y = y_bottom_text, label = paste0("n=", n)), 
              inherit.aes = FALSE, size = 3.5, color = "grey40", fontface = "bold") +
    
    scale_fill_manual(values = PAL) +
    scale_color_manual(values = PAL) +
    coord_cartesian(ylim = c(y_bottom_text - (y_range * 0.05), y_max + (y_range * 0.45))) +
    
    labs(
      title = paste0(y_label, "\n", title_desc),
      subtitle = p_title,
      x = NULL, y = y_label
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13, lineheight = 1.2),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 11, margin = margin(b = 10)),
      axis.text.x = element_text(face = "bold", color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(face = "bold"),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
    )
  
  # 🌟 6. 核心修复：把 emmeans 算出来的显著星号手动加到图上 🌟
  if (nrow(sig_pairs) > 0) {
    # 为每条比较线计算 Y 轴高度
    sig_pairs$y.position <- y_max + y_range * seq(0.15, 0.35, length.out = nrow(sig_pairs))
    
    p <- p + stat_pvalue_manual(
      sig_pairs, 
      label = "p.signif", # 显示星号
      tip.length = 0.015,
      bracket.size = 0.5,
      size = 5
    )
  }
  
  return(p)
}

# ------------------------------ 5. 生成四个子图面板 ------------------------------
p_h <- make_native_adjusted_panel(df, "NK_Cells", "NK Cells", "(higher in C3)")
p_i <- make_native_adjusted_panel(df, "MPV", "MPV", "(lower in C3)")
p_j <- make_native_adjusted_panel(df, "Total_T", "Total T Cells", "(higher in C3)")
p_k <- make_native_adjusted_panel(df, "Albumin", "Albumin", "(lower in C3)")

# ------------------------------ 6. 无缝拼接与全局标题 (Cowplot) ------------------------------
combined_plot <- plot_grid(
  p_h, p_i, p_j, p_k, 
  ncol = 4, align = 'vh', labels = c("h", "i", "j", "k"), label_size = 20, label_fontface = "plain"
)

title_wrapper <- ggdraw() + 
  draw_label("Figure 5 h–k: Age & Size Adjusted Laboratory Markers by Subtype", 
             fontface = 'bold', size = 16, vjust = 1) +
  draw_label("Adjusted values (residuals + mean). Pairwise comparisons via emmeans (* p<0.05, ** p<0.01, *** p<0.001)", 
             fontface = 'plain', size = 11, color = "grey30", vjust = 4)

final_output <- plot_grid(title_wrapper, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

save_path <- file.path(OUTPUT_DIR, "Figure_5_h_k_FullyAdjusted_Fixed.pdf")
ggsave(save_path, plot = final_output, width = 10, height = 8, device = cairo_pdf)

message("✅ 完美复现！已使用 emmeans 计算多变量两两比较，并成功标注星号：", save_path)# ==============================================================================
# Step 7-Clin Addon v3: Publication Sankey (Subtype Flow) + Raincloud plots
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. Config
# ------------------------------------------------------------------------------
CONFIG_PATH <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/00_config.R"
if (file.exists(CONFIG_PATH)) {
  source(CONFIG_PATH)
  setwd(PROJECT_DIR)
} else {
  warning("00_config.R not found. Using local defaults.")
  PROJECT_DIR <- getwd()
  FIG_DIR     <- file.path(PROJECT_DIR, "Figures")
  TABLE_DIR   <- file.path(PROJECT_DIR, "Tables")
  RDATA_DIR   <- file.path(PROJECT_DIR, "RData")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(forcats)
  library(scales)
  library(stringr)
  library(ggalluvial)
  library(rlang)
  library(ggnewscale)
})

HAS_GGDIST   <- requireNamespace("ggdist", quietly = TRUE)
HAS_GGSIGNIF <- requireNamespace("ggsignif", quietly = TRUE)

if (HAS_GGDIST)   suppressPackageStartupMessages(library(ggdist))
if (HAS_GGSIGNIF) suppressPackageStartupMessages(library(ggsignif))

OUT_FIG_DIR <- file.path(FIG_DIR, "Step7_Clinical")
OUT_TAB_DIR <- file.path(TABLE_DIR, "Step7_Clinical")
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Theme / save helpers
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

if (!exists("theme_publication")) {
  theme_publication <- function(base_size = 11) {
    theme_classic(base_size = base_size) +
      theme(
        axis.line        = element_line(linewidth = 0.5, colour = "black"),
        axis.text        = element_text(color = "black", size = base_size - 1),
        axis.title       = element_text(color = "black", size = base_size),
        legend.text      = element_text(size = base_size - 1),
        legend.title     = element_text(size = base_size, face = "bold"),
        plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = base_size),
        panel.grid       = element_blank(),
        plot.margin      = margin(8, 10, 8, 8)
      )
  }
}

save_pdf_safe <- function(plot_obj, filename, width = 7, height = 5) {
  out <- file.path(OUT_FIG_DIR, filename)
  ggsave(out, plot = plot_obj, width = width, height = height,
         device = cairo_pdf, bg = "white")
  message("✅ Saved: ", out)
}

format_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# ------------------------------------------------------------------------------
# 2. Robust Column Data Cleaning (Regex Supported)
# ------------------------------------------------------------------------------
resolve_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (length(hit) == 0 || is.na(hit)) {
    if (required) stop("❌ Missing column. Tried: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  hit
}

clean_yes_no <- function(x) {
  x <- trimws(as.character(x))
  case_when(
    grepl("有|Yes|Y|1|丰富|富", x, ignore.case = TRUE) ~ "Yes",
    grepl("无|No|N|0|不丰富|否", x, ignore.case = TRUE) ~ "No",
    TRUE ~ NA_character_
  )
}

clean_texture <- function(x) {
  x <- trimws(as.character(x))
  out <- case_when(
    grepl("软|Soft", x, ignore.case = TRUE) ~ "Soft",
    grepl("中|一般|Medium", x, ignore.case = TRUE) ~ "Medium",
    grepl("硬|Hard", x, ignore.case = TRUE) ~ "Hard",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("Soft", "Medium", "Hard"), ordered = TRUE)
}

clean_koos <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("级", "", x)
  x <- gsub("^KOOS\\s*", "", x, ignore.case = TRUE)
  x <- toupper(x)
  out <- case_when(
    grepl("IV|4", x) ~ "IV",
    grepl("III|3", x) ~ "III",
    grepl("II|2", x) ~ "II",
    grepl("I|1", x) ~ "I",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("I", "II", "III", "IV"), ordered = TRUE)
}

clean_resection <- function(x) {
  x <- trimws(as.character(x))
  out <- case_when(
    grepl("近全|NTR", x, ignore.case = TRUE) ~ "NTR",
    grepl("次全|STR", x, ignore.case = TRUE) ~ "STR",
    grepl("部分|大部分|PR", x, ignore.case = TRUE) ~ "PR",
    grepl("全切|GTR", x, ignore.case = TRUE) ~ "GTR",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("GTR", "NTR", "STR", "PR"))
}

detect_subtype_col <- function(df) {
  subtype_col <- resolve_col(df, c("Subtype", "group", "Group", "Unnamed: 0", "...1"), required = FALSE)
  if (!is.na(subtype_col)) return(subtype_col)
  
  detect_cols <- sapply(df, function(col) {
    vals <- trimws(as.character(col))
    mean(vals %in% c("C1", "C2", "C3"), na.rm = TRUE)
  })
  detect_cols[sapply(df, is.numeric)] <- 0
  
  best_col_idx <- which.max(detect_cols)
  if (length(best_col_idx) == 1 && detect_cols[best_col_idx] >= 0.5) {
    return(colnames(df)[best_col_idx])
  }
  stop("❌ Subtype column not found.")
}

# ------------------------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------------------------
df <- NULL
clinical_rdata <- file.path(RDATA_DIR, "Clinical_with_Subtype.RData")
clinical_xlsx1 <- file.path(PROJECT_DIR, "VS 手术分析.xlsx")
clinical_xlsx2 <- file.path(PROJECT_DIR, "VS 手术分析(1).xlsx")
clinical_xlsx3 <- file.path(PROJECT_DIR, "cl_bulk_with_Subtype.xlsx")

if (is.null(df) && file.exists(clinical_xlsx1)) df <- readxl::read_xlsx(clinical_xlsx1, sheet = 1)
if (is.null(df) && file.exists(clinical_xlsx2)) df <- readxl::read_xlsx(clinical_xlsx2, sheet = 1)
if (is.null(df) && file.exists(clinical_xlsx3)) df <- readxl::read_xlsx(clinical_xlsx3, sheet = 1)
if (is.null(df) && file.exists(clinical_rdata)) {
  load(clinical_rdata)
  if (exists("cl_data_with_subtype")) df <- cl_data_with_subtype
}
if (is.null(df)) stop("❌ No clinical data found.")

colnames(df) <- trimws(stringr::str_squish(colnames(df)))
empty_names <- which(colnames(df) == "" | is.na(colnames(df)))
if (length(empty_names) > 0) colnames(df)[empty_names] <- paste0("X_empty_", seq_along(empty_names))

sample_col    <- resolve_col(df, c("SampleID", "sample_id", "Name"))
subtype_col   <- detect_subtype_col(df)
koos_col      <- resolve_col(df, c("KOOS", "Koos", "koos"))
blood_col     <- resolve_col(df, c("血供丰富", "血供"))
texture_col   <- resolve_col(df, c("肿瘤质地"))
brain_col     <- resolve_col(df, c("脑干面粘连"))
resection_col <- resolve_col(df, c("切除情况", "切除", "切除状态"))

df_plot <- df %>%
  filter(!is.na(.data[[sample_col]]), .data[[sample_col]] != "") %>%
  mutate(
    SampleID = as.character(.data[[sample_col]]),
    Subtype  = factor(trimws(as.character(.data[[subtype_col]])), levels = c("C1", "C2", "C3")),
    KOOS_clean = clean_koos(.data[[koos_col]]),
    BloodSupply = factor(clean_yes_no(.data[[blood_col]]), levels = c("No", "Yes")),
    Texture = clean_texture(.data[[texture_col]]),
    BrainstemAdhesion = factor(clean_yes_no(.data[[brain_col]]), levels = c("No", "Yes")),
    ResectionClean = clean_resection(.data[[resection_col]])
  ) %>%
  filter(!is.na(Subtype))

write.csv(df_plot, file.path(OUT_TAB_DIR, "Clinical_Cleaned_for_ClinicalPlots.csv"), row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# 4. Sankey with Flow tracking Subtype (C1/C2/C3)
# ------------------------------------------------------------------------------
sankey_case <- df_plot %>%
  filter(
    !is.na(KOOS_clean), !is.na(BloodSupply), !is.na(BrainstemAdhesion),
    !is.na(ResectionClean), !is.na(Subtype)
  ) %>%
  mutate(
    KOOS_clean = factor(KOOS_clean, levels = c("I", "II", "III", "IV"), ordered = TRUE),
    BloodSupply = factor(BloodSupply, levels = c("No", "Yes")),
    BrainstemAdhesion = factor(BrainstemAdhesion, levels = c("No", "Yes")),
    ResectionClean = factor(ResectionClean, levels = c("GTR", "NTR", "STR", "PR")),
    Subtype = factor(Subtype, levels = c("C1", "C2", "C3"))
  )

# ================= 顶刊配色方案 (Top-Tier Journal Colors) =================
# 核心数据流（Subtype）：使用高辨识度 NPG Nature 配色
SUBTYPE_COLORS <- c(C1 = "#00A087", C2 = "#4DBBD5", C3 = "#E64B35")

# 节点颜色：使用莫兰迪低饱和度配色，烘托核心流线，不抢眼
# KOOS: 渐变冷灰蓝 (Lancet-inspired)
KOOS_NODE_COLORS <- c("I" = "#E2E9F3", "II" = "#B5C8DF", "III" = "#82A0C2", "IV" = "#51779F")

# 血供 (Blood): 灰绿 (No) vs 柔和灰粉 (Yes)
BLOOD_NODE_COLORS <- c("No" = "#AABCAE", "Yes" = "#C6878F")

# 脑干面粘连 (Adhesion): 灰蓝 (No) vs 灰紫 (Yes)
ADHESION_NODE_COLORS <- c("No" = "#B8C5D6", "Yes" = "#988B9E")

# 切除 (Resection): 纯净灰度，使汇集在此处的 C1/C2/C3 颜色最为鲜明
RESECTION_NODE_COLORS <- c("GTR" = "#EAEAEA", "NTR" = "#D0D0D0", "STR" = "#A6A6A6", "PR" = "#7D7D7D")
# ==============================================================================

make_stratum_df <- function(data, var, x_pos) {
  out <- data %>% count(.data[[var]], name = "Freq")
  colnames(out)[1] <- "stratum"
  out$x <- x_pos
  out
}

stratum_koos      <- make_stratum_df(sankey_case, "KOOS_clean", 1)
stratum_blood     <- make_stratum_df(sankey_case, "BloodSupply", 2)
stratum_brain     <- make_stratum_df(sankey_case, "BrainstemAdhesion", 3)
stratum_resection <- make_stratum_df(sankey_case, "ResectionClean", 4)

p_sankey <- ggplot(
  sankey_case,
  aes(
    axis1 = KOOS_clean,
    axis2 = BloodSupply,
    axis3 = BrainstemAdhesion,
    axis4 = ResectionClean,
    y = 1
  )
) +
  # 1. 数据流层：映射到 Subtype (C1/C2/C3)
  geom_alluvium(aes(fill = Subtype), width = 0.18, alpha = 0.68, knot.pos = 0.35, curve_type = "cubic") +
  scale_fill_manual(values = SUBTYPE_COLORS, name = "Molecular Subtype") +
  
  # 2. KOOS 节点
  ggnewscale::new_scale_fill() +
  geom_stratum(data = stratum_koos, aes(x = x, stratum = stratum, y = Freq, fill = stratum), width = 0.22, color = "grey40", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(values = KOOS_NODE_COLORS, guide = "none") +
  
  # 3. Blood 节点
  ggnewscale::new_scale_fill() +
  geom_stratum(data = stratum_blood, aes(x = x, stratum = stratum, y = Freq, fill = stratum), width = 0.22, color = "grey40", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(values = BLOOD_NODE_COLORS, guide = "none") +
  
  # 4. Brainstem 节点
  ggnewscale::new_scale_fill() +
  geom_stratum(data = stratum_brain, aes(x = x, stratum = stratum, y = Freq, fill = stratum), width = 0.22, color = "grey40", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(values = ADHESION_NODE_COLORS, guide = "none") +
  
  # 5. Resection 节点
  ggnewscale::new_scale_fill() +
  geom_stratum(data = stratum_resection, aes(x = x, stratum = stratum, y = Freq, fill = stratum), width = 0.22, color = "grey40", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(values = RESECTION_NODE_COLORS, guide = "none") +
  
  # 标签绘制
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4.2, family = "sans", fontface = "plain") +
  
  scale_x_discrete(
    limits = c("Pre-op KOOS", "Blood supply", "Brainstem adhesion", "Resection Outcome"),
    expand = c(0.05, 0.05)
  ) +
  labs(x = NULL, y = "Number of Cases") +
  theme_publication(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", size = 12, color = "grey20"),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid = element_blank()
  )

save_pdf_safe(p_sankey, "Fig_Clinical_Sankey_SubtypeFlow.pdf", width = 11.0, height = 7.0)


# ------------------------------------------------------------------------------
# 5. Generic stacked bar 
# ------------------------------------------------------------------------------
plot_stacked_bar <- function(data, var, fill_values, filename, legend_title = NULL) {
  var_sym <- rlang::sym(var)
  
  tab <- data %>%
    filter(!is.na(!!var_sym)) %>%
    count(Subtype, !!var_sym, name = "n") %>%
    group_by(Subtype) %>%
    mutate(
      n_total = sum(n),
      prop = n / n_total,
      label = paste0(round(prop * 100), "%")
    ) %>%
    ungroup()
  
  n_lab <- tab %>% distinct(Subtype, n_total)
  p_global <- tryCatch(fisher.test(table(data$Subtype, data[[var]]))$p.value, error = function(e) NA_real_)
  
  p <- ggplot(tab, aes(x = Subtype, y = prop, fill = !!var_sym)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.35, position = "stack") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 4.6) +
    geom_text(
      data = n_lab, aes(x = Subtype, y = 1.045, label = paste0("n=", n_total)),
      inherit.aes = FALSE, size = 4, color = "grey20"
    ) +
    annotate("text", x = 2, y = 1.10, label = paste0("Fisher p = ", format_p(p_global)), size = 4.2, fontface = "bold") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.12), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = fill_values, name = legend_title %||% var) +
    labs(x = NULL, y = "Proportion") +
    theme_publication(base_size = 12) +
    theme(legend.position = "right", axis.title.x = element_blank())
  
  save_pdf_safe(p, filename, width = 6.2, height = 5.8)
  invisible(tab)
}

texture_cols <- c(Soft = "#7FCDBB", Medium = "#FDBB84", Hard = "#B30000")
brain_cols <- c(No = "#B8C5D6", Yes = "#988B9E")  # Matched the new node colors

tab_texture <- plot_stacked_bar(df_plot, "Texture", texture_cols, "Fig_Texture_StackedBar_by_Subtype.pdf", "Texture")
tab_brain <- plot_stacked_bar(df_plot, "BrainstemAdhesion", brain_cols, "Fig_BrainstemAdhesion_StackedBar_by_Subtype.pdf", "Brainstem adhesion")







# ==============================================================================
# Step 7-Clin Addon v12: Bulletproof Texture & Adhesion (Native Engine)
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. 核心包加载 (去除了不稳定的 ggdist 和 patchwork，改用最稳健的 gridExtra)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
  library(gridExtra) # 🌟 核心修复：使用最底层的拼图包
})

# ------------------------------------------------------------------------------
# 2. 顶刊级配色方案
# ------------------------------------------------------------------------------
PAL_SUBTYPE <- c("C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35")
PAL_TEXTURE <- c("Soft" = "#B2E2E2", "Medium" = "#FEB24C", "Hard" = "#E31A1C")
PAL_ADHESION <- c("No" = "#A6CEE3", "Yes" = "#FB9A99")

# ------------------------------------------------------------------------------
# 3. 稳健的数据读取与全量清洗
# ------------------------------------------------------------------------------
PROJECT_DIR <- getwd()
FIG_DIR     <- file.path(PROJECT_DIR, "Figures", "Step7_Clinical")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

df <- NULL
if (file.exists("VS 手术分析.xlsx")) df <- readxl::read_xlsx("VS 手术分析.xlsx", sheet = 1)
if (is.null(df) && file.exists("cl_bulk_with_Subtype.xlsx")) df <- readxl::read_xlsx("cl_bulk_with_Subtype.xlsx", sheet = 1)
if (is.null(df) && file.exists("RData/Clinical_with_Subtype.RData")) {
  load("RData/Clinical_with_Subtype.RData")
  if (exists("cl_data_with_subtype")) df <- cl_data_with_subtype
}
if (is.null(df)) stop("❌ 未找到临床数据！")

colnames(df) <- trimws(stringr::str_squish(colnames(df)))

resolve_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) stop("Missing column: ", paste(candidates, collapse = ", "))
  return(hit)
}

detect_subtype_col <- function(df) {
  candidates <- c("Subtype", "group", "Group", "Unnamed: 0", "...1")
  hit <- candidates[candidates %in% colnames(df)][1]
  if (!is.na(hit)) return(hit)
  detect_cols <- sapply(df, function(col) mean(trimws(as.character(col)) %in% c("C1", "C2", "C3"), na.rm = TRUE))
  best_idx <- which.max(detect_cols)
  if (length(best_idx) == 1 && detect_cols[best_idx] >= 0.3) return(colnames(df)[best_idx])
  stop("❌ 无法识别 Subtype 列！")
}

clean_yes_no <- function(x) case_when(grepl("有|Yes|Y|1|丰富|富", x, ignore.case=T) ~ "Yes", TRUE ~ "No")
clean_texture <- function(x) {
  x <- trimws(as.character(x))
  out <- case_when(
    grepl("软|Soft", x, ignore.case = TRUE) ~ "Soft",
    grepl("中|一般|Medium", x, ignore.case = TRUE) ~ "Medium",
    grepl("硬|Hard", x, ignore.case = TRUE) ~ "Hard",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("Soft", "Medium", "Hard"), ordered = TRUE)
}

sub_col <- detect_subtype_col(df)

df_plot <- df %>%
  mutate(
    Subtype = factor(trimws(as.character(.data[[sub_col]])), levels=c("C1","C2","C3")),
    Texture = clean_texture(.data[[resolve_col(df, c("肿瘤质地", "质地"))]]),
    BrainstemAdhesion = factor(clean_yes_no(.data[[resolve_col(df, c("脑干面粘连"))]]), levels=c("No","Yes"))
  ) %>%
  filter(!is.na(Subtype))

# ------------------------------------------------------------------------------
# 4. 绘图函数定义 (原生组件，绝对防撞车)
# ------------------------------------------------------------------------------
plot_clinical_bar <- function(data, var, colors, title, tag_label) {
  clean_data <- data %>% filter(!is.na(Subtype) & !is.na(!!sym(var)))
  p_val <- fisher.test(table(clean_data$Subtype, clean_data[[var]]))$p.value
  
  clean_data %>%
    count(Subtype, !!sym(var)) %>%
    group_by(Subtype) %>%
    mutate(prop = n / sum(n)) %>%
    ggplot(aes(x = Subtype, y = prop, fill = !!sym(var))) +
    geom_col(width = 0.7, color = "white", linewidth = 0.5) + 
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(title = title, subtitle = paste0("Fisher's p = ", sprintf("%.3f", p_val)), 
         x = NULL, y = "Proportion", tag = tag_label) + # 原生加入 ABCD 标签
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right", 
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      plot.tag = element_text(size = 18, face = "bold"), # 控制 ABCD 标签
      axis.text.x = element_text(face = "bold", color = "black")
    )
}

# 🌟 核心修复：完全使用原生 ggplot2 组件，丢弃 ggdist
plot_clinical_violin <- function(data, var, val_map, colors, title, y_lab, tag_label) {
  df_rc <- data %>%
    filter(!is.na(Subtype) & !is.na(!!sym(var))) %>%
    mutate(value = as.numeric(factor(!!sym(var), levels = names(val_map))))
  
  kw_p <- kruskal.test(value ~ Subtype, data = df_rc)$p.value
  
  ggplot(df_rc, aes(x = Subtype, y = value, fill = Subtype, color = Subtype)) +
    # 原生小提琴图：adjust=1.5 用于平滑离散数据，彻底解决 bandwidth 报错
    geom_violin(trim = FALSE, adjust = 1.5, alpha = 0.4, color = NA) +
    # 原生箱线图：置于中心，纯白填充黑框
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 1, color = "black", fill = "white", linewidth = 0.5) +
    # 原生抖动散点
    geom_jitter(width = 0.1, alpha = 0.6, size = 1.8) +
    
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(breaks = 1:length(val_map), labels = names(val_map)) +
    labs(title = title, subtitle = paste0("Kruskal-Wallis p = ", sprintf("%.3f", kw_p)), 
         x = NULL, y = y_lab, tag = tag_label) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none", 
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      plot.tag = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(face = "bold", color = "black")
    )
}

# ------------------------------------------------------------------------------
# 5. 执行绘图并拼图保存 (GridExtra)
# ------------------------------------------------------------------------------
# A. 质地 (Texture)
p1 <- plot_clinical_bar(df_plot, "Texture", PAL_TEXTURE, "Tumor Texture Composition", "A")
p2 <- plot_clinical_violin(df_plot, "Texture", 
                           c("Soft"=1, "Medium"=2, "Hard"=3), PAL_SUBTYPE, 
                           "Texture Grade Distribution", "Texture", "B")

# B. 粘连 (Adhesion)
p3 <- plot_clinical_bar(df_plot, "BrainstemAdhesion", PAL_ADHESION, "Brainstem Adhesion Rate", "C")
p4 <- plot_clinical_violin(df_plot, "BrainstemAdhesion", 
                           c("No"=1, "Yes"=2), PAL_SUBTYPE, 
                           "Adhesion Tendency", "Adhesion", "D")

# C. 拼合 (使用原生底层的 gridExtra::arrangeGrob，无视任何版本冲突)
combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2)

save_path <- file.path(FIG_DIR, "Fig_Texture_Adhesion_Combined_Native.pdf")
ggsave(save_path, plot = combined_plot, width = 12, height = 10, device = cairo_pdf)

message("✅ 完美修复！质地与粘连图表已使用底层原生包生成，请查看: ", save_path)




# ==============================================================================
# Step 7-Clin Addon v16: True Native Source-Driven Sankey (Zero Mismatch)
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. 包加载与环境配置
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
  library(ggalluvial)
})

PROJECT_DIR <- getwd()
FIG_DIR     <- file.path(PROJECT_DIR, "Figures", "Step7_Clinical")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. 终极超级调色盘 (Master Palette)
# ------------------------------------------------------------------------------
# 将所有可能出现的类别合并到一个调色盘中，让原生引擎自己去抓取颜色
MASTER_PALETTE <- c(
  # 分子亚型
  "C1" = "#00A087", "C2" = "#4DBBD5", "C3" = "#E64B35",
  # KOOS 分级
  "I" = "#DEEBF7", "II" = "#9ECAE1", "III" = "#4292C6", "IV" = "#084594",
  # 肿瘤质地
  "Soft" = "#B2E2E2", "Medium" = "#FEB24C", "Hard" = "#E31A1C",
  # 血供 (带后缀区分)
  "No (Blood)" = "#F0F0F0", "Yes (Blood)" = "#FC9272",
  # 粘连 (带后缀区分)
  "No (Adhesion)" = "#F0F0F0", "Yes (Adhesion)" = "#9E9AC8",
  # 切除结局
  "GTR" = "#8DD3C7", "NTR" = "#BEBADA", "STR" = "#FB8072", "PR" = "#D9D9D9"
)

# ------------------------------------------------------------------------------
# 3. 稳健的数据读取与清洗
# ------------------------------------------------------------------------------
df <- NULL
if (file.exists("VS 手术分析.xlsx")) df <- readxl::read_xlsx("VS 手术分析.xlsx", sheet = 1)
if (is.null(df) && file.exists("cl_bulk_with_Subtype.xlsx")) df <- readxl::read_xlsx("cl_bulk_with_Subtype.xlsx", sheet = 1)
if (is.null(df) && file.exists("RData/Clinical_with_Subtype.RData")) {
  load("RData/Clinical_with_Subtype.RData")
  if (exists("cl_data_with_subtype")) df <- cl_data_with_subtype
}
if (is.null(df)) stop("❌ 未找到临床数据！")

# 修复空列名
colnames(df) <- trimws(stringr::str_squish(colnames(df)))
empty_names <- which(colnames(df) == "" | is.na(colnames(df)))
if (length(empty_names) > 0) colnames(df)[empty_names] <- paste0("X_empty_", seq_along(empty_names))

resolve_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) stop("Missing column: ", paste(candidates, collapse = ", "))
  return(hit)
}

detect_subtype_col <- function(df) {
  candidates <- c("Subtype", "group", "Group", "Unnamed: 0", "...1", "X_empty_1")
  hit <- candidates[candidates %in% colnames(df)][1]
  if (!is.na(hit)) return(hit)
  detect_cols <- sapply(df, function(col) mean(trimws(as.character(col)) %in% c("C1", "C2", "C3"), na.rm = TRUE))
  best_idx <- which.max(detect_cols)
  if (length(best_idx) == 1 && detect_cols[best_idx] >= 0.3) return(colnames(df)[best_idx])
  stop("❌ 无法自动识别 Subtype 列！")
}

clean_yes_no <- function(x) case_when(grepl("有|Yes|Y|1|丰富|富", x, ignore.case=T) ~ "Yes", TRUE ~ "No")
clean_koos   <- function(x) {
  x <- toupper(gsub("级|^KOOS\\s*", "", x, ignore.case=T))
  factor(case_when(grepl("IV|4", x)~"IV", grepl("III|3", x)~"III", grepl("II|2", x)~"II", TRUE~"I"), levels=c("I","II","III","IV"))
}
clean_texture <- function(x) {
  x <- trimws(as.character(x))
  factor(case_when(grepl("软|Soft", x, ignore.case = TRUE) ~ "Soft",
                   grepl("中|一般|Medium", x, ignore.case = TRUE) ~ "Medium",
                   grepl("硬|Hard", x, ignore.case = TRUE) ~ "Hard", TRUE ~ "Medium"), 
         levels = c("Soft", "Medium", "Hard"))
}
clean_resec <- function(x) {
  factor(case_when(grepl("近全|NTR", x)~"NTR", grepl("次全|STR", x)~"STR", grepl("部分|大部分|PR", x)~"PR", TRUE~"GTR"), levels=c("GTR","NTR","STR","PR"))
}

# 🌟 核心技巧：对血供和粘连的 Yes/No 加上后缀，防止颜色冲突
df_plot <- df %>%
  mutate(
    Subtype = factor(trimws(as.character(.data[[detect_subtype_col(df)]])), levels=c("C1","C2","C3")),
    KOOS_clean = clean_koos(.data[[resolve_col(df, c("KOOS", "Koos"))]]),
    Texture = clean_texture(.data[[resolve_col(df, c("肿瘤质地", "质地"))]]),
    # 强制加上隐藏后缀，这样原生引擎就能给血供和粘连分别涂上橙色和紫色
    BloodSupply = factor(paste0(clean_yes_no(.data[[resolve_col(df, c("血供丰富", "血供"))]]), " (Blood)"), levels=c("No (Blood)","Yes (Blood)")),
    BrainstemAdhesion = factor(paste0(clean_yes_no(.data[[resolve_col(df, c("脑干面粘连"))]]), " (Adhesion)"), levels=c("No (Adhesion)","Yes (Adhesion)")),
    ResectionClean = clean_resec(.data[[resolve_col(df, c("切除情况", "切除", "切除状态"))]])
  ) %>%
  filter(!is.na(Subtype), !is.na(KOOS_clean), !is.na(ResectionClean))

# ------------------------------------------------------------------------------
# 4. 终极原生绘图 (100% 杜绝错位)
# ------------------------------------------------------------------------------
p <- ggplot(df_plot,
            aes(axis1 = Subtype, axis2 = KOOS_clean, axis3 = Texture, 
                axis4 = BloodSupply, axis5 = BrainstemAdhesion, axis6 = ResectionClean)) +
  
  # A. 流动色带: 水流根据 Subtype 自动平滑追踪
  geom_alluvium(aes(fill = Subtype), width = 0.22, alpha = 0.55, knot.pos = 0.4) +
  
  # B. 原生方块: 完美贴合水流！全列采用统一的高级黑线边框
  geom_stratum(aes(fill = after_stat(stratum)), color = "black", linewidth = 0.6, width = 0.22) +
  
  # C. 智能避让文本标签 (提取原生文本并去处之前的防冲突后缀)
  geom_text(stat = "stratum", 
            aes(
              # 根据轴的位置，智能决定文字是在左侧、中间还是右侧
              x = after_stat(x) + case_when(after_stat(x) == 1 ~ -0.15, after_stat(x) == 6 ~ 0.15, TRUE ~ 0),
              hjust = case_when(after_stat(x) == 1 ~ 1, after_stat(x) == 6 ~ 0, TRUE ~ 0.5),
              # 正则表达式把 " (Blood)" 这种后缀从图表上抹掉，只显示 "Yes" 或 "No"
              label = gsub(" \\(.*\\)", "", after_stat(stratum))
            ),
            fontface = "bold", size = 4.2, color = "black") +
  
  # D. 颜色映射与图例控制 (只显示 C1/C2/C3 图例)
  scale_fill_manual(values = MASTER_PALETTE, name = "Molecular Subtype Origin", breaks = c("C1", "C2", "C3")) +
  
  # E. 全局高级排版
  scale_x_discrete(limits = c("Molecular Subtype", "KOOS Stage", "Tumor Texture", "Blood Supply", "Adhesion", "Resection Outcome"), 
                   expand = c(0.12, 0.12)) + # 左右加宽，防止外侧文本被切断
  labs(y = "Number of Cases", x = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "bottom", legend.box = "horizontal", panel.grid = element_blank()
  )

print(p)

save_path <- file.path(FIG_DIR, "Fig_Native_Source_Driven_Sankey.pdf")
ggsave(save_path, p, width = 14.5, height = 7, device = cairo_pdf)
message("✅ 完美重构完毕！水流与方块已100%对齐，绝不错位，图片已保存至: ", save_path)





# ==============================================================================
# Step 8 V2.3: Immune-Clinical Phenotype Correlation Analysis
# (Fixed Cartesian Product / Point Explosion Bug)
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. 核心包加载
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(readxl)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 2. 路径设置与配色方案
# ------------------------------------------------------------------------------
file_blood   <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/cl_bulk_with_Subtype.xlsx"
file_surgery <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/VS 手术分析.xlsx"

output_dir <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/Step8_Blood_Clinical_Cor/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

PAL_TEXTURE  <- c("Soft" = "#B2E2E2", "Medium" = "#FEB24C", "Hard" = "#E31A1C")
PAL_ADHESION <- c("No" = "#A6CEE3", "Yes" = "#FB9A99")

# ------------------------------------------------------------------------------
# 3. 🌟 核心修复：智能多表读取与 1对1 强制合并 🌟
# ------------------------------------------------------------------------------
df_blood <- read_excel(file_blood)
if(any(grepl("--|×10\\^", as.character(df_blood[1, ])))) {
  df_blood <- df_blood[-1, ]
}
colnames(df_blood) <- trimws(stringr::str_squish(colnames(df_blood)))
df_raw <- df_blood

candidates_texture <- c("肿瘤质地", "质地", "Texture", "texture")
if (!any(candidates_texture %in% colnames(df_raw))) {
  message("⚠️ 外周血表中未找到【肿瘤质地】，正在自动合并【手术记录表】...")
  
  if (file.exists(file_surgery)) {
    df_surgery <- read_excel(file_surgery)
    colnames(df_surgery) <- trimws(stringr::str_squish(colnames(df_surgery)))
    
    # 智能寻找真正的 患者ID 列
    id_candidates <- c("SampleID", "sample_id", "Name", "ID", "患者编号", "样本编号", "编号", "Unnamed: 0", "...1")
    
    id_blood <- id_candidates[id_candidates %in% colnames(df_blood)][1]
    if (is.na(id_blood)) id_blood <- colnames(df_blood)[1] # 兜底方案
    
    id_surg <- id_candidates[id_candidates %in% colnames(df_surgery)][1]
    if (is.na(id_surg)) id_surg <- colnames(df_surgery)[1] # 兜底方案
    
    message("🔍 锁定匹配列: 血液表 [", id_blood, "] <---> 手术表 [", id_surg, "]")
    
    # 强制去重，防止任何笛卡尔积（数据爆炸）
    df_blood_uniq <- df_blood %>% distinct(!!sym(id_blood), .keep_all = TRUE)
    df_surg_uniq  <- df_surgery %>% distinct(!!sym(id_surg), .keep_all = TRUE)
    
    df_raw <- left_join(df_blood_uniq, df_surg_uniq, by = setNames(id_surg, id_blood))
    message("✅ 成功完成 1对1 合并！当前总样本量: ", nrow(df_raw), " 例")
  } else {
    stop("❌ 找不到手术记录表，无法获取质地数据！")
  }
}

# ------------------------------------------------------------------------------
# 4. 智能模糊匹配字典
# ------------------------------------------------------------------------------
blood_dict <- list(
  "NK Cells"      = c("自然杀伤细胞", "NK细胞", "NK"),
  "Albumin"       = c("白蛋白", "ALB", "Albumin"),
  "Age"           = c("年龄", "Age", "age", "AGE"),
  "MPV"           = c("平均血小板体积", "MPV", "平均血小板容积"),
  "Total T Cells" = c("总T淋巴细胞", "总T细胞", "T淋巴细胞", "Total T"),
  "eGFR"          = c("估算肾小球滤过率(CKD-EPI)", "估算肾小球滤过率", "eGFR", "肾小球滤过率"),
  "A/G Ratio"     = c("白/球比", "白球比", "A/G", "A/G Ratio")
)

found_cn_cols <- c()
found_en_cols <- c()

for (en_name in names(blood_dict)) {
  candidates <- blood_dict[[en_name]]
  hit <- candidates[candidates %in% colnames(df_raw)][1]
  if (!is.na(hit)) {
    found_cn_cols <- c(found_cn_cols, hit)
    found_en_cols <- c(found_en_cols, en_name)
  } else {
    message("⚠️ 警告: 未找到指标 [", en_name, "]，将自动跳过。")
  }
}

# ------------------------------------------------------------------------------
# 5. 数据清洗与特征工程
# ------------------------------------------------------------------------------
resolve_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) stop("Missing column: ", paste(candidates, collapse = ", "))
  return(hit)
}

clean_numeric <- function(x) as.numeric(trimws(gsub("<|>", "", as.character(x))))
clean_yes_no <- function(x) case_when(grepl("有|Yes|Y|1|丰富|富", x, ignore.case=T) ~ "Yes", TRUE ~ "No")
clean_texture <- function(x) {
  x <- trimws(as.character(x))
  factor(case_when(grepl("软|Soft", x, ignore.case = T) ~ "Soft",
                   grepl("中|一般|Medium", x, ignore.case = T) ~ "Medium",
                   grepl("硬|Hard", x, ignore.case = T) ~ "Hard", TRUE ~ NA_character_), 
         levels = c("Soft", "Medium", "Hard"), ordered = TRUE)
}

col_texture <- resolve_col(df_raw, c("肿瘤质地", "质地", "Texture"))
col_adhesion <- resolve_col(df_raw, c("脑干面粘连", "脑干粘连", "粘连"))

df_clean <- df_raw %>%
  mutate(across(any_of(found_cn_cols), clean_numeric)) %>%
  rename_with(~ found_en_cols, all_of(found_cn_cols))

if ("NK Cells" %in% colnames(df_clean) && "Total T Cells" %in% colnames(df_clean)) {
  df_clean <- df_clean %>% mutate(`NK/Total T Ratio` = `NK Cells` / `Total T Cells`)
  final_blood_vars <- c(found_en_cols, "NK/Total T Ratio")
} else {
  final_blood_vars <- found_en_cols
}

df_clean <- df_clean %>%
  mutate(
    Texture = clean_texture(.data[[col_texture]]),
    BrainstemAdhesion = factor(clean_yes_no(.data[[col_adhesion]]), levels = c("No", "Yes")),
    Texture_Score = as.numeric(Texture),             
    Adhesion_Score = ifelse(BrainstemAdhesion == "Yes", 1, 0) 
  )

# ------------------------------------------------------------------------------
# 6. 分析 A：全局 Spearman 相关性气泡图 
# ------------------------------------------------------------------------------
message("⏳ 正在计算外周血与临床表型的 Spearman 相关性...")

cor_results <- data.frame()
for (clin in c("Texture_Score", "Adhesion_Score")) {
  for (blood in final_blood_vars) {
    tmp <- df_clean %>% dplyr::select(all_of(c(clin, blood))) %>% drop_na()
    if (nrow(tmp) > 10) { 
      res <- cor.test(tmp[[1]], tmp[[2]], method = "spearman", exact = FALSE)
      cor_results <- rbind(cor_results, data.frame(
        Clinical_Trait = ifelse(clin == "Texture_Score", "Tumor Texture", "Brainstem Adhesion"),
        Blood_Marker = blood,
        Rho = res$estimate,
        P_value = res$p.value
      ))
    }
  }
}

if (nrow(cor_results) > 0) {
  cor_results <- cor_results %>%
    mutate(
      Significance = case_when(P_value < 0.001 ~ "***", P_value < 0.01 ~ "**", P_value < 0.05 ~ "*", TRUE ~ ""),
      Blood_Marker = factor(Blood_Marker, levels = rev(sort(unique(Blood_Marker))))
    )
  
  p_cor <- ggplot(cor_results, aes(x = Clinical_Trait, y = Blood_Marker)) +
    geom_point(aes(size = -log10(P_value), fill = Rho), shape = 21, color = "grey30", stroke = 0.5) +
    geom_text(aes(label = Significance), vjust = 0.7, hjust = 0.5, size = 5, color = "black", fontface = "bold") +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0, 
                         name = "Spearman Rho\n(Correlation)", limits = c(-max(abs(cor_results$Rho)), max(abs(cor_results$Rho)))) +
    scale_size_continuous(name = "-log10(P-value)", range = c(3, 10)) +
    labs(x = NULL, y = NULL, title = "Correlation: Peripheral Blood vs. Surgical Phenotype") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  ggsave(paste0(output_dir, "Fig1_Correlation_Bubble_Map.pdf"), p_cor, width = 7.5, height = 6.5, device = cairo_pdf)
}

# ------------------------------------------------------------------------------
# 7. 分析 B & C：临床分组差异小提琴图 
# ------------------------------------------------------------------------------
message("⏳ 正在绘制临床分组的小提琴图...")

plot_long <- df_clean %>%
  dplyr::select(Texture, BrainstemAdhesion, all_of(final_blood_vars)) %>%
  pivot_longer(cols = all_of(final_blood_vars), names_to = "Indicator", values_to = "Value") %>%
  drop_na(Value) %>%
  mutate(Indicator = factor(Indicator, levels = final_blood_vars))

# 【图 B】按脑干粘连分组
p_adhesion <- ggplot(plot_long %>% filter(!is.na(BrainstemAdhesion)), aes(x = BrainstemAdhesion, y = Value, fill = BrainstemAdhesion)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.0, alpha = 0.5, color = "gray20") +
  scale_fill_manual(values = PAL_ADHESION) +
  facet_wrap(~ Indicator, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(), legend.position = "none", panel.grid.major.x = element_blank()
  ) +
  labs(y = "Blood Level / Ratio", title = "Peripheral Blood Markers stratified by Brainstem Adhesion") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("No", "Yes")), label = "p.format", color = "firebrick")

ggsave(paste0(output_dir, "Fig2_Blood_by_Adhesion_Violin.pdf"), p_adhesion, width = 12, height = 8.5, device = cairo_pdf)

# 【图 C】按肿瘤质地分组
p_texture <- ggplot(plot_long %>% filter(!is.na(Texture)), aes(x = Texture, y = Value, fill = Texture)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.0, alpha = 0.5, color = "gray20") +
  scale_fill_manual(values = PAL_TEXTURE) +
  facet_wrap(~ Indicator, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(), legend.position = "none", panel.grid.major.x = element_blank()
  ) +
  labs(y = "Blood Level / Ratio", title = "Peripheral Blood Markers stratified by Tumor Texture") +
  stat_compare_means(method = "kruskal.test", label.y.npc = "bottom", size = 3.5, color = "firebrick") +
  stat_compare_means(comparisons = list(c("Soft", "Hard"), c("Medium", "Hard")), method = "wilcox.test", label = "p.signif", hide.ns = TRUE)

ggsave(paste0(output_dir, "Fig3_Blood_by_Texture_Violin.pdf"), p_texture, width = 12, height = 8.5, device = cairo_pdf)



message("✅ 分析完美完成！数据笛卡尔积漏洞已彻底修复，图表已保存至: ", output_dir)


# ==============================================================================
# Step 8 V5.3: Pre-filtered Correlation Correlogram & Raincloud Plots
# (Ultimate Base R Rename Fix - 100% Bulletproof)
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. 核心包加载
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(stringr)
  library(ggpubr)
  library(ggcorrplot) 
  library(ggdist)     
})

# ------------------------------------------------------------------------------
# 2. 路径设置
# ------------------------------------------------------------------------------
file_blood   <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/cl_bulk_with_Subtype.xlsx"
file_surgery <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/VS 手术分析.xlsx"

output_dir <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/Figures/Step8_Filtered_Visuals/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

PAL_ADHESION <- c("No" = "#4DBBD5", "Yes" = "#E64B35")
PAL_TEXTURE  <- c("Soft" = "#B2E2E2", "Medium" = "#FEB24C", "Hard" = "#E31A1C")

# ------------------------------------------------------------------------------
# 3. 智能多表读取与数据合并
# ------------------------------------------------------------------------------
df_blood <- read_excel(file_blood)
if(any(grepl("--|×10\\^", as.character(df_blood[1, ])))) df_blood <- df_blood[-1, ]
colnames(df_blood) <- trimws(stringr::str_squish(colnames(df_blood)))
df_raw <- df_blood

find_col <- function(df, patterns) {
  for (p in patterns) {
    hit <- grep(p, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  return(NA_character_)
}

if (is.na(find_col(df_raw, c("肿瘤质地", "质地", "Texture")))) {
  df_surgery <- read_excel(file_surgery)
  colnames(df_surgery) <- trimws(stringr::str_squish(colnames(df_surgery)))
  
  id_candidates <- c("SampleID", "sample_id", "Name", "ID", "患者编号", "样本编号", "编号")
  id_blood <- find_col(df_blood, id_candidates) %||% colnames(df_blood)[1]
  id_surg  <- find_col(df_surgery, id_candidates) %||% colnames(df_surgery)[1]
  
  df_raw <- left_join(
    df_blood %>% distinct(!!sym(id_blood), .keep_all = TRUE), 
    df_surgery %>% distinct(!!sym(id_surg), .keep_all = TRUE), 
    by = setNames(id_surg, id_blood)
  )
}

# ------------------------------------------------------------------------------
# 4. 数据清洗与原表特征提取 
# ------------------------------------------------------------------------------
blood_dict <- list(
  "NK Cells"      = c("自然杀伤细胞", "NK细胞", "NK"),
  "Albumin"       = c("白蛋白", "ALB", "Albumin"),
  "Age"           = c("年龄", "Age", "age"),
  "MPV"           = c("平均血小板体积", "MPV"),
  "Total T Cells" = c("总T淋巴细胞", "总T细胞", "T淋巴细胞"),
  "eGFR"          = c("估算肾小球滤过率(CKD-EPI)", "eGFR"),
  "A/G Ratio"     = c("白/球比", "白球比", "A/G Ratio")
)

found_cn <- c(); found_en <- c()
for (en_name in names(blood_dict)) {
  hit <- find_col(df_raw, blood_dict[[en_name]])
  if (!is.na(hit)) { found_cn <- c(found_cn, hit); found_en <- c(found_en, en_name) }
}

clean_numeric <- function(x) as.numeric(trimws(gsub("<|>", "", as.character(x))))
clean_yes_no <- function(x) ifelse(grepl("有|Yes|Y|1|丰富|富", as.character(x), ignore.case=T), "Yes", "No")

col_texture <- find_col(df_raw, c("肿瘤质地", "质地", "Texture"))
col_adhesion <- find_col(df_raw, c("脑干面粘连", "脑干粘连", "粘连"))

if (is.na(col_texture) || is.na(col_adhesion)) stop("❌ 找不到肿瘤质地或脑干粘连列，请检查 Excel 表头！")

df_clean <- df_raw %>%
  mutate(across(any_of(found_cn), clean_numeric))

# 🔴 防御：使用 base R 强制改名
for (i in seq_along(found_cn)) {
  colnames(df_clean)[colnames(df_clean) == found_cn[i]] <- found_en[i]
}

df_clean <- df_clean %>%
  mutate(
    Texture = factor(case_when(grepl("软|Soft", .data[[col_texture]]) ~ "Soft", grepl("中", .data[[col_texture]]) ~ "Medium", grepl("硬|Hard", .data[[col_texture]]) ~ "Hard", TRUE ~ NA_character_), levels = c("Soft", "Medium", "Hard"), ordered = TRUE),
    BrainstemAdhesion = factor(clean_yes_no(.data[[col_adhesion]]), levels = c("No", "Yes")),
    Adhesion_Score = ifelse(BrainstemAdhesion == "Yes", 1, 0),
    Texture_Score = as.numeric(Texture)
  )

candidate_blood_vars <- found_en

# ------------------------------------------------------------------------------
# 5. P 值漏斗 (统计学前置筛查)
# ------------------------------------------------------------------------------
message("⏳ 正在进行统计学筛查，寻找显著相关指标...")

p_vals <- data.frame()
for (blood in candidate_blood_vars) {
  tmp <- df_clean %>% dplyr::select(all_of(c("Texture_Score", "Adhesion_Score", blood))) %>% drop_na()
  if (nrow(tmp) > 10) {
    p_tex <- cor.test(tmp[["Texture_Score"]], tmp[[blood]], method = "spearman", exact = FALSE)$p.value
    p_adh <- cor.test(tmp[["Adhesion_Score"]], tmp[[blood]], method = "spearman", exact = FALSE)$p.value
    p_vals <- rbind(p_vals, data.frame(Marker = blood, Min_P = min(p_tex, p_adh, na.rm = T)))
  }
}

sig_markers <- p_vals %>% filter(Min_P < 0.1) %>% pull(Marker)

if (length(sig_markers) == 0) {
  message("⚠️ 警告：未发现显著相关指标。将强制展示 P 值最小的 Top 3 指标。")
  sig_markers <- p_vals %>% arrange(Min_P) %>% head(3) %>% pull(Marker)
}

message("🎯 最终筛选出进入画图环节的指标：", paste(sig_markers, collapse = ", "))

# ------------------------------------------------------------------------------
# 6. 绘图 1：下三角相关性矩阵热图 
# ------------------------------------------------------------------------------
df_cor <- df_clean %>% 
  dplyr::select(all_of(c("Texture_Score", "Adhesion_Score", sig_markers))) %>%
  drop_na()

# 🔴 终极防御：完全抛弃 rename()，使用底层的 colnames() 强制改名
colnames(df_cor)[colnames(df_cor) == "Texture_Score"] <- "Tumor Texture"
colnames(df_cor)[colnames(df_cor) == "Adhesion_Score"] <- "Adhesion"

cor_matrix <- cor(df_cor, method = "spearman")
p_matrix <- cor_pmat(df_cor, method = "spearman")

p_heat <- ggcorrplot(
  cor_matrix,
  method = "square",       
  type = "lower",          
  p.mat = p_matrix,        
  sig.level = 0.05,        
  insig = "pch",           
  pch.col = "grey50",
  outline.col = "white",   
  colors = c("#B2182B", "white", "#2166AC"), 
  lab = TRUE,              
  lab_size = 3.5,
  tl.col = "black",        
  tl.srt = 45              
) +
  ggtitle("Correlogram: Selected Markers & Surgical Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

ggsave(paste0(output_dir, "Plot1_Filtered_Correlogram.pdf"), p_heat, width = 7.5, height = 7.5)



# ------------------------------------------------------------------------------
# 7. 绘图 2：针对筛选指标的原生特写图 (Native Violin + Boxplot, 100% 防崩溃)
# ------------------------------------------------------------------------------
message("⏳ 正在绘制差异分布图...")

df_long <- df_clean %>%
  dplyr::select(all_of(c("BrainstemAdhesion", "Texture", sig_markers))) %>%
  pivot_longer(cols = all_of(sig_markers), names_to = "Indicator", values_to = "Value") %>%
  drop_na(Value) %>%
  mutate(Indicator = factor(Indicator, levels = sig_markers))

# 🌟 原生小提琴图 (针对粘连) - 彻底抛弃 ggdist
p_rain_adh <- ggplot(df_long %>% filter(!is.na(BrainstemAdhesion)), 
                     aes(x = BrainstemAdhesion, y = Value, fill = BrainstemAdhesion, color = BrainstemAdhesion)) +
  # 1. 原生小提琴图 (替代云朵，adjust=1.2 让曲线更平滑)
  geom_violin(trim = FALSE, adjust = 1.2, alpha = 0.35, color = NA) +
  # 2. 原生箱线图 (纯白填充，黑边框)
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 1, color = "black", fill = "white", linewidth = 0.5) +
  # 3. 抖动散点
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = PAL_ADHESION) +
  scale_color_manual(values = PAL_ADHESION) +
  facet_wrap(~ Indicator, scales = "free_y") +
  theme_classic(base_size = 13) +
  labs(y = "Marker Level", x = "Brainstem Adhesion", title = "Significant Markers by Adhesion Status") +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", color = "black")
  ) +
  # 🌟 修复警告：加入 exact = FALSE
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(exact = FALSE), 
                     comparisons = list(c("No", "Yes")), 
                     label = "p.format", color = "firebrick")

ggsave(paste0(output_dir, "Plot2_Violin_Adhesion.pdf"), p_rain_adh, width = length(sig_markers)*2.5 + 2, height = 5)

# 🌟 原生小提琴图 (针对质地)
p_rain_tex <- ggplot(df_long %>% filter(!is.na(Texture)), 
                     aes(x = Texture, y = Value, fill = Texture, color = Texture)) +
  geom_violin(trim = FALSE, adjust = 1.2, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 1, color = "black", fill = "white", linewidth = 0.5) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = PAL_TEXTURE) +
  scale_color_manual(values = PAL_TEXTURE) +
  facet_wrap(~ Indicator, scales = "free_y") +
  theme_classic(base_size = 13) +
  labs(y = "Marker Level", x = "Tumor Texture", title = "Significant Markers by Texture") +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", color = "black")
  ) +
  # 🌟 Kruskal-Wallis 本身不需要 exact，Wilcoxon 加入 exact = FALSE
  stat_compare_means(method = "kruskal.test", label.y.npc = "bottom", size = 3.5, color = "firebrick") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(exact = FALSE),
                     comparisons = list(c("Soft", "Hard"), c("Medium", "Hard")), 
                     label = "p.signif", hide.ns = TRUE)

ggsave(paste0(output_dir, "Plot3_Violin_Texture.pdf"), p_rain_tex, width = length(sig_markers)*2.5 + 2, height = 5)

message("✅ 图表更新完毕！已使用原生引擎彻底解决崩溃问题！")





# ------------------------------------------------------------------------------
# 6. 绘图 1：下三角相关性矩阵热图 (顶级定制版：展示 R 值 + 星号 + 具体的 P 值)
# ------------------------------------------------------------------------------
df_cor <- df_clean %>% 
  dplyr::select(all_of(c("Texture_Score", "Adhesion_Score", sig_markers))) %>%
  drop_na()

# 使用底层 colnames() 强制改名
colnames(df_cor)[colnames(df_cor) == "Texture_Score"] <- "Tumor Texture"
colnames(df_cor)[colnames(df_cor) == "Adhesion_Score"] <- "Adhesion"

# 计算相关系数和 P 值矩阵
cor_matrix <- cor(df_cor, method = "spearman")
p_matrix <- cor_pmat(df_cor, method = "spearman")

# 将上三角和对角线的数据设为 NA (只保留下三角)
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
p_matrix[upper.tri(p_matrix, diag = TRUE)] <- NA

# 🌟 终极防御：完全抛弃 as.table 和 rename，使用 pivot_longer 强行指定新列名
cor_df <- as.data.frame(cor_matrix) %>%
  rownames_to_column(var = "Var1") %>%
  pivot_longer(cols = -Var1, names_to = "Var2", values_to = "Rho") %>%
  drop_na()

p_df <- as.data.frame(p_matrix) %>%
  rownames_to_column(var = "Var1") %>%
  pivot_longer(cols = -Var1, names_to = "Var2", values_to = "Pval") %>%
  drop_na()

plot_df <- left_join(cor_df, p_df, by = c("Var1", "Var2"))

# 🌟 核心定制：生成包含 R值、星号、P值的完美文本标签
plot_df <- plot_df %>%
  mutate(
    # 计算显著性星号
    Signif = case_when(Pval < 0.001 ~ "***", Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ ""),
    # 格式化 P 值 (小于0.001显示<0.001，否则保留3位小数)
    P_format = ifelse(Pval < 0.001, "<0.001", sprintf("%.3f", Pval)),
    # 拼接最终文字标签
    Label = sprintf("%.2f%s\n(p=%s)", Rho, Signif, P_format),
    # 控制坐标轴的阶梯状排列顺序
    Var1 = factor(Var1, levels = rev(colnames(cor_matrix))),
    Var2 = factor(Var2, levels = colnames(cor_matrix))
  )

# 使用原生 ggplot2 绘制高级热图
p_heat <- ggplot(plot_df, aes(x = Var2, y = Var1, fill = Rho)) +
  geom_tile(color = "white", linewidth = 1.2) + # 白色网格线分割
  # 智能文字颜色：当相关系数绝对值 > 0.45（颜色很深）时，字变白；否则字为黑
  geom_text(aes(label = Label, color = abs(Rho) > 0.45), size = 3.8, fontface = "bold", lineheight = 0.9) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  # 经典的红白蓝渐变配色
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limit = c(-1, 1), 
                       space = "Lab", name = "Spearman Rho\nCorrelation") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15))
  ) +
  coord_fixed() + # 保证画出来的必定是正方形方块
  ggtitle("Correlogram: Selected Markers & Surgical Phenotypes")

ggsave(paste0(output_dir, "Plot1_Filtered_Correlogram_with_Pvalue.pdf"), p_heat, width = 8, height = 8)
message("✅ 带具体 P 值和星号的高级相关性热图已生成！")





# ==============================================================================
# Multivariate_Age_Adjustment.R
# 功能：合并手术临床数据与围手术期实验室数据，
#       控制年龄后验证亚型对核心结局的独立预测作用
# ==============================================================================

library(tidyverse)
library(readxl)
library(broom)

# ── 路径配置（与你现有代码一致）──────────────────────────────────────────────
source("/media/desk16/iy5111/VS.sc/vs.paper/VS.4/00_config.R")

input_lab  <- "/media/desk16/iy5111/VS.sc/vs.paper/VS.4/cl_bulk_with_Subtype.xlsx"
OUT_DIR    <- file.path(TABLE_DIR, "Multivariate_AgeAdjust")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. 加载实验室数据（与 38_blood_R.txt 逻辑相同）──────────────────────────
clean_numeric <- function(x) {
  as.numeric(gsub("<|>", "", trimws(as.character(x))))
}

df_lab <- read_excel(input_lab)

# # 去掉参考值行（如果存在）
# if (any(grepl("--|×10\\^", as.character(df_lab[1, ])))) {
#   df_lab <- df_lab[-1, ]
# }
# 
# # 提取所需列，重命名为英文
# df_lab_clean <- df_lab %>%
#   select(
#     SampleID    = SampleID,
#     Subtype     = Subtype,
#     Age         = Age,
#     NK_Cells    = `自然杀伤细胞`,
#     Albumin     = `白蛋白`,
#     MPV         = `平均血小板体积`,
#     TotalT      = `总T淋巴细胞`,
#     eGFR        = `估算肾小球滤过率(CKD-EPI)`
#   ) %>%
#   filter(Subtype %in% c("C1", "C2", "C3")) %>%
#   mutate(
#     across(c(Age, NK_Cells, Albumin, MPV, TotalT, eGFR), clean_numeric),
#     Subtype = factor(Subtype, levels = c("C1", "C2", "C3"))
#   )


# 去掉参考值行（如果存在）
if (any(grepl("--|×10\\^", as.character(df_lab[1, ])))) {
  df_lab <- df_lab[-1, ]
}

# 提取所需列，重命名为英文（注意这里加上了 dplyr::）
df_lab_clean <- df_lab %>%
  dplyr::select(
    SampleID    = SampleID,
    Subtype     = Subtype,
    Age         = Age,
    NK_Cells    = `自然杀伤细胞`,
    Albumin     = `白蛋白`,
    MPV         = `平均血小板体积`,
    TotalT      = `总T淋巴细胞`,
    eGFR        = `估算肾小球滤过率(CKD-EPI)`
  ) %>%
  dplyr::filter(Subtype %in% c("C1", "C2", "C3")) %>%
  dplyr::mutate(
    dplyr::across(c(Age, NK_Cells, Albumin, MPV, TotalT, eGFR), clean_numeric),
    Subtype = factor(Subtype, levels = c("C1", "C2", "C3"))
  )

cat("✅ 实验室数据加载完成，共", nrow(df_lab_clean), "例\n")
cat("各亚型：\n"); print(table(df_lab_clean$Subtype))

# ── 2. 加载手术临床数据（来自 术中.R 的 df_plot）────────────────────────────
# 如果你的手术表保存在 cl_bulk_with_Subtype.xlsx 同一文件里，直接用 df_lab
# 如果单独在 VS 手术分析.xlsx 里，则读取并合并

# 检查手术相关列是否已在 df_lab_clean
surgery_cols <- c("BrainstemAdhesion", "Texture", "脑干面粘连", "肿瘤质地")
present <- intersect(surgery_cols, colnames(df_lab))
cat("\n检测到的手术相关列:", paste(present, collapse = ", "), "\n")

# 如果没有，从手术分析文件加载并合并
if (length(present) == 0) {
  surgery_files <- c(
    file.path(PROJECT_DIR, "VS 手术分析.xlsx"),
    file.path(PROJECT_DIR, "VS 手术分析(1).xlsx")
  )
  surg_file <- surgery_files[file.exists(surgery_files)][1]
  
  if (!is.na(surg_file)) {
    df_surg_raw <- read_excel(surg_file, sheet = 1)
    colnames(df_surg_raw) <- trimws(colnames(df_surg_raw))
    
    # 清洗手术数据
    clean_yn <- function(x) {
      case_when(
        grepl("有|Yes|Y|1|丰富|粘连", x, ignore.case = TRUE) ~ "Yes",
        grepl("无|No|N|0|否", x, ignore.case = TRUE)         ~ "No",
        TRUE ~ NA_character_
      )
    }
    clean_tex <- function(x) {
      case_when(
        grepl("软|Soft", x, ignore.case = TRUE)   ~ 1L,
        grepl("中|Medium", x, ignore.case = TRUE) ~ 2L,
        grepl("硬|Hard", x, ignore.case = TRUE)   ~ 3L,
        TRUE ~ NA_integer_
      )
    }
    
    # 自动识别样本ID列
    id_col <- intersect(c("SampleID", "sample_id", "Name"), colnames(df_surg_raw))[1]
    adh_col <- intersect(c("脑干面粘连", "BrainstemAdhesion"), colnames(df_surg_raw))[1]
    tex_col <- intersect(c("肿瘤质地", "Texture"), colnames(df_surg_raw))[1]
    
    df_surg <- df_surg_raw %>%
      transmute(
        SampleID          = as.character(.data[[id_col]]),
        BrainstemAdhesion = factor(clean_yn(.data[[adh_col]]), levels = c("No","Yes")),
        Texture_num       = clean_tex(.data[[tex_col]])
      ) %>%
      filter(!is.na(SampleID))
    
    # 合并
    df_merged <- df_lab_clean %>%
      left_join(df_surg, by = "SampleID")
    
    cat("✅ 手术数据合并完成，共", nrow(df_merged), "例\n")
    
  } else {
    stop("❌ 未找到手术分析文件，请检查路径")
  }
} else {
  # 手术数据已在同一文件中
  clean_yn <- function(x) {
    case_when(
      grepl("有|Yes|Y|1|粘连", x, ignore.case = TRUE) ~ "Yes",
      grepl("无|No|N|0|否", x, ignore.case = TRUE)    ~ "No",
      TRUE ~ NA_character_
    )
  }
  
  adh_col_raw <- intersect(c("脑干面粘连","BrainstemAdhesion"), colnames(df_lab))[1]
  tex_col_raw <- intersect(c("肿瘤质地","Texture"), colnames(df_lab))[1]
  
  df_merged <- df_lab_clean %>%
    mutate(
      BrainstemAdhesion = factor(clean_yn(df_lab[[adh_col_raw]]), levels = c("No","Yes")),
      Texture_num = case_when(
        grepl("软|Soft", df_lab[[tex_col_raw]], ignore.case = TRUE) ~ 1L,
        grepl("中|Medium", df_lab[[tex_col_raw]], ignore.case = TRUE) ~ 2L,
        grepl("硬|Hard", df_lab[[tex_col_raw]], ignore.case = TRUE) ~ 3L,
        TRUE ~ NA_integer_
      )
    )
}

# ── 3. 多因素校正分析 ──────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════\n")
cat("多因素分析（控制年龄后，亚型的独立预测作用）\n")
cat("══════════════════════════════════════════════════\n\n")

results_list <- list()

# ── 3a. 脑干面粘连（二分类 → Logistic Regression）────────────────────────────
cat("【1】脑干面粘连 | Logistic Regression (Age + Subtype)\n")
df_adh <- df_merged %>% filter(!is.na(BrainstemAdhesion), !is.na(Age), !is.na(Subtype))

fit_adh <- glm(BrainstemAdhesion ~ Age + Subtype, data = df_adh, family = binomial)
res_adh <- tidy(fit_adh, conf.int = TRUE, exponentiate = TRUE) %>%  # OR
  mutate(Outcome = "BrainstemAdhesion", Model = "Logistic")
print(res_adh)
results_list[["BrainstemAdhesion"]] <- res_adh

# ── 3b. 肿瘤质地（有序数值 → Linear Regression）─────────────────────────────
cat("\n【2】肿瘤质地（Soft=1/Medium=2/Hard=3）| Linear Regression\n")
df_tex <- df_merged %>% filter(!is.na(Texture_num), !is.na(Age), !is.na(Subtype))

fit_tex <- lm(Texture_num ~ Age + Subtype, data = df_tex)
res_tex <- tidy(fit_tex, conf.int = TRUE) %>%
  mutate(Outcome = "Texture", Model = "Linear")
print(res_tex)
results_list[["Texture"]] <- res_tex

# ── 3c. NK 细胞（连续变量 → Linear Regression）──────────────────────────────
cat("\n【3】NK 细胞 | Linear Regression (Age + Subtype)\n")
df_nk <- df_merged %>% filter(!is.na(NK_Cells), !is.na(Age), !is.na(Subtype))

fit_nk <- lm(NK_Cells ~ Age + Subtype, data = df_nk)
res_nk <- tidy(fit_nk, conf.int = TRUE) %>%
  mutate(Outcome = "NK_Cells", Model = "Linear")
print(res_nk)
results_list[["NK_Cells"]] <- res_nk

# ── 3d. 白蛋白 ──────────────────────────────────────────────────────────────
cat("\n【4】白蛋白（Albumin）| Linear Regression\n")
df_alb <- df_merged %>% filter(!is.na(Albumin), !is.na(Age), !is.na(Subtype))

fit_alb <- lm(Albumin ~ Age + Subtype, data = df_alb)
res_alb <- tidy(fit_alb, conf.int = TRUE) %>%
  mutate(Outcome = "Albumin", Model = "Linear")
print(res_alb)
results_list[["Albumin"]] <- res_alb

# ── 3e. eGFR ────────────────────────────────────────────────────────────────
cat("\n【5】eGFR | Linear Regression\n")
df_egfr <- df_merged %>% filter(!is.na(eGFR), !is.na(Age), !is.na(Subtype))

fit_egfr <- lm(eGFR ~ Age + Subtype, data = df_egfr)
res_egfr <- tidy(fit_egfr, conf.int = TRUE) %>%
  mutate(Outcome = "eGFR", Model = "Linear")
print(res_egfr)
results_list[["eGFR"]] <- res_egfr

# ── 4. 汇总导出 ─────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════\n")
cat("汇总：年龄校正后各指标亚型效应\n")
cat("══════════════════════════════════════════════════\n\n")

# all_results <- bind_rows(results_list) %>%
#   filter(grepl("Subtype", term)) %>%
#   select(Outcome, term, estimate, conf.low, conf.high, p.value) %>%
#   mutate(
#     Significant = ifelse(p.value < 0.05, "Yes *", "No"),
#     p_fmt = case_when(
#       p.value < 0.001 ~ sprintf("%.2e", p.value),
#       p.value < 0.01  ~ sprintf("%.4f", p.value),
#       TRUE            ~ sprintf("%.3f", p.value)
#     )
#   )
all_results <- dplyr::bind_rows(results_list) %>%
  dplyr::filter(grepl("Subtype", term)) %>%
  dplyr::select(Outcome, term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate(
    Significant = ifelse(p.value < 0.05, "Yes *", "No"),
    p_fmt = dplyr::case_when(
      p.value < 0.001 ~ sprintf("%.2e", p.value),
      p.value < 0.01  ~ sprintf("%.4f", p.value),
      TRUE            ~ sprintf("%.3f", p.value)
    )
  )
print(all_results)
write.csv(all_results, file.path(OUT_DIR, "Multivariate_AgeAdjust_Summary.csv"), row.names = FALSE)
cat("\n✅ 结果已保存至:", file.path(OUT_DIR, "Multivariate_AgeAdjust_Summary.csv"), "\n")












# 