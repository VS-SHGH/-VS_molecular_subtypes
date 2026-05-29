# ################################################################################
# #           前庭神经鞘瘤 (VS) - GSE216781 专项补全脚本 (RNA-seq)                #
# # 功能：
# # 1. 仅下载/重新处理 GSE216781
# # 2. 增强样本名匹配逻辑，找回丢失的样本
# # 3. 保持原有路径和输出文件名不变
# ################################################################################
# 
# rm(list = ls())
# gc()
# 
# # =============================================================================
# # 1. 环境与路径设置 (保持不变)
# # =============================================================================
# work_dir <- "/media/desk16/iy5111/VS.sc/vs.paper"
# if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
# setwd(work_dir)
# 
# # 加载必要包
# pkgs <- c("GEOquery", "limma", "edgeR", "dplyr", "tidyr", "stringr", "data.table", 
#           "AnnotationDbi", "org.Hs.eg.db")
# for (pkg in pkgs) { if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update =F, ask=F) }
# lapply(pkgs, require, character.only = TRUE)
# 
# # 读取元数据
# meta_filename <- "bulk_annotations_filtered.csv"
# if (!file.exists(meta_filename)) {
#   stop("❌ 错误: 在当前目录下找不到 'bulk_annotations_filtered.csv'。请确保该文件存在。")
# }
# meta <- fread(meta_filename, data.table = FALSE)
# 
# # 输出路径设置
# out_root <- "Bulk_Data_Cleaned"
# raw_data_dir <- "raw_data_cache"
# if (!dir.exists(out_root)) dir.create(out_root)
# if (!dir.exists(raw_data_dir)) dir.create(raw_data_dir)
# 
# # =============================================================================
# # 2. GSE216781 核心处理逻辑
# # =============================================================================
# message("\n>>> [Target Fix] 正在重新处理 GSE216781 ...")
# 
# # 1. 确定目标样本 (从元数据中提取)
# target_samples <- meta %>% filter(Study == "Present") %>% pull(Sample_geo_accession)
# message(paste0("📊 元数据中包含该数据集样本数: ", length(target_samples)))
# 
# # 2. 准备保存目录
# save_dir <- file.path(out_root, "GSE216781_Present")
# if (!dir.exists(save_dir)) dir.create(save_dir)
# 
# # 3. 下载/读取数据 (GSE216781 属于 Supplementary Files)
# getGEOSuppFiles("GSE216781", baseDir = raw_data_dir, makeDirectory = TRUE)
# supp_path <- file.path(raw_data_dir, "GSE216781")
# 
# # 寻找 Count Matrix 文件
# count_file <- list.files(supp_path, pattern = "count|matrix", full.names = TRUE)
# # 优先选 raw_counts 相关的 csv/txt
# count_file <- count_file[grep("csv|txt|tsv", count_file)][1]
# 
# if (!is.na(count_file)) {
#   message("📖 正在读取原始矩阵: ", basename(count_file))
#   
#   # !!! 关键修改: check.names = FALSE 防止 R 将 '-' 自动变为 '.' 导致匹配失败
#   raw_counts <- fread(count_file, data.table = FALSE, check.names = FALSE) 
#   
#   # 处理行名 (第一列通常是 Gene ID)
#   rownames(raw_counts) <- raw_counts[,1]
#   raw_counts <- raw_counts[,-1]
#   
#   # 4. 处理列名映射 (Internal ID -> GSM ID)
#   # GSE216781 通常附带一个 mapping 文件，或者是列名本身需要清洗
#   map_file <- list.files(supp_path, pattern = "IDs|mapping", full.names = TRUE)[1]
#   
#   if (!is.na(map_file)) {
#     message("🔄 检测到 ID 映射文件，正在校正样本名...")
#     mapping <- fread(map_file, data.table = FALSE)
#     
#     # 获取当前矩阵的列名
#     curr_cols <- colnames(raw_counts)
#     new_cols <- curr_cols
#     
#     # 建立映射字典 (假设 mapping 文件有 ID_1=GSM, ID_2=Internal)
#     # 这里的逻辑稍微放宽，防止部分字符匹配不上
#     for(i in 1:nrow(mapping)) {
#       geo_id <- mapping[i, 1] # 假设第一列是 GSM
#       int_id <- mapping[i, 2] # 假设第二列是 TL-xx
#       
#       # 精确匹配
#       idx <- which(curr_cols == int_id)
#       
#       # 如果精确匹配失败，尝试模糊匹配 (处理 R 可能读入的差异)
#       if(length(idx) == 0) {
#         # 尝试将 - 换成 . 再找
#         int_id_dot <- gsub("-", "\\.", int_id)
#         idx <- which(curr_cols == int_id_dot)
#       }
#       
#       if(length(idx) > 0) {
#         new_cols[idx] <- geo_id
#       }
#     }
#     colnames(raw_counts) <- new_cols
#   }
#   
#   # 5. 样本交集与丢失样本诊断 (Debug 环节)
#   valid_samps <- intersect(colnames(raw_counts), target_samples)
#   
#   # --- 诊断信息 ---
#   missing_samps <- setdiff(target_samples, colnames(raw_counts))
#   if(length(missing_samps) > 0) {
#     message("⚠️  警告: 以下样本在元数据中存在，但在表达矩阵中未找到 (请检查是否因名称不一致导致):")
#     print(missing_samps)
#   } else {
#     message("✅ 所有目标样本均已匹配成功。")
#   }
#   # ----------------
#   
#   if(length(valid_samps) > 0) {
#     counts_sub <- raw_counts[, valid_samps, drop=FALSE]
#     message(paste0("⚙️  正在处理 ", ncol(counts_sub), " 个样本的基因注释与标准化..."))
#     
#     # 6. Ensembl ID -> Symbol
#     # 去除版本号 (例如 ENSG000001.3 -> ENSG000001)
#     ens_ids <- sub("\\..*", "", rownames(counts_sub)) 
#     
#     # 转换
#     syms <- mapIds(org.Hs.eg.db, keys = ens_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
#     
#     # 将 Symbol 赋给矩阵
#     counts_sub$Symbol <- syms
#     
#     # 去除无法注释的基因 (NA)
#     counts_sub <- counts_sub[!is.na(counts_sub$Symbol), ]
#     
#     # 7. 处理重复基因 Symbol (取平均值)
#     # 注意: 先提取数值列进行 avereps
#     mat_only <- counts_sub[, -ncol(counts_sub)]
#     counts_final <- limma::avereps(mat_only, ID = counts_sub$Symbol)
#     
#     # 8. 标准化 (EdgeR: TMM -> Log2 CPM)
#     dge <- DGEList(counts = counts_final)
#     dge <- calcNormFactors(dge) # TMM 标准化因子
#     expr_log2 <- cpm(dge, log = TRUE, prior.count = 2) # Log2 转换
#     
#     # 9. 保存结果
#     out_file <- file.path(save_dir, "GSE216781_Expression_Log2.csv")
#     write.csv(expr_log2, out_file)
#     
#     message("----------------------------------------------------------------")
#     message("✅ GSE216781 修正完成！")
#     message("📁 文件已保存至: ", out_file)
#     message("📊 最终矩阵维度: ", nrow(expr_log2), " genes x ", ncol(expr_log2), " samples")
#     message("----------------------------------------------------------------")
#     
#   } else {
#     stop("❌ 严重错误: 样本匹配后数量为 0。请检查 mapping 文件列名是否反了，或者 raw counts 的列名格式。")
#   }
#   
# } else {
#   stop("❌ 错误: 在 raw_data_cache/GSE216781 中未找到 counts 或 matrix 文件。下载可能未完成。")
# }






################################################################################
#           前庭神经鞘瘤 (VS) - GSE216781 最终修正 (强制找回丢失样本版)       #
################################################################################

rm(list = ls())
gc()

# =============================================================================
# 1. 设置路径与环境
# =============================================================================
work_dir <- "/media/desk16/iy5111/VS.sc/vs.paper"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

pkgs <- c("GEOquery", "limma", "edgeR", "dplyr", "data.table", "AnnotationDbi", "org.Hs.eg.db")
lapply(pkgs, require, character.only = TRUE)

# 读取元数据
meta <- fread("bulk_annotations_filtered.csv", data.table = FALSE)
target_samples <- meta %>% filter(Study == "Present") %>% pull(Sample_geo_accession)
message("🎯 目标样本总数: ", length(target_samples), " (包含 GSM6692737)")

# =============================================================================
# 2. 读取并清洗数据 (FeatureCounts 格式处理)
# =============================================================================
raw_data_dir <- "raw_data_cache"
supp_path <- file.path(raw_data_dir, "GSE216781")
count_file <- list.files(supp_path, pattern = "count|matrix", full.names = TRUE)[1]

message("📖 读取文件: ", basename(count_file))
raw_df <- fread(count_file, data.table = FALSE, check.names = FALSE)

# 设置行名 (GeneID)
rownames(raw_df) <- raw_df[, 1]

# 剔除 FeatureCounts 的注释列 (Chr, Start, End, Strand, Length)
# 以及第一列 Geneid 本身
cols_to_remove <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
is_annotation <- colnames(raw_df) %in% cols_to_remove
counts_matrix <- raw_df[, !is_annotation]

message("📊 原始矩阵清洗后维度: ", nrow(counts_matrix), " genes x ", ncol(counts_matrix), " samples")

# =============================================================================
# 3. 智能匹配与“强制修正” (核心步骤)
# =============================================================================
map_file <- list.files(supp_path, pattern = "IDs|mapping", full.names = TRUE)[1]
mapping <- fread(map_file, data.table = FALSE)

# 1. 构建字典 (Internal ID -> GSM ID)
dict_keys <- gsub("-", ".", mapping[,2]) # 将 - 转为 . 以适应 R 的列名
dict_values <- mapping[,1]
id_map <- setNames(dict_values, dict_keys)

# 2. 尝试初步匹配
current_cols <- colnames(counts_matrix)
matched_cols <- id_map[current_cols] # 这里会有一列是 NA (因为匹配不上)

# 3. 执行强制替换逻辑
# 找到那个匹配不上的列 (NA)
na_indices <- which(is.na(matched_cols))

# 找到那个元数据里有、但还没匹配上的目标 GSM ID
# 我们先看看哪些已经匹配上了
matched_gsms <- na.omit(matched_cols)
missing_target <- setdiff(target_samples, matched_gsms)

if (length(na_indices) == 1 && length(missing_target) == 1) {
  # 场景确认：正好有一个列不知道是谁，也正好缺一个目标样本
  unknown_col_name <- current_cols[na_indices]
  target_gsm <- missing_target[1] # 应该是 GSM6692737
  
  message("\n-------------------------------------------------------------")
  message("🛠️  触发强制匹配机制 (User Request):")
  message("   检测到未知列名: [ ", unknown_col_name, " ]")
  message("   检测到缺失目标: [ ", target_gsm, " ] (对应 TL-19-550BBF)")
  message("   ✅ 已将该列强制重命名为 ", target_gsm)
  message("-------------------------------------------------------------")
  
  # 填补空缺
  matched_cols[na_indices] <- target_gsm
  
} else if (length(na_indices) > 0) {
  # 如果有多个匹配不上，或者逻辑对不上，就不敢乱改，打印警告
  message("⚠️ 警告: 无法自动执行一对一强制匹配。")
  message("未知列数量: ", length(na_indices))
  message("缺失目标数量: ", length(missing_target))
}

# 4. 赋值新的列名 (现在应该是完整的 GSM 编号了)
colnames(counts_matrix) <- matched_cols

# 5. 最终筛选 (Double Check)
# 这一步理论上会保留所有 22 个样本
final_samples <- intersect(colnames(counts_matrix), target_samples)
counts_final <- counts_matrix[, final_samples, drop=FALSE]

message("✅ 最终保留样本数: ", ncol(counts_final), " / ", length(target_samples))
if(ncol(counts_final) == length(target_samples)) {
  message("🎉 完美！所有 22 个样本（包括找回的那个）都已就位。")
}

# =============================================================================
# 4. 标准化与保存
# =============================================================================
# Ensembl -> Symbol
ens_ids <- sub("\\..*", "", rownames(counts_final))
syms <- mapIds(org.Hs.eg.db, keys = ens_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
counts_final$Symbol <- syms
counts_final <- counts_final[!is.na(counts_final$Symbol), ]

# 重复基因取平均
mat_only <- counts_final[, -ncol(counts_final)]
counts_clean <- limma::avereps(mat_only, ID = counts_final$Symbol)

# Log2 CPM 标准化
dge <- DGEList(counts = counts_clean)
dge <- calcNormFactors(dge)
expr_log2 <- cpm(dge, log = TRUE, prior.count = 2)

# 保存
save_dir <- "Bulk_Data_Cleaned/GSE216781_Present"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

out_file <- file.path(save_dir, "GSE216781_Expression_Log2_22.csv")
write.csv(expr_log2, out_file)

message("\n💾 文件已保存至: ", out_file)