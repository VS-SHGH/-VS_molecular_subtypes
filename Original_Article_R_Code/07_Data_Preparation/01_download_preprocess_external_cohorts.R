################################################################################
#           前庭神经鞘瘤 (VS) 公共数据整理脚本 (仅下载与清洗)                 #
# 功能：
# 1. 下载 GSE216781 (RNA-seq), GSE141801 (Array), GSE39645 (Array)
# 2. 统一转换为 Gene Symbol 行名
# 3. 统一进行 Log2 标准化
# 4. 保存为整洁的 CSV 文件，供后续分析使用
################################################################################

# =============================================================================
# 1. 环境准备与文件检查
# =============================================================================
rm(list = ls())
gc()

# 设置工作目录
work_dir <- "/media/desk16/iy5111/VS.sc/vs.paper"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# --- 自动检查并“抢救”注释文件 ---
meta_filename <- "bulk_annotations_filtered.csv"
if (!file.exists(meta_filename)) {
  message("⚠️ 当前目录下未找到 ", meta_filename, "，正在尝试从用户主目录搜索...")
  # 尝试在 iy5111 用户目录下搜索（限制搜索深度以防太慢）
  found <- list.files("/media/desk16/iy5111", pattern = meta_filename, recursive = TRUE, full.names = TRUE)
  # 排除掉缓存或回收站里的文件
  found <- found[!grepl("Trash|cache", found)]
  
  if (length(found) > 0) {
    message("✅ 找到文件: ", found[1])
    file.copy(found[1], paste0(work_dir, "/", meta_filename))
    message("已自动将其复制到工作目录。")
  } else {
    stop("\n❌ 严重错误: 无法自动找到 'bulk_annotations_filtered.csv'。\n请您手动将该文件上传或移动到: ", work_dir)
  }
}

# 加载包
pkgs <- c("GEOquery", "limma", "edgeR", "dplyr", "tidyr", "stringr", "data.table", 
          "AnnotationDbi", "org.Hs.eg.db", "hgu133plus2.db", "hugene10sttranscriptcluster.db", "hgu219.db")
for (pkg in pkgs) { if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update =F, ask=F) }
lapply(pkgs, require, character.only = TRUE)

# 读取元数据
meta <- fread(meta_filename, data.table = FALSE)

# 创建输出总目录
out_root <- "Bulk_Data_Cleaned"
if (!dir.exists(out_root)) dir.create(out_root)
raw_data_dir <- "raw_data_cache" # 下载缓存目录，避免重复下载
if (!dir.exists(raw_data_dir)) dir.create(raw_data_dir)

# =============================================================================
# 2. 通用函数定义
# =============================================================================

# 芯片探针转基因 Symbol 函数
process_array_data <- function(gset, target_samples, db_package) {
  # 提取表达矩阵
  ex <- exprs(gset)
  
  # 1. 样本筛选
  common <- intersect(colnames(ex), target_samples)
  if(length(common) == 0) return(NULL)
  ex <- ex[, common, drop=FALSE]
  
  # 2. Log2 转换检查 (自动判断是否需要取Log)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex <- log2(ex + 1) 
    message("   -> 检测到原始数值，已执行 Log2 转换。")
  } else {
    message("   -> 检测到数据已 Log2 化，跳过转换。")
  }
  
  # 3. 探针 ID 转 Symbol
  require(db_package, character.only = TRUE)
  probe_ids <- rownames(ex)
  symbols <- mapIds(get(db_package), keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  
  # 移除无注释探针
  valid_idx <- !is.na(symbols) & symbols != ""
  ex <- ex[valid_idx, ]
  symbols <- symbols[valid_idx]
  
  # 4. 重复基因取平均
  ex_clean <- limma::avereps(ex, ID = symbols)
  return(ex_clean)
}

# =============================================================================
# 3. 处理数据集 1: GSE216781 (Present / RNA-seq)
# =============================================================================
message("\n>>> [1/3] 正在处理 GSE216781 (Present Study)...")
target_samples <- meta %>% filter(Study == "Present") %>% pull(Sample_geo_accession)
save_dir <- file.path(out_root, "GSE216781_Present")
if (!dir.exists(save_dir)) dir.create(save_dir)

# 下载并读取 Supplementary Files
getGEOSuppFiles("GSE216781", baseDir = raw_data_dir, makeDirectory = TRUE)
supp_path <- file.path(raw_data_dir, "GSE216781")
count_file <- list.files(supp_path, pattern = "count|matrix", full.names = TRUE)[1]

if (!is.na(count_file)) {
  raw_counts <- fread(count_file, data.table = FALSE)
  rownames(raw_counts) <- raw_counts[,1]; raw_counts <- raw_counts[,-1]
  
  # 尝试匹配列名 (处理 TL-xx vs GSMxx 的映射)
  map_file <- list.files(supp_path, pattern = "IDs|mapping", full.names = TRUE)[1]
  if (!is.na(map_file)) {
    mapping <- fread(map_file, data.table = FALSE)
    curr_cols <- colnames(raw_counts)
    new_cols <- curr_cols
    for(i in 1:nrow(mapping)) {
      # 替换逻辑：文件里的内部ID (如 TL-19-001) -> GEO ID (GSM...)
      pattern <- gsub("-", "\\.", mapping[i, "ID_2"]) # 处理 R 读入时将 - 变 . 的情况
      idx <- which(str_detect(curr_cols, pattern))
      if(length(idx) > 0) new_cols[idx] <- mapping[i, "ID_1"]
    }
    colnames(raw_counts) <- new_cols
  }
  
  # 筛选样本
  valid_samps <- intersect(colnames(raw_counts), target_samples)
  if(length(valid_samps) > 0) {
    counts_sub <- raw_counts[, valid_samps, drop=FALSE]
    
    # Ensembl ID -> Symbol
    ens_ids <- sub("\\..*", "", rownames(counts_sub)) # 去版本号
    syms <- mapIds(org.Hs.eg.db, keys = ens_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    
    # 聚合
    counts_sub$Symbol <- syms
    counts_sub <- counts_sub[!is.na(counts_sub$Symbol), ]
    counts_final <- limma::avereps(counts_sub[, -ncol(counts_sub)], ID = counts_sub$Symbol)
    
    # 标准化 (Log2 CPM)
    dge <- DGEList(counts = counts_final)
    dge <- calcNormFactors(dge)
    expr_log2 <- cpm(dge, log = TRUE, prior.count = 2)
    
    # 保存
    write.csv(expr_log2, file.path(save_dir, "GSE216781_Expression_Log2.csv"))
    message("✅ GSE216781 处理完成！")
  } else {
    message("⚠️ 警告: GSE216781 样本匹配失败，请检查 Mapping 文件逻辑。")
  }
} else {
  message("⚠️ 警告: 未找到 GSE216781 的 Count 矩阵文件。")
}

# =============================================================================
# 4. 处理数据集 2: GSE141801 (Gugel / Array)
# =============================================================================
message("\n>>> [2/3] 正在处理 GSE141801 (Gugel et al)...")
target_samples <- meta %>% filter(Study == "Gugel et al 2020") %>% pull(Sample_geo_accession)
save_dir <- file.path(out_root, "GSE141801_Gugel")
if (!dir.exists(save_dir)) dir.create(save_dir)

tryCatch({
  gset <- getGEO("GSE141801", destdir = raw_data_dir, getGPL = FALSE)[[1]]
  # 使用 Affy U219 注释包
  expr_clean <- process_array_data(gset, target_samples, "hgu219.db")
  
  if(!is.null(expr_clean)){
    write.csv(expr_clean, file.path(save_dir, "GSE141801_Expression_Log2.csv"))
    message("✅ GSE141801 处理完成！")
  }
}, error = function(e) { message("❌ GSE141801 处理出错: ", e$message) })

# =============================================================================
# 5. 处理数据集 3: GSE39645 (Torres / Array)
# =============================================================================
message("\n>>> [3/3] 正在处理 GSE39645 (Torres-Martin et al)...")
target_samples <- meta %>% filter(Study == "Torres-Martin et al 2013") %>% pull(Sample_geo_accession)
save_dir <- file.path(out_root, "GSE39645_Torres")
if (!dir.exists(save_dir)) dir.create(save_dir)

tryCatch({
  gset <- getGEO("GSE39645", destdir = raw_data_dir, getGPL = FALSE)[[1]]
  # 使用 HuGene 1.0 ST 注释包
  expr_clean <- process_array_data(gset, target_samples, "hugene10sttranscriptcluster.db")
  
  if(!is.null(expr_clean)){
    write.csv(expr_clean, file.path(save_dir, "GSE39645_Expression_Log2.csv"))
    message("✅ GSE39645 处理完成！")
  }
}, error = function(e) { message("❌ GSE39645 处理出错: ", e$message) })

message("\n################################################################")
message("🎉 数据整理全部完成！")
message("结果保存在: ", file.path(work_dir, out_root))
message("包含三个子文件夹，每个里面都有整理好的 'Expression_Log2.csv'")
message("################################################################")

################################################################################
#           前庭神经鞘瘤 (VS) 公共数据整理脚本 (下载、清洗 + 保存临床信息)      #
################################################################################

# =============================================================================
# 1. 环境准备与文件检查
# =============================================================================
rm(list = ls())
gc()

# 设置工作目录
work_dir <- "/media/desk16/iy5111/VS.sc/vs.paper"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# --- 自动检查并“抢救”注释文件 ---
meta_filename <- "bulk_annotations_filtered.csv"
if (!file.exists(meta_filename)) {
  message("⚠️ 当前目录下未找到 ", meta_filename, "，正在尝试从用户主目录搜索...")
  found <- list.files("/media/desk16/iy5111", pattern = meta_filename, recursive = TRUE, full.names = TRUE)
  found <- found[!grepl("Trash|cache", found)]
  
  if (length(found) > 0) {
    message("✅ 找到文件: ", found[1])
    file.copy(found[1], paste0(work_dir, "/", meta_filename))
    message("已自动将其复制到工作目录。")
  } else {
    stop("\n❌ 严重错误: 无法自动找到 'bulk_annotations_filtered.csv'。\n请您手动将该文件上传或移动到: ", work_dir)
  }
}

# 加载包
pkgs <- c("GEOquery", "limma", "edgeR", "dplyr", "tidyr", "stringr", "data.table", 
          "AnnotationDbi", "org.Hs.eg.db", "hgu133plus2.db", "hugene10sttranscriptcluster.db", "hgu219.db")
for (pkg in pkgs) { if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update =F, ask=F) }
lapply(pkgs, require, character.only = TRUE)

# 读取元数据
meta <- fread(meta_filename, data.table = FALSE)

# 创建输出总目录
out_root <- "Bulk_Data_Cleaned"
if (!dir.exists(out_root)) dir.create(out_root)
raw_data_dir <- "raw_data_cache" 
if (!dir.exists(raw_data_dir)) dir.create(raw_data_dir)

# =============================================================================
# 2. 通用函数定义
# =============================================================================

# 芯片探针转基因 Symbol 函数
process_array_data <- function(gset, target_samples, db_package) {
  ex <- exprs(gset)
  
  # 1. 样本筛选
  common <- intersect(colnames(ex), target_samples)
  if(length(common) == 0) return(NULL)
  ex <- ex[, common, drop=FALSE]
  
  # 2. Log2 转换检查
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex <- log2(ex + 1) }
  
  # 3. 探针 ID 转 Symbol
  require(db_package, character.only = TRUE)
  probe_ids <- rownames(ex)
  symbols <- mapIds(get(db_package), keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  
  # 移除无注释探针
  valid_idx <- !is.na(symbols) & symbols != ""
  ex <- ex[valid_idx, ]
  symbols <- symbols[valid_idx]
  
  # 4. 重复基因取平均
  ex_clean <- limma::avereps(ex, ID = symbols)
  return(ex_clean)
}

# 辅助函数：保存临床数据
save_clinical_data <- function(expr_mat, full_meta, save_path) {
  # 找出当前矩阵里的样本
  samples_in_matrix <- colnames(expr_mat)
  # 从总表中筛选这些样本
  clinical_subset <- full_meta %>% filter(Sample_geo_accession %in% samples_in_matrix)
  # 确保顺序一致（可选，但推荐）
  clinical_subset <- clinical_subset[match(samples_in_matrix, clinical_subset$Sample_geo_accession), ]
  
  write.csv(clinical_subset, file.path(save_path, "Clinical_Data.csv"), row.names = FALSE)
  message("   -> 对应的临床数据已保存。")
}

# =============================================================================
# 3. 处理数据集 1: GSE216781 (Present / RNA-seq)
# =============================================================================
message("\n>>> [1/3] 正在处理 GSE216781 (Present Study)...")
target_samples <- meta %>% filter(Study == "Present") %>% pull(Sample_geo_accession)
save_dir <- file.path(out_root, "GSE216781_Present")
if (!dir.exists(save_dir)) dir.create(save_dir)

getGEOSuppFiles("GSE216781", baseDir = raw_data_dir, makeDirectory = TRUE)
supp_path <- file.path(raw_data_dir, "GSE216781")
count_file <- list.files(supp_path, pattern = "count|matrix", full.names = TRUE)[1]

if (!is.na(count_file)) {
  raw_counts <- fread(count_file, data.table = FALSE)
  rownames(raw_counts) <- raw_counts[,1]; raw_counts <- raw_counts[,-1]
  
  map_file <- list.files(supp_path, pattern = "IDs|mapping", full.names = TRUE)[1]
  if (!is.na(map_file)) {
    mapping <- fread(map_file, data.table = FALSE)
    curr_cols <- colnames(raw_counts)
    new_cols <- curr_cols
    for(i in 1:nrow(mapping)) {
      pattern <- gsub("-", "\\.", mapping[i, "ID_2"])
      idx <- which(str_detect(curr_cols, pattern))
      if(length(idx) > 0) new_cols[idx] <- mapping[i, "ID_1"]
    }
    colnames(raw_counts) <- new_cols
  }
  
  valid_samps <- intersect(colnames(raw_counts), target_samples)
  if(length(valid_samps) > 0) {
    counts_sub <- raw_counts[, valid_samps, drop=FALSE]
    ens_ids <- sub("\\..*", "", rownames(counts_sub))
    syms <- mapIds(org.Hs.eg.db, keys = ens_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    counts_sub$Symbol <- syms
    counts_sub <- counts_sub[!is.na(counts_sub$Symbol), ]
    counts_final <- limma::avereps(counts_sub[, -ncol(counts_sub)], ID = counts_sub$Symbol)
    
    dge <- DGEList(counts = counts_final)
    dge <- calcNormFactors(dge)
    expr_log2 <- cpm(dge, log = TRUE, prior.count = 2)
    
    # 保存表达矩阵
    write.csv(expr_log2, file.path(save_dir, "GSE216781_Expression_Log2.csv"))
    # 【新增】保存临床数据
    save_clinical_data(expr_log2, meta, save_dir)
    
    message("✅ GSE216781 处理完成！")
  } else {
    message("⚠️ 警告: GSE216781 样本匹配失败。")
  }
} else {
  message("⚠️ 警告: 未找到 GSE216781 的 Count 矩阵文件。")
}

# =============================================================================
# 4. 处理数据集 2: GSE141801 (Gugel / Array)
# =============================================================================
message("\n>>> [2/3] 正在处理 GSE141801 (Gugel et al)...")
target_samples <- meta %>% filter(Study == "Gugel et al 2020") %>% pull(Sample_geo_accession)
save_dir <- file.path(out_root, "GSE141801_Gugel")
if (!dir.exists(save_dir)) dir.create(save_dir)

tryCatch({
  gset <- getGEO("GSE141801", destdir = raw_data_dir, getGPL = FALSE)[[1]]
  expr_clean <- process_array_data(gset, target_samples, "hgu219.db")
  
  if(!is.null(expr_clean)){
    write.csv(expr_clean, file.path(save_dir, "GSE141801_Expression_Log2.csv"))
    # 【新增】保存临床数据
    save_clinical_data(expr_clean, meta, save_dir)
    message("✅ GSE141801 处理完成！")
  }
}, error = function(e) { message("❌ GSE141801 处理出错: ", e$message) })

# =============================================================================
# 5. 处理数据集 3: GSE39645 (Torres / Array)
# =============================================================================
message("\n>>> [3/3] 正在处理 GSE39645 (Torres-Martin et al)...")
target_samples <- meta %>% filter(Study == "Torres-Martin et al 2013") %>% pull(Sample_geo_accession)
save_dir <- file.path(out_root, "GSE39645_Torres")
if (!dir.exists(save_dir)) dir.create(save_dir)

tryCatch({
  gset <- getGEO("GSE39645", destdir = raw_data_dir, getGPL = FALSE)[[1]]
  expr_clean <- process_array_data(gset, target_samples, "hugene10sttranscriptcluster.db")
  
  if(!is.null(expr_clean)){
    write.csv(expr_clean, file.path(save_dir, "GSE39645_Expression_Log2.csv"))
    # 【新增】保存临床数据
    save_clinical_data(expr_clean, meta, save_dir)
    message("✅ GSE39645 处理完成！")
  }
}, error = function(e) { message("❌ GSE39645 处理出错: ", e$message) })

message("\n################################################################")
message("🎉 数据整理全部完成！")
message("结果保存在: ", file.path(work_dir, out_root))
message("包含三个子文件夹，每个里面都有 'Expression_Log2.csv' 和 'Clinical_Data.csv'")
message("################################################################")
