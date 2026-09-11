library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readxl)

# Define file paths
cibersort_file <- "/Users/yanchen/Desktop/VS.4/supplement/all/step8_Immune/CIBERSORT_Result.csv"
mcp_file <- "/Users/yanchen/Desktop/VSbulk+单细胞/VS_revision/output/mcp_counter_results.rds"
clinical_file <- "/Users/yanchen/Desktop/VSbulk+单细胞/数据/cl_bulk_with_Subtype.xlsx"
out_dir <- "/Users/yanchen/Desktop/VSbulk+单细胞/VS_revision/figures/"

# 1. Load Data
message("Loading data...")
cibersort <- read.csv(cibersort_file, row.names=1)
mcp <- readRDS(mcp_file)
cl <- read_excel(clinical_file)

# 2. Process Clinical Data
cl_clean <- cl %>% 
  filter(!is.na(SampleID)) %>%
  mutate(SampleID_full = paste0("tpm.", SampleID)) %>%
  select(SampleID_full, Subtype)

# 3. Process MCP-counter Data
# mcp is a matrix with cell types as rows and samples as columns
mcp_df <- as.data.frame(t(mcp))
mcp_df$SampleID_full <- rownames(mcp_df)

# 4. Process CIBERSORT Data
cibersort$SampleID_full <- rownames(cibersort)

# Calculate derived CIBERSORT scores to match MCP-counter
cibersort_derived <- cibersort %>%
  mutate(
    `T cells CD8` = `T.cells.CD8`,
    `B lineage` = `B.cells.naive` + `B.cells.memory`,
    `Monocytic lineage` = `Macrophages.M2`,
    `NK cells` = `NK.cells.resting` + `NK.cells.activated`
  ) %>%
  select(SampleID_full, `T cells CD8`, `B lineage`, `Monocytic lineage`, `NK cells`)

# 5. Merge Data
merged_data <- mcp_df %>%
  select(SampleID_full, `CD8 T cells`, `B lineage`, `Monocytic lineage`, `NK cells`) %>%
  rename(
    MCP_CD8_T_cells = `CD8 T cells`,
    MCP_B_lineage = `B lineage`,
    MCP_Monocytic_lineage = `Monocytic lineage`,
    MCP_NK_cells = `NK cells`
  ) %>%
  inner_join(cibersort_derived, by = "SampleID_full") %>%
  rename(
    CIBERSORT_CD8_T_cells = `T cells CD8`,
    CIBERSORT_B_lineage = `B lineage`,
    CIBERSORT_Monocytic_lineage = `Monocytic lineage`,
    CIBERSORT_NK_cells = `NK cells`
  ) %>%
  inner_join(cl_clean, by = "SampleID_full")

# Make Subtype a factor
merged_data$Subtype <- factor(merged_data$Subtype, levels = c("C1", "C2", "C3"))

# Define colors
subtype_colors <- c("C1" = "#1B9E77", "C2" = "#1F78B4", "C3" = "#D95F02")

# 6. Plotting Function
create_scatter <- function(data, x_col, y_col, x_label, y_label, title) {
  p <- ggscatter(data, x = x_col, y = y_col, color = "Subtype",
            palette = subtype_colors, size = 3, alpha = 0.8,
            add = "reg.line", add.params = list(color = "black", linetype = "dashed")) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(x = x_label, y = y_label, title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right"
    )
  return(p)
}

# 7. Create Plots
p1 <- create_scatter(merged_data, "CIBERSORT_CD8_T_cells", "MCP_CD8_T_cells", 
                     "CIBERSORT (CD8 T cells)", "MCP-counter (CD8 T cells)", 
                     "CD8 T cells Validation")

p2 <- create_scatter(merged_data, "CIBERSORT_B_lineage", "MCP_B_lineage", 
                     "CIBERSORT (B naive + memory)", "MCP-counter (B lineage)", 
                     "B lineage Validation")

p3 <- create_scatter(merged_data, "CIBERSORT_Monocytic_lineage", "MCP_Monocytic_lineage", 
                     "CIBERSORT (Macrophages M2)", "MCP-counter (Monocytic lineage)", 
                     "Monocytic Lineage Validation")

p4 <- create_scatter(merged_data, "CIBERSORT_NK_cells", "MCP_NK_cells", 
                     "CIBERSORT (NK resting + activated)", "MCP-counter (NK cells)", 
                     "NK cells Validation")

# Combine plots
combined_plot <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# 8. Save Plot
output_file <- paste0(out_dir, "Figure2i_CIBERSORT_vs_MCPcounter.pdf")
ggsave(output_file, combined_plot, width = 10, height = 10)
message("Saved plot to: ", output_file)

# Save individual plots just in case
ggsave(paste0(out_dir, "Figure2i_CD8_T.pdf"), p1, width = 5, height = 5)
ggsave(paste0(out_dir, "Figure2i_B_lineage.pdf"), p2, width = 5, height = 5)
ggsave(paste0(out_dir, "Figure2i_Monocytic.pdf"), p3, width = 5, height = 5)
ggsave(paste0(out_dir, "Figure2i_NK_cells.pdf"), p4, width = 5, height = 5)
message("Saved individual plots as well.")
