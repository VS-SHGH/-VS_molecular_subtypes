# =============================================================================
# Generate the current Supplementary Figure 2 from audited outputs.
# Historical proxy decomposition (71.8%) is deliberately not plotted here.
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(42) # Reproducible display jitter; input values are unchanged.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))

out_dir <- file.path(AUDIT_ROOT, "suppfig2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mgs_input <- file.path(AUDIT_ROOT, "mgs", "variance_partitioning_input_official_estimate.csv")
external_input <- file.path(AUDIT_ROOT, "external_validation", "external_mgs_figure_input.csv")
vp_input <- file.path(AUDIT_ROOT, "mgs", "variance_partitioning_steps_official_estimate.csv")
assoc_input <- file.path(AUDIT_ROOT, "mgs", "official_estimate_association_statistics.csv")
required <- c(mgs_input, external_input, vp_input, assoc_input)
if (any(!file.exists(required))) {
  stop("Missing audited Supplementary Figure 2 input(s): ",
       paste(required[!file.exists(required)], collapse = "; "))
}

subtype_cols <- c(C1 = "#2E8B57", C2 = "#20B2AA", C3 = "#CD5C5C")
class_cols <- c("Injury-like" = "#CD5C5C", "nmSC Core" = "#4682B4")
df <- read.csv(mgs_input, check.names = FALSE)
ext <- read.csv(external_input, check.names = FALSE)
vp <- read.csv(vp_input, check.names = FALSE)
assoc <- read.csv(assoc_input, check.names = FALSE)
stopifnot(nrow(df) == 38L)
df$Subtype <- factor(df$Subtype, levels = c("C1", "C2", "C3"))

rho_test <- suppressWarnings(cor.test(df$MGS, df$Tumour_Purity,
                                      method = "spearman", exact = FALSE))
sc_res <- residuals(lm(SC_Score ~ Tumour_Purity, data = df))
mgs_res <- residuals(lm(MGS ~ Tumour_Purity, data = df))
partial_r <- unname(cor(sc_res, mgs_res, method = "pearson"))
partial_p <- assoc$value[assoc$metric ==
                         "SC_vs_MGS_partial_P_controlling_Tumour_Purity"]

panel_theme <- theme_minimal(base_size = 10.5, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        axis.title = element_text(size = 9.5), axis.text = element_text(size = 8.5),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8.5),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))

p_a <- ggplot(df, aes(x = Tumour_Purity, y = MGS, colour = Subtype)) +
  geom_point(size = 2.4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", linewidth = 0.65) +
  scale_colour_manual(values = subtype_cols, drop = FALSE) +
  annotate("text", x = 0.67, y = 1.08,
           label = sprintf("rho = -%.3f\nP < 0.001", abs(unname(rho_test$estimate))),
           hjust = 0, vjust = 1, size = 3.0, colour = "black") +
  labs(x = "ESTIMATE tumour purity", y = "MGS", title = "MGS vs tumour purity") +
  annotate("text", x = 0.57, y = 1.25, label = "a", fontface = "bold", size = 4.2,
           hjust = 0, vjust = 1) +
  panel_theme + theme(legend.position = "none")

df$SC_purity_adjusted <- sc_res
p_b <- ggplot(df, aes(x = MGS, y = SC_purity_adjusted, colour = Subtype)) +
  geom_point(size = 2.4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", linewidth = 0.65) +
  scale_colour_manual(values = subtype_cols, drop = FALSE) +
  annotate("text", x = 0.56, y = 0.92,
           label = sprintf("Partial r = -%.3f\nP = %.4f", abs(partial_r), partial_p),
           hjust = 0, vjust = 1, size = 3.0, colour = "black") +
  labs(x = "MGS", y = "SC score (purity-adjusted)",
       title = "Purity-adjusted SC score") +
  panel_theme + theme(legend.position = "right")

ext$Class <- factor(ext$Class, levels = c("Injury-like", "nmSC Core"))
ext$Dataset <- factor(ext$Dataset, levels = c("GSE141801", "GSE39645"))
make_external_panel <- function(dat, dataset_label) {
  pval <- kruskal.test(MGS ~ Class, data = dat)$p.value
  ggplot(dat, aes(x = Class, y = MGS, fill = Class)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.72, width = 0.6) +
    geom_jitter(width = 0.12, size = 1.3, alpha = 0.55, colour = "black") +
    annotate("text", x = 1.5, y = 0.95,
             label = if (pval < 0.001) "P < 0.001" else sprintf("P = %.2f", pval), size = 3.0) +
    scale_fill_manual(values = class_cols, drop = FALSE) +
    scale_y_continuous(limits = c(-0.05, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = NULL, y = "Molecular Gradient Score", title = dataset_label) +
    panel_theme + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
                         axis.text.x = element_text(size = 8))
}
p_c1 <- make_external_panel(filter(ext, Dataset == "GSE141801"), "GSE141801")
p_c2 <- make_external_panel(filter(ext, Dataset == "GSE39645"), "GSE39645")
p_c1 <- p_c1 + annotate("text", x = 0.65, y = 1.12, label = "b", fontface = "bold", size = 4.2,
                        hjust = 0, vjust = 1)

stromal_r2 <- vp$R2[vp$Model == "Stromal score only"]
total_r2 <- vp$R2[vp$Model == "Official ESTIMATE stromal + immune scores"]
residual <- 1 - total_r2
immune_delta <- total_r2 - stromal_r2
vp_plot <- data.frame(
  Component = factor(c("Null", "+ Stromal score", "+ Immune score", "Residual"),
                     levels = c("Null", "+ Stromal score", "+ Immune score", "Residual")),
  Value = c(0, stromal_r2, total_r2, residual),
  Type = c("Cumulative R2", "Cumulative R2", "Cumulative R2", "Residual")
)
p_d <- ggplot(vp_plot, aes(x = Component, y = Value, fill = Type)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * Value)), vjust = -0.45, size = 3.1) +
  scale_fill_manual(values = c("Cumulative R2" = "#4472C4", "Residual" = "#A5A5A5")) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), limits = c(0, 0.88),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Variance proportion", title = "Official ESTIMATE variance partitioning",
       subtitle = "Sequential increments: stromal 45.4%, immune 31.4%; joint R2 76.7% (rounded); residual 23.3%") +
  panel_theme + theme(legend.position = "bottom", legend.title = element_blank(),
        plot.subtitle = element_text(size = 7.2, hjust = 0),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 8)) +
  annotate("text", x = 0.55, y = 0.88, label = "c", fontface = "bold", size = 4.2,
           hjust = 0, vjust = 1)

row_a <- (p_a + p_b + plot_layout(widths = c(1, 1.05))) +
  plot_annotation(tag_levels = "a")
row_b <- (p_c1 + p_c2 + plot_layout(widths = c(1, 1))) +
  plot_annotation(tag_levels = "a")
row_a <- row_a & theme(plot.tag = element_text(face = "bold", size = 12))
row_b <- row_b & theme(plot.tag = element_text(face = "bold", size = 12))
final_plot <- row_a / row_b / p_d + plot_layout(heights = c(1.0, 1.0, 0.86))

pdf_path <- file.path(out_dir, "SupFig2_current.pdf")
png_path <- file.path(out_dir, "SupFig2_current.png")
svg_path <- file.path(out_dir, "SupFig2_current.svg")
tiff_path <- file.path(out_dir, "SupFig2_current.tiff")
write_vector_pdf <- function(path, plot) {
  cairo_ok <- tryCatch({
    grDevices::cairo_pdf(path, width = 6.69, height = 9.45, family = "Helvetica")
    print(plot)
    grDevices::dev.off()
    TRUE
  }, error = function(e) {
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
    if (file.exists(path)) unlink(path)
    FALSE
  })
  if (!cairo_ok) {
    # This fallback is used on the current macOS runtime because its X11
    # Cairo library is unavailable. The SVG export remains editable/vector.
    grDevices::pdf(path, width = 6.69, height = 9.45, family = "Helvetica")
    print(plot)
    grDevices::dev.off()
  }
}
write_vector_pdf(pdf_path, final_plot)
ggsave(png_path, final_plot, width = 6.69, height = 9.45, dpi = 600, limitsize = FALSE)
svglite::svglite(svg_path, width = 6.69, height = 9.45)
print(final_plot)
dev.off()
ggsave(tiff_path, final_plot, width = 6.69, height = 9.45, dpi = 600, limitsize = FALSE)
writeLines(c(
  "Supplementary Figure 2 current generation",
  sprintf("MGS vs tumour purity Spearman rho = %.10f", unname(rho_test$estimate)),
  sprintf("Purity-adjusted Schwann-cell partial Pearson r = %.10f", partial_r),
  sprintf("Stromal component = %.10f%%", 100 * stromal_r2),
  sprintf("Immune component = %.10f%%", 100 * immune_delta),
  sprintf("Total explained = %.10f%%", 100 * total_r2),
  sprintf("Residual = %.10f%%", 100 * residual),
  "Historical 71.8% proxy decomposition is not plotted."
), file.path(out_dir, "SupFig2_current_generation.txt"))
file.copy(pdf_path, file.path(package_root, "02_Submission_Materials", "02_Supplementary_Figures", "SupFig2.pdf"), overwrite = TRUE)
file.copy(png_path, file.path(package_root, "02_Submission_Materials", "02_Supplementary_Figures", "SupFig2.png"), overwrite = TRUE)
cat("Wrote current Supplementary Figure 2 to:", pdf_path, "\n")
