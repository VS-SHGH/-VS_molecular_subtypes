# =============================================================================
# Revision Analysis 1: NK-cell ROC + Regression Diagnostics
# =============================================================================

library(dplyr)
library(ggplot2)
library(pROC)
library(car)
library(lmtest)

# ---- Load data ----
cl <- readxl::read_excel("/Users/yanchen/Desktop/VSbulk+单细胞/cl_bulk_with_Subtype.xlsx")
# Skip row 0 (reference ranges), keep rows with actual data
cl <- cl[!is.na(cl$Subtype) & cl$Subtype != "", ]
cl <- cl[!is.na(cl$Age), ]

cat("Samples with subtype:", nrow(cl), "\n")
cat("Subtype distribution:\n")
print(table(cl$Subtype))

# ---- Extract variables ----
cl$NK_cells <- as.numeric(cl$自然杀伤细胞)
cl$Age <- as.numeric(cl$Age)
cl$Subtype <- factor(cl$Subtype, levels = c("C1", "C2", "C3"))
cl$C3_vs_other <- ifelse(cl$Subtype == "C3", "C3", "C1_C2")

# ---- 1. NK-cell ROC Analysis ----
cat("\n========== NK-cell ROC Analysis ==========\n")

nk_data <- cl[!is.na(cl$NK_cells) & !is.na(cl$C3_vs_other), ]
cat("Samples with NK data:", nrow(nk_data), "\n")

# ROC
roc_obj <- roc(nk_data$C3_vs_other, nk_data$NK_cells, levels = c("C1_C2", "C3"))
cat(sprintf("ROC AUC: %.3f (95%% CI: %.3f-%.3f)\n",
            auc(roc_obj), ci.auc(roc_obj)[1], ci.auc(roc_obj)[3]))

# Optimal cutoff (Youden index)
coords_best <- coords(roc_obj, "best", best.method = "youden")
cat(sprintf("Best cutoff: NK = %.2f, Sensitivity: %.2f, Specificity: %.2f\n",
            coords_best$threshold, coords_best$sensitivity, coords_best$specificity))

# Cutoff for 80% sensitivity
coords_80 <- coords(roc_obj, x = 0.8, input = "sensitivity")
cat(sprintf("Cutoff for 80%% sensitivity: NK = %.2f, Specificity: %.2f\n",
            coords_80$threshold[1], coords_80$specificity[1]))

# PPV/NPV at best cutoff
nk_data$pred_C3 <- ifelse(nk_data$NK_cells >= coords_best$threshold, "C3", "C1_C2")
tp <- sum(nk_data$pred_C3 == "C3" & nk_data$C3_vs_other == "C3")
fp <- sum(nk_data$pred_C3 == "C3" & nk_data$C3_vs_other == "C1_C2")
tn <- sum(nk_data$pred_C3 == "C1_C2" & nk_data$C3_vs_other == "C1_C2")
fn <- sum(nk_data$pred_C3 == "C1_C2" & nk_data$C3_vs_other == "C3")
ppv <- tp / (tp + fp)
npv <- tn / (tn + fn)
cat(sprintf("PPV: %.2f, NPV: %.2f\n", ppv, npv))

# Save ROC plot
pdf("/Users/yanchen/Desktop/VS_revision/output/NK_ROC_C3_vs_C1C2.pdf", width = 5, height = 5)
plot(roc_obj, main = sprintf("NK cells: C3 vs C1/C2 (AUC=%.2f)", auc(roc_obj)))
dev.off()

# ---- 2. Regression Diagnostics ----
cat("\n========== Regression Diagnostics ==========\n")

# Age-adjusted model
m1 <- lm(NK_cells ~ Subtype + Age, data = cl)
cat("\n--- Model 1: NK ~ Subtype + Age ---\n")
print(summary(m1))

# Age + tumour size model
cl$SizeLarge <- ifelse(cl$SizeGrade == "large", 1, 0)
m2 <- lm(NK_cells ~ Subtype + Age + SizeLarge, data = cl)
cat("\n--- Model 2: NK ~ Subtype + Age + SizeLarge ---\n")
print(summary(m2))
cat(sprintf("EPV (Model 2): %.1f\n", nrow(cl) / length(coef(m2))))

# Residual diagnostics
cat("\n--- Shapiro-Wilk test for normality of residuals ---\n")
cat(sprintf("Model 1: W=%.3f, P=%.4f\n",
            shapiro.test(residuals(m1))$statistic, shapiro.test(residuals(m1))$p.value))
cat(sprintf("Model 2: W=%.3f, P=%.4f\n",
            shapiro.test(residuals(m2))$statistic, shapiro.test(residuals(m2))$p.value))

cat("\n--- Breusch-Pagan test for homoscedasticity ---\n")
cat(sprintf("Model 1: BP=%.2f, P=%.4f\n",
            bptest(m1)$statistic, bptest(m1)$p.value))
cat(sprintf("Model 2: BP=%.2f, P=%.4f\n",
            bptest(m2)$statistic, bptest(m2)$p.value))

# Cook's distance
cooks1 <- cooks.distance(m1)
cooks2 <- cooks.distance(m2)
cat(sprintf("\nModel 1: max Cook's D = %.4f (threshold = %.4f)\n",
            max(cooks1), 4/nrow(cl)))
cat(sprintf("Model 2: max Cook's D = %.4f (threshold = %.4f)\n",
            max(cooks2), 4/nrow(cl)))
cat(sprintf("Influential obs (Cook's D > 4/n): Model 1: %d, Model 2: %d\n",
            sum(cooks1 > 4/nrow(cl)), sum(cooks2 > 4/nrow(cl))))

# VIF for Model 2
cat("\n--- VIF (Model 2) ---\n")
print(vif(m2))

# ---- Diagnostic plots ----
pdf("/Users/yanchen/Desktop/VS_revision/output/Regression_Diagnostics.pdf", width = 10, height = 8)
par(mfrow = c(2, 3))
# Model 1
plot(fitted(m1), residuals(m1), main = "Model 1: Residuals vs Fitted",
     xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
qqnorm(residuals(m1), main = "Model 1: Q-Q Plot")
qqline(residuals(m1), col = "red")
plot(cooks1, type = "h", main = "Model 1: Cook's Distance",
     ylab = "Cook's D", xlab = "Observation")
abline(h = 4/nrow(cl), col = "red", lty = 2)

# Model 2
plot(fitted(m2), residuals(m2), main = "Model 2: Residuals vs Fitted",
     xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
qqnorm(residuals(m2), main = "Model 2: Q-Q Plot")
qqline(residuals(m2), col = "red")
plot(cooks2, type = "h", main = "Model 2: Cook's Distance",
     ylab = "Cook's D", xlab = "Observation")
abline(h = 4/nrow(cl), col = "red", lty = 2)
dev.off()

# ---- 3. NK-cell distribution by subtype (violin plot) ----
pdf("/Users/yanchen/Desktop/VS_revision/output/NK_by_Subtype.pdf", width = 6, height = 5)
ggplot(cl[!is.na(cl$NK_cells), ], aes(x = Subtype, y = NK_cells, fill = Subtype)) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
  labs(y = "NK Cells (%)", title = "Circulating NK Cells by Molecular Subtype") +
  theme_minimal(base_size = 13)
dev.off()

cat("\n========== Analysis Complete ==========\n")
