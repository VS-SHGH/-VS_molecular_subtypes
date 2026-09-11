# Complete categorical tests from the published baseline contingency counts.
# Existing means, standard deviations and continuous-variable tests are retained.
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep('^--file=', commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub('^--file=', '', script_arg[1]) else getwd()
root <- normalizePath(file.path(dirname(script_path), '..', '..'), mustWork = TRUE)
table_path <- if (length(args)) args[1] else file.path(root, '02_Submission_Materials', '03_Supplementary_Tables', 'Supplementary_Table1_Baseline.csv')
tab <- read.csv(table_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(tab) %in% c(22L, 23L), ncol(tab) == 5L)
# The source clinical workbook records one additional C2 IAC category.
# Retain it as recorded rather than merging it into a different category.
postoperative_label <- 'IAC involvement (Expanded, postoperative change), n (%)'
if (!postoperative_label %in% tab[[1]]) {
  added <- tab[12, , drop = FALSE]
  added[1, ] <- c(postoperative_label, '0 (0.0%)', '1 (8.3%)', '0 (0.0%)', '')
  tab <- rbind(tab[1:12, ], added, tab[13:nrow(tab), ])
}

labels <- c('Tumour consistency (Solid), n (%)' = 'Tumour morphology (Solid), n (%)',
            'Tumour consistency (Cystic), n (%)' = 'Tumour morphology (Cystic), n (%)',
            'NK cell count (%), mean ± SD' = 'NK-cell percentage (%), mean ± SD',
            'Total T lymphocyte count (%), mean ± SD' = 'Total T-cell percentage (%), mean ± SD',
            'Total B lymphocyte count (%), mean ± SD' = 'Total B-cell percentage (%), mean ± SD')
for (old in names(labels)) tab[[1]][tab[[1]] == old] <- labels[[old]]

groups <- list(Sex = 2L, Koos_grade = 3:5, Tumour_size = 6:7,
               Tumour_morphology = 8:9, IAC_involvement = 10:13, Recurrence = 14:15)
results <- list()
for (group_name in names(groups)) {
  rr <- groups[[group_name]]
  if (group_name == 'Sex') {
    male <- as.numeric(sub('M([0-9]+).*', '\\1', as.character(tab[rr, 2:4])))
    female <- as.numeric(sub('.*F([0-9]+).*', '\\1', as.character(tab[rr, 2:4])))
    counts <- rbind(male, female)
  } else {
    counts <- apply(as.matrix(tab[rr, 2:4]), 2, function(x) as.numeric(sub(' .*', '', x)))
  }
  stopifnot(identical(as.numeric(colSums(counts)), c(13, 12, 13)))
  p <- fisher.test(counts)$p.value
  tab[rr, 5] <- ''
  tab[rr[1], 5] <- sprintf('%.3f', p)
  results[[group_name]] <- data.frame(Variable = group_name, P_value = p,
                                     Test = 'Fisher exact test', N = sum(counts))
}
write.csv(tab, table_path, row.names = FALSE, na = '', fileEncoding = 'UTF-8')
out_dir <- file.path(root, 'outputs', 'baseline_table')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, results), file.path(out_dir, 'categorical_fisher_tests.csv'), row.names = FALSE)
cat('Completed baseline table: existing descriptive values retained.\n')
print(do.call(rbind, results))
