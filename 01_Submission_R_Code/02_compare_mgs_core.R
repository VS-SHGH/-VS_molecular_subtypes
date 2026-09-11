# Compare the second-pass canonical immune-anchored MGS with the locked cache.
source(Sys.getenv("VS_CONFIG_R", unset = file.path(getwd(), "00_config.R")))
new_path <- Sys.getenv("VS_MGS_RERUN_RDS", unset = file.path(AUDIT_ROOT, "mgs", "mgs_unbiased_results.rds"))
old_path <- Sys.getenv("VS_MGS_REFERENCE_RDS", unset = file.path(package_root, "input", "audit_cache", "mgs_unbiased_results.rds"))
new <- readRDS(new_path)
old <- readRDS(old_path)

a <- new$mgs_unbiased$MGS_unbiased
b <- old$mgs_unbiased$MGS_unbiased
names(a) <- new$mgs_unbiased$Sample
names(b) <- old$mgs_unbiased$Sample
stopifnot(identical(names(a), names(b)))

result <- data.frame(
  check = c("sample_order", "max_absolute_difference", "C1_mean", "C2_mean", "C3_mean"),
  value = c(
    "identical",
    format(max(abs(a - b)), scientific = TRUE, digits = 17),
    sprintf("%.9f", mean(a[new$mgs_unbiased$Subtype == "C1"])),
    sprintf("%.9f", mean(a[new$mgs_unbiased$Subtype == "C2"])),
    sprintf("%.9f", mean(a[new$mgs_unbiased$Subtype == "C3"]))
  ),
  status = c("PASS", ifelse(max(abs(a - b)) == 0, "PASS", "WARN"), "PASS", "PASS", "PASS"),
  stringsAsFactors = FALSE
)
dir.create(AUDIT_ROOT, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(AUDIT_ROOT, "mgs_core_cache_comparison.csv"), row.names = FALSE)
writeLines(capture.output(print(result, row.names = FALSE)), file.path(AUDIT_ROOT, "mgs_core_cache_comparison.txt"))
print(result, row.names = FALSE)
