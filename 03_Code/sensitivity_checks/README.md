# Supplementary Data 1: methodological sensitivity analyses

These scripts and 72 numerical CSV files come from the completed 6–7 September 2026 methodological reanalysis, without adding patients. Numerical CSVs are copied without alteration. Input-location and optional-library defaults were made portable; analytical procedures were not changed in this copy.

The code covers the original 15 public single-cell samples, subsets of 11 fresh samples and 8 untreated fresh samples, the 38-patient bulk discovery cohort, and the same 88 external bulk samples. Patient identifiers in output tables are study codes, not names. The data-sharing scope remains subject to the authors' ethics and consent requirements.

## Dependencies and order

Place this directory at `03_Code/sensitivity_checks/` within the analysis package or set `VS_REPAIR_ROOT` and `VS_BASELINE_ROOT` explicitly. The baseline root must provide the primary analysis outputs referenced by the scripts. Large QS/RDS inputs and newly generated RDS intermediates are not included here.

`code/00_config.R` defines input and output variables. In particular, `VS_SERVER_ROOT` contains the `RData/` folder; `VS_SC_QS`, `VS_TPM_XLSX` and `VS_EXTERNAL_ROOT` point to the single-cell, bulk and external inputs. The external root contains `GSE141801_Gugel/` and `GSE39645_Torres/`, each with its expression CSV and `Clinical_Data.csv`.

The clinical script reads `cl_bulk_with_Subtype.xlsx` and `VS 手术分析.xlsx` directly under `VS_SERVER_ROOT`. This is the original workbook filename used by the sensitivity workflow, distinct from the primary workflow's configurable filename. Public numerical tables are not a replacement for these original clinical workbooks.

Run scripts as separate R processes, in this order:

```sh
Rscript code/01_inspect_singlecell.R
Rscript code/02_repair_singlecell.R
Rscript code/02b_patient_paired_and_pathways.R
Rscript code/03_repair_clinical.R
Rscript code/04_repair_external.R
Rscript code/05_repair_cellchat_fresh.R
Rscript code/06_repair_mgs_clustering.R
Rscript code/07_audit_assertions.R
```

`00_install_missing.R` records historical dependency preparation; it is not run automatically. Review the requested installations and package compatibility before using it. The old GSVA API may require a separate library selected with `VS_GSVA_COMPAT_LIBRARY`.

The canonical results are the original completed CSV outputs, not this editorial review's syntax check. The key limitations are normalisation-sensitive scoring, no retained patient-aware pathway FDR support after excluding defining genes, and the absence of a C3-high-specific whole-network hub in the fresh-sample analysis. These findings are reported rather than hidden by the retained Figure 3 graphics.
