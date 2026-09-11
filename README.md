# VS molecular subtypes — current analysis code

This repository contains the current downstream numerical-verification and sensitivity-analysis code associated with the Discover Oncology submission package.

## Repository contents

- `03_Code/verified_numeric_reruns/`: downstream reruns for bulk clustering, MGS, external validation, single-cell/CellChat numerical checks, clinical models, tables and updated figure panels.
- `03_Code/sensitivity_checks/`: methodological sensitivity analyses, including patient-aware myeloid contrasts, fresh-sample CellChat, corrected external annotations and MGS sensitivity checks.
- `03_Code/figure_editing/`: graphical-abstract generation script.
- `03_Code/analysis_package_versions.csv`: package-version snapshot used for the audit.

The sensitivity results are included as CSV evidence under `03_Code/sensitivity_checks/results/`. The repository does not include raw FASTQ, raw 10X files, private clinical source files or large Seurat/CellChat objects.

## Inputs and reproducibility

The scripts are not a self-contained download-and-run data bundle. Required inputs vary by workflow and include the bulk expression matrix, clinical tables, supplied RData/Seurat/CellChat objects, CIBERSORT estimates and external-cohort files. Place them under a local `input/` directory or set the `VS_*` environment variables documented in `03_Code/verified_numeric_reruns/00_config.R`.

After supplying compatible inputs and R packages, the primary downstream runner is:

```sh
cd 03_Code/verified_numeric_reruns
Rscript run_verified_downstream.R
```

Selected reporting and figure-generation scripts are run separately:

```sh
Rscript 12_generate_figure4e_4h_updated.R
Rscript 13_generate_supfig2_current.R
Rscript 15_complete_baseline_table.R
Rscript 16_complete_crossmethod_table.R
python3 ../figure_editing/generate_graphical_abstract.py
```

See `03_Code/README_code.md` and `03_Code/sensitivity_checks/README.md` for workflow details. Syntax checks and selected reruns do not constitute a new end-to-end replication when required source inputs are unavailable.

## Software notes

`analysis_package_versions.csv` is an extract of direct package references matched to the supplied server-environment snapshot; it is not a complete lockfile. Original analyses reported R 4.4.1, whereas the numerical-verification environment used R 4.2.3. Historical GSVA and Seurat/CellChat code may require compatible package libraries.

The current manuscript and Figure 4h use the official ESTIMATE component model; the historical 71.8% proxy analysis is retained only for provenance in the accompanying audit package.
