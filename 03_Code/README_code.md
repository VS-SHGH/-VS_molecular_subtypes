# Analysis code and numerical evidence

## Two distinct workflows

`verified_numeric_reruns/` contains the downstream workflow previously used to check the existing figures. Its single-cell scripts reproduce the original expression scale and cached CellChat results. Cache agreement alone does not establish correct normalisation.

`sensitivity_checks/` contains the completed methodological reanalysis: log-normalisation, patient-aware myeloid contrasts, fresh-sample CellChat, clinical-model sensitivity, corrected external annotations, and MGS sensitivity. These results constrain the manuscript's interpretation. They are not replacements for the unchanged original-scale Figure 3 panels.

The current editorial review reran script 10 with R 4.2.3 and MCPcounter 1.2.0, compared its values with the earlier computed outputs, and executed scripts 15 and 16 for table reporting. It also re-exported Supplementary Figure 2 for a title-layout correction. All supplied R scripts passed syntax parsing. Syntax checks and selected reruns are not a new end-to-end replication.

## Inputs

Run primary scripts from `verified_numeric_reruns/`, so they can locate `00_config.R`. The configuration defaults to a package-level `input/` folder; environmental overrides are documented in that file. Required inputs vary by script:

- `TPM.xlsx`: 38-sample bulk expression matrix, with `gene_name` and `tpm.<SampleID>` columns.
- `cl_bulk_with_Subtype.xlsx` and `VS_surgery_analysis.xlsx`: clinical and operative variables.
- `RData/`: the original Step1/Step2/Step3 caches and the relevant supplied Seurat/CellChat objects.
- `sc_nc_raw.qs`: the original supplied single-cell object, with counts and annotations.
- `CIBERSORT_Result.csv`: supplied CIBERSORT estimates; script 10 recalculates MCP-counter, not CIBERSORT deconvolution.
- External expression and clinical CSVs for GSE141801 and GSE39645.

These raw or large inputs are not inside the upload archives. The code is not a self-contained download-and-run data bundle. Do not replace an unavailable input with simulated data or treat existing output files as independent source data.

## Execution

After supplying the inputs and compatible packages, the existing primary runner is:

```sh
cd 03_Code/verified_numeric_reruns
Rscript run_verified_downstream.R
```

This runner follows the original downstream dependencies and includes the historical single-cell/cache checks. It does not execute the corrected sensitivity workflow. Figures and table reporting scripts are separate:

```sh
Rscript 12_generate_figure4e_4h_updated.R
Rscript 13_generate_supfig2_current.R
Rscript 15_complete_baseline_table.R
Rscript 16_complete_crossmethod_table.R
```

Script 15 uses the baseline contingency counts and documents the restored source IAC category; original continuous-variable summaries are retained. Script 16 uses the calculated cross-method result CSV, not the earlier rounded submission values. XLSX files are formatted presentations of these CSV tables.

For the new graphical abstract, use Python with matplotlib:

```sh
python3 03_Code/figure_editing/generate_graphical_abstract.py
```

It reads the Figure 4h variance CSV and draws only text, shapes and a numerical bar.

## Software

`analysis_package_versions.csv` is an extract of direct package references matched to the supplied server environment snapshot. It is not a complete dependency lockfile or a statement that every listed version was used in the current rerun. Original analyses reported R 4.4.1; the numerical-verification environment used R 4.2.3. Historical GSVA and Seurat/CellChat code may require compatible package libraries. Use `VS_AUX_R_LIBRARY`, `VS_GSVA_COMPAT_LIBRARY` and the sensitivity configuration's optional library variables rather than overwriting a system library.

The submission-facing files are organised under `02_Submission_Materials/`. The full audit outputs remain under the package-level `outputs/` directory because the R configuration and sensitivity scripts use that location by default.

## Historical outputs

The 71.8% proxy analysis is retained only for provenance. The current manuscript and Figure 4h use the official ESTIMATE component model, joint R² 76.7%. Older output names containing `verified`, `validation` or `PASS` describe specific numerical comparisons, not a blanket judgment on biological validity.
