# VS Molecular Subtyping — Analysis Code

## Environment
R 4.2.3 | Packages: ConsensusClusterPlus, limma, GSVA, ESTIMATE, MCPcounter, CellChat, Seurat, princurve, pROC, clusterProfiler, ggplot2, dplyr, pheatmap

---

## Code by Manuscript Section

### 01_Bulk_Clustering (Figure 1)
| Script | What it does |
|--------|-------------|
| `01_consensus_clustering_subtype_discovery.R` | Consensus clustering (k=3), PCA, DEGs, GO/KEGG enrichment, Hallmark validation |
| `02_clustering_sensitivity_PAC.R` | PAC scores across gene counts (1,000-10,000), k=2-6 |
| `03_nf2_mutation_external_cohorts.R` | NF2 mutation analysis in GSE141801, GSE39645 |

### 02_Immune_Microenvironment (Figure 2)
| Script | What it does |
|--------|-------------|
| `01_ESTIMATE_immune_signatures.R` | ESTIMATE, TME heatmaps, checkpoint genes, immune recruitment scores |
| `02_comprehensive_immune_scoring.R` | ESTIMATE + custom signatures, MCPcounter, CIBERSORT, correlation matrices |
| `03_MCPcounter_deconvolution.R` | MCP-counter deconvolution (10 cell types × 38 samples) |
| `04_CIBERSORT_MCP_cross_validation.R` | Cross-method comparison scatter plots, Spearman correlations |

### 03_SingleCell_CellChat (Figure 3)
| Script | What it does |
|--------|-------------|
| `01_single_cell_analysis.R` | UMAP, module scoring, GMM (C3-high vs C3-low), CellChat communication |

### 04_Gradient_MGS (Figure 4)
| Script | What it does |
|--------|-------------|
| `01_MGS_construction_original.R` | PCA + principal curve MGS (original subtype-anchored version) |
| `02_MGS_unbiased_immune_anchoring.R` | **Immune-anchored MGS** — eliminates circularity, PC1 loading analysis |
| `03_MGS_external_validation_projection.R` | MGS projected to GSE141801, GSE39645 |
| `04_MGS_independent_per_cohort.R` | Independent PCA + principal curve in each cohort |
| `05_variance_partitioning_MGS.R` | Hierarchical variance partitioning (71.8% composition) |
| `06_schwann_cell_differentiation_analysis.R` | SC marker expression by subtype, compositional vs cell-intrinsic |
| `07_MGS_sensitivity_bootstrap.R` | Permutation, bootstrap, LOCO validation |

### 05_Clinical_Phenotypes (Figure 5)
| Script | What it does |
|--------|-------------|
| `01_surgical_phenotypes_perioperative.R` | Alluvial diagram, texture/adhesion analysis, perioperative heatmaps |
| `02_NK_cell_ROC_regression.R` | NK-cell ROC (AUC=0.874), regression diagnostics, Shapiro-Wilk, Cook's D |
| `03_clinical_outcomes_HB_EOR.R` | Subtype vs HB grade + extent of resection (25 matched samples) |

### 06_Revision_Analyses (Supplementary)
| Script | What it does |
|--------|-------------|
| `01_generate_new_figure_panels.R` | Generate new panels: Fig 1i, 2h, 4a-inset, 4k, 4l, Sup Figs 6/7/10/11 |
| `02_verify_all_results_reproducibility.R` | Re-run all key analyses, verify manuscript values (21/25 PASS) |

### 07_Data_Preparation
| Script | What it does |
|--------|-------------|
| `01_download_preprocess_external_cohorts.R` | Download GSE141801, GSE39645; probe-to-symbol; log2 normalization |
| `02_download_preprocess_GSE216781.R` | Download GSE216781 RNA-seq counts; TPM normalization |

---

## Key Metrics (Verified)

| Metric | Value |
|--------|-------|
| MGS C1 / C2 / C3 | 0.13 / 0.47 / 0.82 (ρ=0.929, P=1.16×10⁻⁷) |
| NK AUC (C3 vs C1/C2) | 0.874 (95% CI: 0.763-0.985) |
| C3-NK full model | β=3.251, P=0.053 (tumor size β=5.696, P=0.0003) |
| MGS variance (composition) | 71.8% (purity 27.2%, immune 41.2%, fibroblast 1.4%, other 2.0%) |
| MCP Fibroblasts C2 | 444 vs C1=199 |
| MCP Monocytes C3 | 178 vs C1=86, P<0.0001 |
| CIBERSORT-MCP CD8T ρ | 0.833, P<0.001 |
| NF2 C1 enrichment | 41.2% vs 7.4-11.1% (observed 21.6% vs expected ~65%) |

## Data Files
- `TPM.xlsx` — 31,626 genes × 46 samples WITH gene_name
- `cl_bulk_with_Subtype.xlsx` — 38-sample clinical data
- `CIBERSORT_Result.csv` — CIBERSORT fractions (in VS.4/supplement/all/step8_Immune/)
- GSE141801, GSE39645, GSE216784 — GEO public datasets
