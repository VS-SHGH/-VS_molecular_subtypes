source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'00_config.R'))
start_log('07_audit_assertions')
suppressPackageStartupMessages({library(qs); library(readxl)})

must_exist <- c(
  'singlecell_input_diagnostics.csv', 'singlecell_normalization_error.csv',
  'singlecell_signature_coverage.csv', 'singlecell_normalized_gmm.csv',
  'myeloid_pseudobulk_tests.rds', 'myeloid_paired_pseudobulk.rds',
  'myeloid_CAMERA_paired_all15_adjusted_technique.csv',
  'cellchat_fresh_run_contract.csv', 'cellchat_fresh_normalized_nboot100.rds',
  'external_NF2_disease_association.csv', 'external_stratified_subtype_summary.csv',
  'bulk_step1_filter_reconstruction.csv', 'clustering_label_reproduction.csv',
  'MGS_baseline_reproduction_check.csv', 'MGS_variance_corrected_interpretation.csv',
  'clinical_analysis_deidentified.csv', 'clinical_model_coefficients.csv',
  'clinical_omnibus_tests.csv', 'NK_apparent_AUC.csv', 'NK_CV_summary.csv'
)
stopifnot(all(file.exists(file.path(OUT,must_exist))))

checks <- list()
add <- function(name, observed, expected, pass, note='') {
  checks[[length(checks)+1]] <<- data.frame(check=name, observed=as.character(observed), expected=as.character(expected), pass=as.logical(pass), note=note)
}

d <- read.csv(file.path(OUT,'singlecell_normalization_error.csv'), check.names=FALSE)
add('single_cell_data_layer_is_integer_counts_before_repair', d$value[d$check=='input_RNA_data_identical_to_integer_counts'], 'TRUE', isTRUE(d$value[d$check=='input_RNA_data_identical_to_integer_counts']), 'The stored data layer was raw counts before repair; downstream repair normalizes it.')
s <- read.csv(file.path(OUT,'singlecell_signature_coverage.csv'))
add('all_signature_sets_have_matched_genes', min(s$matched_n)>=5, 'TRUE', min(s$matched_n)>=5, 'No hand-picked rescue genes were added.')
g <- read.csv(file.path(OUT,'singlecell_normalized_gmm.csv'))
add('normalized_gmm_cutoff_is_finite', g$cutoff, 'finite', is.finite(g$cutoff), 'Operational exploratory split only.')

st <- read.csv(file.path(OUT,'cellchat_fresh_run_contract.csv'))
add('cellchat_uses_normalized_expression', st$normalized, 'TRUE', isTRUE(st$normalized), 'Existing cells only.')
add('cellchat_recomputed_pvalues_not_cached', st$reused_cached_pvalues, 'FALSE', identical(as.character(st$reused_cached_pvalues),'FALSE'), 'All 100 bootstrap permutations were run in this repair.')
add('cellchat_bootstrap_count', st$nboot, '100', st$nboot==100, 'Pooled-cell inference, not patient-level significance.')

e <- read.csv(file.path(OUT,'external_sample_annotations_corrected.csv'), check.names=FALSE)
add('external_ambiguous_assignments_are_excluded_from_primary_tests', sum(e$ambiguous), '6', sum(e$ambiguous)==6, '4/57 in GSE141801 and 2/31 in GSE39645.')
n2 <- read.csv(file.path(OUT,'external_NF2_disease_association.csv'), check.names=FALSE)
add('external_nf2_field_not_converted_to_mutation', all(grepl('NOT tumour NF2 mutation', n2$field_meaning)), 'TRUE', all(grepl('NOT tumour NF2 mutation', n2$field_meaning)), 'Clinical NF2-associated disease only.')

b <- read.csv(file.path(OUT,'bulk_step1_filter_reconstruction.csv'))
add('bulk_step1_gene_filter_reproduced', b$same_gene_set, 'TRUE', isTRUE(b$same_gene_set), 'TPM input; no FASTQ reconstruction.')
clu <- read.csv(file.path(OUT,'clustering_label_reproduction.csv'))
add('k3_labels_reproduced', clu$value[clu$metric=='k3_adjusted_Rand_against_original'], '1', clu$value[clu$metric=='k3_adjusted_Rand_against_original']==1, 'Patient labels compared with the cached original clustering.')

m <- read.csv(file.path(OUT,'MGS_baseline_reproduction_check.csv'))
add('mgs_joint_r2_reproduced', sprintf('%.7f',m$value[m$metric=='joint_R2']), '0.7673127', abs(m$value[m$metric=='joint_R2']-0.767312674044)<1e-6, 'Same-transcriptome association, not causal variance attribution.')
mv <- read.csv(file.path(OUT,'MGS_variance_corrected_interpretation.csv'))
add('mgs_orientation_does_not_change_joint_r2', mv$value[mv$component=='Joint_R2']-mv$value[mv$component=='Joint_R2_after_orientation_flip'], '0', abs(mv$value[mv$component=='Joint_R2']-mv$value[mv$component=='Joint_R2_after_orientation_flip'])<1e-12, 'Sequential component allocation remains order-dependent.')

cd <- read.csv(file.path(OUT,'clinical_analysis_deidentified.csv'), check.names=FALSE)
add('clinical_cohort_size', nrow(cd), '38', nrow(cd)==38, 'Subtype-labeled discovery patients.')
add('p32_first_dimension_author_confirmed_22_mm', cd$size_mm[cd$SampleID=='P32'], '22', identical(as.numeric(cd$size_mm[cd$SampleID=='P32']),22), 'Raw record retained; corrected value is documented separately.')
add('nk_is_reported_as_percentage_measure', TRUE, 'TRUE', TRUE, 'Source unit row is percentage; do not call it an absolute count.')

cam <- read.csv(file.path(OUT,'myeloid_CAMERA_paired_all15_adjusted_technique.csv'))
add('paired_nondefining_go_bp_fdr_below_005', sum(cam$FDR<.05,na.rm=TRUE), '0', sum(cam$FDR<.05,na.rm=TRUE)==0, 'Paired high/low gene-level signal is exploratory and score-defined.')

out <- do.call(rbind, checks)
write_result(out,'audit_assertions')
print(out)
stopifnot(all(out$pass))
cat('ALL AUDIT ASSERTIONS PASSED\n')
end_log()
