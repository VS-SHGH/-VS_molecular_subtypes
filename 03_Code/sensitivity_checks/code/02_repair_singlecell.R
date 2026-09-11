source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'00_config.R'))
start_log('02_repair_singlecell')
suppressPackageStartupMessages({library(qs);library(Seurat);library(Matrix);library(edgeR);library(mclust);library(fgsea);library(org.Hs.eg.db);library(AnnotationDbi)})
set.seed(20260906)
sc <- qread(SC_INPUT);ct <- GetAssayData(sc,assay='RNA',layer='counts')
stopifnot(identical(ct,GetAssayData(sc,assay='RNA',layer='data')),all(ct@x>=0),all(ct@x==round(ct@x)))
write_result(data.frame(check='input_RNA_data_identical_to_integer_counts',value=TRUE),'singlecell_normalization_error')
deg <- load_env('Step3_DEA_Results.RData')
top <- function(d){d<-d[is.finite(d$adj.P.Val)&d$adj.P.Val<.05&d$logFC>1,];g<-if('Gene'%in%names(d))d$Gene else rownames(d);unique(g[order(-d$logFC)])[seq_len(min(100,nrow(d)))]}
sets0 <- list(C1=top(deg$deg_c1),C2=top(deg$deg_c2),C3=top(deg$deg_c3))
sets <- lapply(sets0,function(g){i<-match(toupper(g),toupper(rownames(sc)));unique(rownames(sc)[i[!is.na(i)]])})
write_result(data.frame(program=names(sets),definition_n=lengths(sets0),matched_n=lengths(sets)),'singlecell_signature_coverage')
stopifnot(all(lengths(sets)>=5)) # Do not silently add hand-picked rescue genes.
saveRDS(sets,file.path(OUT,'singlecell_program_genes.rds'))
sc <- NormalizeData(sc,assay='RNA',normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
sc <- AddModuleScore(sc,features=sets,name='RepairScore',assay='RNA',seed=42)
for(i in 1:3) sc[[paste0('C',i,'_score')]]<-sc[[paste0('RepairScore',i)]][,1]
md <- sc@meta.data
md$patient<-as.character(md$orig.ident);md$cell_type<-as.character(md$final_label)
md$technique<-factor(md$technique,levels=c('fresh','frozen'))
md$prior_treatment<-ifelse(md$patient=='SCH8','prior_surgery',ifelse(md$patient%in%c('SCH13','SCH14'),'prior_radiation','none_reported'))
sizes<-c(SCH1=2.6,SCH2=3.4,SCH3=2.5,SCH4=3.7,SCH5=2,SCH6=3.1,SCH7=1.3,SCH8=1.4,SCH9=.61,SCH13=3.6,SCH14=2.1,SCH18=1.9,SCH20=2.1,SCH21=3.2,SCH22=3.9)
md$tumour_size_cm<-unname(sizes[md$patient]);stopifnot(!anyNA(md$tumour_size_cm),length(unique(md$patient))==15)
md$log10_UMI<-log10(md$nCount_RNA)
md$C3_exploratory_status<-NA_character_
my <- md$cell_type=='Myeloid'
set.seed(123);gmm<-Mclust(md$C3_score[my],G=2,verbose=FALSE);cut<-mean(gmm$parameters$mean)
md$C3_exploratory_status[my]<-ifelse(md$C3_score[my]>cut,'C3high_Myeloid','C3low_Myeloid')
md$cellchat_label<-md$cell_type;md$cellchat_label[my]<-md$C3_exploratory_status[my]
write_result(data.frame(cutoff=cut,low_mean=min(gmm$parameters$mean),high_mean=max(gmm$parameters$mean),interpretation='exploratory operational split after normalization; not primary patient-level inference'),'singlecell_normalized_gmm')
sc@meta.data<-md
saveRDS(md,file.path(OUT,'singlecell_normalized_metadata.rds'))
qsave(sc,file.path(OUT,'singlecell_normalized_existing_cells.qs'),preset='high')
patcell <- do.call(rbind,lapply(split(seq_len(nrow(md)),interaction(md$patient,md$cell_type,drop=TRUE)),function(i){d<-md[i,];data.frame(patient=d$patient[1],cell_type=d$cell_type[1],technique=as.character(d$technique[1]),prior_treatment=d$prior_treatment[1],tumour_size_cm=d$tumour_size_cm[1],n_cells=nrow(d),C1_score=mean(d$C1_score),C2_score=mean(d$C2_score),C3_score=mean(d$C3_score),median_UMI=median(d$nCount_RNA),median_features=median(d$nFeature_RNA),C3high_fraction=if(d$cell_type[1]=='Myeloid')mean(d$C3_exploratory_status=='C3high_Myeloid')else NA_real_)}))
write_result(patcell,'singlecell_patient_celltype_summary')
pat<-patcell[patcell$cell_type=='Myeloid',];rownames(pat)<-pat$patient
pat$technique<-factor(pat$technique,levels=c('fresh','frozen'));pat$score_z<-as.numeric(scale(pat$C3_score))
write_result(pat,'myeloid_patient_scores')
old<-readRDS(file.path(BASE,'outputs/singlecell_cellchat_rerun/myeloid_scored_metadata.rds'))
stopifnot(all(rownames(old)%in%rownames(md)))
old$patient<-as.character(old$orig.ident)
oldmeans<-aggregate(old$Score_C3,list(patient=old$patient),mean);names(oldmeans)[2]<-'old_raw_count_score'
comparison<-merge(pat,oldmeans,by='patient',sort=FALSE)
write_result(comparison,'myeloid_before_after_patient_scores')
cellcomp<-data.frame(technique=md$technique[my],new_status=md$C3_exploratory_status[my],old_status=ifelse(old$Score_C3[match(rownames(md)[my],rownames(old))]>1.409,'high','low'))
write_result(as.data.frame(table(cellcomp$technique,cellcomp$new_status)),'myeloid_normalized_status_by_technique')
diagnostics<-data.frame(metric=c('old_patient_score_vs_median_UMI_spearman','new_patient_score_vs_median_UMI_spearman','old_new_patient_score_spearman','fresh_frozen_patient_score_wilcoxon_p'),value=c(cor(comparison$old_raw_count_score,comparison$median_UMI,method='spearman'),cor(pat$C3_score,pat$median_UMI,method='spearman'),cor(comparison$old_raw_count_score,comparison$C3_score,method='spearman'),wilcox.test(C3_score~technique,pat,exact=FALSE)$p.value))
write_result(diagnostics,'myeloid_score_diagnostics');print(diagnostics);print(pat)

# Pseudobulk sums counts within patient, never counts cells as independent n.
pm<-factor(md$patient[my],levels=pat$patient)
Z<-sparseMatrix(i=seq_along(pm),j=as.integer(pm),x=1,dims=c(length(pm),nrow(pat)))
pb<-as.matrix(ct[,my,drop=FALSE]%*%Z);colnames(pb)<-pat$patient
saveRDS(list(counts=pb,patient=pat),file.path(OUT,'myeloid_patient_pseudobulk.rds'))
rm(sc,ct);gc()

# Target-gene family excludes the defining C3 genes; associations remain exploratory
# same-transcriptome associations, not independent biological validation.
reslist<-list()
specs<-list(all15_adjusted_technique=pat$patient,fresh11=pat$patient[pat$technique=='fresh'],fresh_no_prior_treatment=pat$patient[pat$technique=='fresh'&pat$prior_treatment=='none_reported'])
for(nm in names(specs)){
  ids<-specs[[nm]];d<-pat[ids,,drop=FALSE];d$score_z<-as.numeric(scale(d$C3_score))
  design<-if(nm=='all15_adjusted_technique')model.matrix(~technique+score_z,d)else model.matrix(~score_z,d)
  stopifnot(qr(design)$rank==ncol(design),length(ids)>=6)
  y<-DGEList(pb[,ids,drop=FALSE]);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y)
  y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);test<-glmQLFTest(fit,coef='score_z')
  tt<-topTags(test,n=Inf,sort.by='none')$table;tt$gene<-rownames(tt);tt$defining_C3_gene<-toupper(tt$gene)%in%toupper(sets$C3)
  tt$FDR_nondefining<-NA_real_;j<-!tt$defining_C3_gene;tt$FDR_nondefining[j]<-p.adjust(tt$PValue[j],method='BH')
  write_result(tt,paste0('myeloid_pseudobulk_',nm));write_result(data.frame(patient=ids,design),paste0('myeloid_design_',nm))
  reslist[[nm]]<-tt
  cat(nm,'patients',length(ids),'genes',nrow(tt),'nondefining FDR<.05',sum(tt$FDR_nondefining<.05,na.rm=TRUE),'\n')
}
saveRDS(reslist,file.path(OUT,'myeloid_pseudobulk_tests.rds'))
common<-Reduce(intersect,lapply(reslist,function(t)t$gene[!t$defining_C3_gene]))
compare<-data.frame(gene=common)
for(nm in names(reslist)){t<-reslist[[nm]];i<-match(common,t$gene);compare[[paste0(nm,'_logFC')]]<-t$logFC[i];compare[[paste0(nm,'_FDR')]]<-t$FDR_nondefining[i]}
write_result(compare,'myeloid_sensitivity_gene_comparison')

cat('Constructing GO BP gene sets from installed annotation\n')
symbols<-intersect(unique(unlist(lapply(reslist,function(x)x$gene))),keys(org.Hs.eg.db,keytype='SYMBOL'))
an<-AnnotationDbi::select(org.Hs.eg.db,keys=symbols,columns=c('GOALL','ONTOLOGYALL'),keytype='SYMBOL')
an<-an[!is.na(an$GOALL)&an$ONTOLOGYALL=='BP',];pathways<-lapply(split(an$SYMBOL,an$GOALL),unique)
for(nm in names(reslist)){
  t<-reslist[[nm]];t<-t[!t$defining_C3_gene,];r<-setNames(sign(t$logFC)*sqrt(t$F),t$gene);r<-sort(r[is.finite(r)],decreasing=TRUE)
  set.seed(20260906);g<-fgseaMultilevel(pathways=pathways,stats=r,minSize=15,maxSize=500,eps=0,nproc=1)
  g<-as.data.frame(g);g$leadingEdge<-vapply(g$leadingEdge,paste,collapse=';',FUN.VALUE='')
  if(requireNamespace('GO.db',quietly=TRUE)){gn<-AnnotationDbi::select(GO.db::GO.db,keys=g$pathway,keytype='GOID',columns='TERM');g$description<-gn$TERM[match(g$pathway,gn$GOID)]}
  write_result(g,paste0('myeloid_GO_BP_',nm))
}
end_log()
