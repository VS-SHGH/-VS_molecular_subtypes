source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'00_config.R'))
start_log('02b_patient_paired_and_pathways')
suppressPackageStartupMessages({library(qs);library(SeuratObject);library(Matrix);library(edgeR);library(limma);library(org.Hs.eg.db);library(AnnotationDbi);library(GO.db)})
sets<-readRDS(file.path(OUT,'singlecell_program_genes.rds'));base<-readRDS(file.path(OUT,'myeloid_patient_pseudobulk.rds'))
pb<-base$counts;pat<-base$patient
symbols<-intersect(rownames(pb),keys(org.Hs.eg.db,keytype='SYMBOL'))
an<-AnnotationDbi::select(org.Hs.eg.db,keys=symbols,columns=c('GOALL','ONTOLOGYALL'),keytype='SYMBOL');an<-an[!is.na(an$GOALL)&an$ONTOLOGYALL=='BP',]
pathways<-lapply(split(an$SYMBOL,an$GOALL),unique);saveRDS(pathways,file.path(OUT,'GO_BP_gene_sets_audit.rds'))
camera_test<-function(counts,design,coef,name){
 y<-DGEList(counts);y<-y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE];y<-calcNormFactors(y)
 v<-voom(y,design=design,plot=FALSE);keep<-!toupper(rownames(v))%in%toupper(sets$C3);v<-v[keep,]
 inds<-lapply(pathways,function(g)which(rownames(v)%in%g));inds<-inds[lengths(inds)>=15&lengths(inds)<=500]
 ca<-camera(v,index=inds,design=design,contrast=coef,inter.gene.cor=NA);ca$GO_ID<-rownames(ca)
 terms<-AnnotationDbi::select(GO.db,keys=ca$GO_ID,keytype='GOID',columns='TERM');ca$description<-terms$TERM[match(ca$GO_ID,terms$GOID)]
 write_result(ca,paste0('myeloid_CAMERA_',name));cat(name,'CAMERA FDR<.05:',sum(ca$FDR<.05),'/',nrow(ca),'\n')
}
specs<-list(all15_adjusted_technique=pat$patient,fresh11=pat$patient[pat$technique=='fresh'],fresh_no_prior_treatment=pat$patient[pat$technique=='fresh'&pat$prior_treatment=='none_reported'])
for(nm in names(specs)){
 d<-pat[specs[[nm]],];d$score_z<-as.numeric(scale(d$C3_score));design<-if(nm=='all15_adjusted_technique')model.matrix(~technique+score_z,d)else model.matrix(~score_z,d)
 camera_test(pb[,d$patient,drop=FALSE],design,which(colnames(design)=='score_z'),nm)
}
# Explicit depth-adjusted sensitivity, not chosen according to significance.
d<-pat;d$score_z<-as.numeric(scale(d$C3_score));d$depth_z<-as.numeric(scale(log10(d$median_UMI)));design<-model.matrix(~technique+depth_z+score_z,d)
y<-DGEList(pb[,d$patient]);y<-y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE)
fit<-glmQLFit(y,design,robust=TRUE);tt<-topTags(glmQLFTest(fit,coef='score_z'),n=Inf)$table;tt$gene<-rownames(tt);tt$defining_C3_gene<-toupper(tt$gene)%in%toupper(sets$C3);j<-!tt$defining_C3_gene;tt$FDR_nondefining<-NA_real_;tt$FDR_nondefining[j]<-p.adjust(tt$PValue[j],'BH');write_result(tt,'myeloid_depth_adjusted_sensitivity')

# Paired pseudobulk answers the original high-vs-low WITHIN-patient question,
# distinct from the primary BETWEEN-patient continuous-score association.
sc<-qread(file.path(OUT,'singlecell_normalized_existing_cells.qs'));md<-sc@meta.data
my<-md$cell_type=='Myeloid';md<-md[my,];ct<-GetAssayData(sc,layer='counts')[,rownames(md),drop=FALSE]
group<-interaction(md$patient,md$C3_exploratory_status,drop=TRUE,sep='|');Z<-sparseMatrix(i=seq_along(group),j=as.integer(group),x=1,dims=c(length(group),nlevels(group)))
paired<-as.matrix(ct%*%Z);colnames(paired)<-levels(group)
info<-md[match(levels(group),as.character(group)),c('patient','technique','C3_exploratory_status','prior_treatment')]
info$group<-levels(group);info$n_cells<-as.numeric(table(group));info$state<-factor(ifelse(info$C3_exploratory_status=='C3high_Myeloid','high','low'),levels=c('low','high'))
saveRDS(list(counts=paired,metadata=info),file.path(OUT,'myeloid_paired_pseudobulk.rds'));rm(sc,ct);gc()
for(nm in names(specs)){
 ii<-info$patient%in%specs[[nm]]&info$n_cells>=10;d<-info[ii,];complete_patients<-names(which(table(d$patient)==2));d<-d[d$patient%in%complete_patients,];d$patient<-factor(d$patient)
 design<-model.matrix(~patient+state,d);stopifnot(qr(design)$rank==ncol(design));counts<-paired[,d$group,drop=FALSE]
 y<-DGEList(counts);y<-y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 t<-topTags(glmQLFTest(fit,coef='statehigh'),n=Inf)$table;t$gene<-rownames(t);t$defining_C3_gene<-toupper(t$gene)%in%toupper(sets$C3);j<-!t$defining_C3_gene;t$FDR_nondefining<-NA_real_;t$FDR_nondefining[j]<-p.adjust(t$PValue[j],'BH')
 write_result(t,paste0('myeloid_paired_DEG_',nm));write_result(data.frame(d,design),paste0('myeloid_paired_design_',nm));camera_test(counts,design,which(colnames(design)=='statehigh'),paste0('paired_',nm))
 cat(nm,'paired patients',length(complete_patients),'nondefining FDR<.05',sum(t$FDR_nondefining<.05,na.rm=TRUE),'\n')
}
end_log()
