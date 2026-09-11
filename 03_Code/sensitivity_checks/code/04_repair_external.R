source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'00_config.R'))
start_log('04_repair_external')
gsva_compat <- Sys.getenv('VS_GSVA_COMPAT_LIBRARY',unset='')
if(nzchar(gsva_compat)) .libPaths(c(gsva_compat,.libPaths()))
suppressPackageStartupMessages({library(GSVA);library(limma)})
gs<-load_env('Step3_DEA_Results.RData')$gene_sets
readexpr<-function(f){d<-read.csv(f,check.names=FALSE);m<-as.matrix(d[,-1]);storage.mode(m)<-'numeric';rownames(m)<-as.character(d[[1]]);stopifnot(!anyNA(m),!anyDuplicated(rownames(m)),!anyDuplicated(colnames(m)));m}
folders<-c(GSE141801='GSE141801_Gugel',GSE39645='GSE39645_Torres')
expr<-lapply(names(folders),function(n)readexpr(file.path(EXTERNAL_ROOT,folders[[n]],paste0(n,'_Expression_Log2.csv'))));names(expr)<-names(folders)
call_scores<-function(z){pred<-max.col(t(z),ties.method='first');margin<-apply(z,2,function(x){r<-sort(x,decreasing=TRUE);r[1]-r[2]});data.frame(SampleID=colnames(z),top_label=rownames(z)[pred],delta=margin,ambiguous=margin<.1,assignment=ifelse(margin<.1,'Ambiguous',rownames(z)[pred]))}
all<-list();tests<-list();nft<-list();summary<-list()
for(n in names(expr)){
 e<-expr[[n]];clin<-read.csv(file.path(EXTERNAL_ROOT,folders[[n]],'Clinical_Data.csv'),check.names=FALSE)
 stopifnot(!anyDuplicated(clin$ID),all(colnames(e)%in%clin$ID));clin<-clin[match(colnames(e),clin$ID),]
 s<-gsva(e,gs,method='ssgsea',kcdf='Gaussian',verbose=FALSE);z<-t(scale(t(s)));d<-call_scores(z)
 d$Dataset<-n;d$NF2_disease_status<-ifelse(clin$NF2=='Y','NF2-associated',ifelse(clin$NF2=='N','sporadic','unknown'))
 d$prior_irradiation<-ifelse(clin$Prior_SRS=='Y','yes',ifelse(clin$Prior_SRS=='N','no','unknown'))
 d$prior_surgery<-ifelse(is.na(clin$Prior_Sx),'unknown',ifelse(clin$Prior_Sx=='Y','yes','no'))
 d$NF2_sequence_mutation<-NA_character_;d$chr22q_LOH<-NA_character_;d$NF2_biallelic_hits<-NA_character_
 d$primary_nonirradiated_sporadic<-d$NF2_disease_status=='sporadic'&d$prior_irradiation=='no'
 # NA molecular annotations are intentional: no unverified disease-to-mutation conversion.
 write_result(data.frame(SampleID=colnames(s),t(s)),paste0(n,'_ssgsea_raw_scores'))
 write_result(data.frame(SampleID=colnames(z),t(z)),paste0(n,'_ssgsea_cohort_z_scores'))
 groups<-list(all_original=rep(TRUE,nrow(d)),nonirradiated_sporadic=d$primary_nonirradiated_sporadic,nonirradiated_NF2=d$NF2_disease_status=='NF2-associated'&d$prior_irradiation=='no',irradiated=d$prior_irradiation=='yes')
 for(g in names(groups)){
  ii<-which(groups[[g]]);dd<-d[ii,,drop=FALSE]
  summary[[paste(n,g)]]<-data.frame(Dataset=n,stratum=g,n=nrow(dd),ambiguous=sum(dd$ambiguous),C1=sum(grepl('^C1',dd$assignment)),C2=sum(grepl('^C2',dd$assignment)),C3=sum(grepl('^C3',dd$assignment)))
  ii<-ii[!d$ambiguous[ii]]
  for(gene in c('MKI67','CD8A','LAG3','HAVCR2')){
   p<-if(length(ii)>=6&&length(unique(d$assignment[ii]))>=2&&gene%in%rownames(e))kruskal.test(e[gene,ii]~factor(d$assignment[ii]))$p.value else NA_real_
   tests[[length(tests)+1]]<-data.frame(Dataset=n,stratum=g,gene=gene,n=length(ii),p=p,used_in_signature=gene%in%unlist(gs),interpretation='signature consistency, not an independent target-label validation')
  }
 }
 for(scope in c('all_unambiguous','nonirradiated_unambiguous')){
  ii<-!d$ambiguous&d$NF2_disease_status!='unknown';if(scope=='nonirradiated_unambiguous')ii<-ii&d$prior_irradiation=='no'
  tab<-table(d$NF2_disease_status[ii],d$assignment[ii]);p<-if(all(dim(tab)>1))fisher.test(tab)$p.value else NA_real_
  nft[[length(nft)+1]]<-data.frame(Dataset=n,scope=scope,n=sum(ii),NF2_associated=sum(d$NF2_disease_status[ii]=='NF2-associated'),p=p,field_meaning='clinical NF2-associated disease; NOT tumour NF2 mutation')
 }
 all[[n]]<-d
}
d<-do.call(rbind,all);st<-do.call(rbind,summary);tt<-do.call(rbind,tests);tt$q_within_dataset_stratum<-ave(tt$p,interaction(tt$Dataset,tt$stratum),FUN=function(p)p.adjust(p,'BH'))
write_result(d,'external_sample_annotations_corrected');write_result(st,'external_stratified_subtype_summary');write_result(tt,'external_marker_tests_corrected');write_result(do.call(rbind,nft),'external_NF2_disease_association')
print(st);print(do.call(rbind,nft))

# Additional transfer sensitivity using ONLY the same existing 88 patients.
# Define a common measurable gene universe without using external outcomes;
# freeze discovery score means/SDs. No external-cohort re-centring or re-fitting.
# This rank-score sensitivity is distinct from the historical ssGSEA method.
disc<-load_env('Step2_Clustering.RData');train<-disc$log_tpm
u<-Reduce(intersect,c(list(rownames(train)),lapply(expr,rownames)))
rankscore<-function(m){r<-apply(m[u,,drop=FALSE],2,rank,ties.method='average')/length(u);s<-vapply(gs,function(g)colMeans(r[intersect(g,u),,drop=FALSE]),numeric(ncol(r)));t(s)}
ts<-rankscore(train);rownames(ts)<-names(gs);colnames(ts)<-colnames(train);mu<-rowMeans(ts);sig<-apply(ts,1,sd)
stopifnot(all(sig>0));write_result(data.frame(signature=names(gs),matched_genes=vapply(gs,function(g)length(intersect(g,u)),integer(1)),mean=mu,sd=sig,universe_n=length(u)),'external_frozen_rank_training_parameters')
for(n in names(expr)){s<-rankscore(expr[[n]]);rownames(s)<-names(gs);colnames(s)<-colnames(expr[[n]]);z<-sweep(sweep(s,1,mu,'-'),1,sig,'/');pred<-call_scores(z);pred$Dataset<-n;write_result(pred,paste0(n,'_frozen_rank_transfer_sensitivity'))}
saveRDS(list(universe=u,gene_sets=gs,mean=mu,sd=sig,rule='within-sample gene-rank means; discovery-only Z scaling; ambiguity top-two delta<0.1; no external retuning'),file.path(OUT,'external_frozen_rank_model.rds'))
end_log()
