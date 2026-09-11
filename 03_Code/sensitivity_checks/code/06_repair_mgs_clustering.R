source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'00_config.R'))
start_log('06_repair_mgs_clustering')
suppressPackageStartupMessages({library(readxl);library(princurve);library(ConsensusClusterPlus);library(cluster);library(mclust)})
cl<-read.csv(file.path(OUT,'clinical_analysis_deidentified.csv'));tpm<-as.data.frame(read_excel(TPM_INPUT))
keep<-!is.na(tpm$gene_name)&nzchar(tpm$gene_name)&!duplicated(tpm$gene_name)
cols<-grep('^tpm\\.',names(tpm),value=TRUE);cols<-cols[sub('^tpm\\.','',cols)%in%cl$SampleID]
mat<-as.matrix(tpm[keep,cols,drop=FALSE]);storage.mode(mat)<-'numeric';rownames(mat)<-tpm$gene_name[keep];colnames(mat)<-sub('^tpm\\.','',cols)
ids<-colnames(mat);stopifnot(ncol(mat)==38,!anyNA(mat));cl<-cl[match(ids,cl$SampleID),];expr<-log2(mat+1)
ref<-load_env('Step2_Clustering.RData');step1<-load_env('Step1_Clean_Data_VS3.RData')
filter_genes<-rowSums(mat>1)>=ceiling(.2*ncol(mat));clean<-mat[filter_genes,,drop=FALSE]
cg<-intersect(rownames(clean),rownames(step1$clean_tpm));clean_err<-max(abs(clean[cg,ids]-step1$clean_tpm[cg,ids]))
write_result(data.frame(raw_gene_named_rows=nrow(mat),included_patients=38,TPM_threshold=1,operator='>',min_positive_samples=8,rederived_filtered_genes=nrow(clean),Step1_cached_genes=nrow(step1$clean_tpm),same_gene_set=setequal(rownames(clean),rownames(step1$clean_tpm)),max_abs_difference=clean_err),'bulk_step1_filter_reconstruction')
est<-read.csv(file.path(BASE,'outputs/estimate_official/estimate_official_scores.csv'));est<-est[match(ids,est$SampleID),];stopifnot(identical(ids,est$SampleID))
comp<-model.matrix(~Stromal_Score+Immune_Score,est)
anchor_genes<-c('CD2','CD3D','CD3E','CD8A','GZMA','PRF1','CTLA4','HAVCR2','LAG3','PDCD1','CD68','CD163','CD74','HLA-DRA','CCL5','CXCL9','CXCL10')
anchor<-colMeans(expr[intersect(anchor_genes,rownames(expr)),,drop=FALSE])
resid_expr<-t(qr.resid(qr(comp),t(expr)))
get_top<-function(e){v<-apply(e,1,mad);names(sort(v[v>0&is.finite(v)],decreasing=TRUE))[seq_len(min(3000,sum(v>0&is.finite(v))))]}
variants<-list(original=expr,anchor_genes_excluded=expr[!rownames(expr)%in%anchor_genes,,drop=FALSE],composition_residualized=resid_expr)
tops<-lapply(variants,get_top)
fit_mgs<-function(e,a){sdg<-apply(e,1,sd);e<-e[is.finite(sdg)&sdg>1e-10,,drop=FALSE];p<-prcomp(t(e),center=TRUE,scale.=TRUE);f<-principal_curve(as.matrix(p$x[,1:5]),smoother='smooth_spline',stretch=0);r<-suppressWarnings(cor(f$lambda,a,method='spearman'));stopifnot(is.finite(r));x<-(f$lambda-min(f$lambda))/diff(range(f$lambda));if(r<0)x<-1-x;list(score=setNames(x,colnames(e)),pca=p,curve=f,anchor_rho=r)}
fits<-lapply(names(variants),function(nm)fit_mgs(variants[[nm]][tops[[nm]],],anchor));names(fits)<-names(variants)
score_df<-data.frame(SampleID=ids,Subtype=cl$Subtype,Stromal_Score=est$Stromal_Score,Immune_Score=est$Immune_Score)
for(nm in names(fits))score_df[[nm]]<-fits[[nm]]$score[ids]
write_result(score_df,'MGS_repaired_sensitivity_scores')
orig<-score_df$original
varmodel<-lm(orig~Stromal_Score+Immune_Score,est);r2<-summary(varmodel)$r.squared
r2flip<-summary(lm(I(1-orig)~Stromal_Score+Immune_Score,est))$r.squared
part1<-summary(lm(orig~Stromal_Score,est))$r.squared
write_result(data.frame(component=c('Stromal_first','Immune_increment_second','Residual','Joint_R2','Joint_R2_after_orientation_flip'),value=c(part1,r2-part1,1-r2,r2,r2flip),interpretation='internal same-transcriptome association; sequential allocation order-dependent; not causation'),'MGS_variance_corrected_interpretation')
base<-readRDS(file.path(BASE,'outputs/mgs/mgs_unbiased_results.rds'))$mgs_unbiased
baseids<-sub('^tpm\\.','',base$Sample);stopifnot(all(ids%in%baseids));basevalue<-base$MGS_unbiased[match(ids,baseids)]
write_result(data.frame(metric=c('max_abs_difference_from_v5_MGS','joint_R2','joint_R2_flip_difference'),value=c(max(abs(orig-basevalue)),r2,abs(r2-r2flip))),'MGS_baseline_reproduction_check')
compres<-do.call(rbind,lapply(names(fits),function(nm){v<-score_df[[nm]];data.frame(variant=nm,spearman_with_original=cor(v,orig,method='spearman'),kendall_with_original=cor(v,orig,method='kendall'),C1_mean=mean(v[cl$Subtype=='C1']),C2_mean=mean(v[cl$Subtype=='C2']),C3_mean=mean(v[cl$Subtype=='C3']),subtype_KW_p=kruskal.test(v~cl$Subtype)$p.value)}))
write_result(compres,'MGS_variant_comparison');print(compres)

# Conditional stability given each full-data selected gene set. In the residual
# variant the composition regression is refitted within each bootstrap draw.
boot<-list();B<-500L
for(nm in names(fits)){
 set.seed(20260905);rows<-vector('list',B)
 for(b in seq_len(B)){
  ii<-sample.int(38,replace=TRUE);pp<-ids[ii]
  e<-if(nm=='composition_residualized')t(qr.resid(qr(comp[ii,,drop=FALSE]),t(expr[tops[[nm]],ii,drop=FALSE])))else variants[[nm]][tops[[nm]],ii,drop=FALSE]
  rows[[b]]<-tryCatch({f<-fit_mgs(e,anchor[ii]);v<-tapply(f$score,pp,median);tau<-cor(v,fits[[nm]]$score[names(v)],method='kendall');data.frame(variant=nm,replicate=b,tau=tau,n_unique=length(v),status='ok')},error=function(e)data.frame(variant=nm,replicate=b,tau=NA_real_,n_unique=length(unique(pp)),status=conditionMessage(e)))
 }
 boot[[nm]]<-do.call(rbind,rows);cat('bootstrap completed',nm,'\n')
}
boot<-do.call(rbind,boot);write_result(boot,'MGS_variant_bootstrap_replicates')
bs<-do.call(rbind,lapply(split(boot,boot$variant),function(z)data.frame(variant=z$variant[1],successful=sum(is.finite(z$tau)),attempted=nrow(z),mean_tau=mean(z$tau,na.rm=TRUE),q025=quantile(z$tau,.025,na.rm=TRUE),q975=quantile(z$tau,.975,na.rm=TRUE),definition='conditional subject bootstrap; selected genes fixed; interval of tau distribution, not CI for a biological effect')))
write_result(bs,'MGS_variant_bootstrap_summary');saveRDS(list(fits=fits,selected_genes=tops,variance_model=varmodel),file.path(OUT,'MGS_repaired_fits.rds'));print(bs)

# Refit consensus clustering for k=2--6 and report, rather than force, stability.
e<-ref$clust_input;set.seed(123456)
cc<-ConsensusClusterPlus(e,maxK=6,reps=1000,pItem=.8,pFeature=1,clusterAlg='hc',distance='pearson',seed=123456,title=file.path(OUT,'ccp_no_plots'),plot=NULL,verbose=FALSE)
distance<-as.dist(1-cor(e));ks<-lapply(2:6,function(k){cm<-cc[[k]]$consensusMatrix;lab<-cc[[k]]$consensusClass;v<-cm[lower.tri(cm)];sil<-silhouette(lab,distance);data.frame(k=k,PAC_01_09=mean(v>.1&v<.9),silhouette_mean=mean(sil[,'sil_width']),cluster_sizes=paste(as.integer(table(lab)),collapse='/'),C2_membership_note='labels and clinical endpoints not used to optimize k')})
write_result(do.call(rbind,ks),'clustering_k2_to_k6_recomputed');saveRDS(cc,file.path(OUT,'clustering_k2_to_k6_recomputed.rds'))
labs<-data.frame(SampleID=names(cc[[3]]$consensusClass),original_Subtype=ref$meta_step2$Subtype[match(names(cc[[3]]$consensusClass),ref$meta_step2$SampleID)])
for(k in 2:6)labs[[paste0('k',k)]]<-cc[[k]]$consensusClass[labs$SampleID]
write_result(labs,'clustering_k2_to_k6_patient_labels')
write_result(data.frame(metric='k3_adjusted_Rand_against_original',value=adjustedRandIndex(labs$k3,labs$original_Subtype)),'clustering_label_reproduction')
end_log()
