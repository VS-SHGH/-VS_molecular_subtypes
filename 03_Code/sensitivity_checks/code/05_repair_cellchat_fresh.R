source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'00_config.R'))
start_log('05_repair_cellchat_fresh')
suppressPackageStartupMessages({library(qs);library(SeuratObject);library(Matrix);library(CellChat);library(future)})
future::plan('sequential');options(future.globals.maxSize=20*1024^3);pbapply::pboptions(type='none')
sc<-qread(file.path(OUT,'singlecell_normalized_existing_cells.qs'))
md<-sc@meta.data;keep<-md$technique=='fresh';md<-md[keep,,drop=FALSE]
mat<-GetAssayData(sc,assay='RNA',layer='data')[,rownames(md),drop=FALSE]
stopifnot(length(unique(md$patient))==11L,max(mat@x)<20,!identical(GetAssayData(sc,layer='counts')[,rownames(md)],mat))
meta<-data.frame(labels=md$cellchat_label,samples=md$patient,row.names=rownames(md))
write_result(as.data.frame(table(meta$samples,meta$labels)),'cellchat_fresh_patient_group_counts')
cellchat<-createCellChat(object=mat,meta=meta,group.by='labels')
cellchat@DB<-CellChatDB.human;cellchat<-subsetData(cellchat)
cat('Rebuilding overexpressed-gene screening on fresh normalized expression',format(Sys.time()),'\n')
cellchat<-identifyOverExpressedGenes(cellchat,do.fast=FALSE)
cellchat<-identifyOverExpressedInteractions(cellchat)
cat('Computing all 100 permutations; no cached significance is read',format(Sys.time()),'\n')
cellchat<-computeCommunProb(cellchat,type='truncatedMean',trim=.1,nboot=100,seed.use=1L,population.size=FALSE)
cellchat<-filterCommunication(cellchat,min.cells=10);cellchat<-computeCommunProbPathway(cellchat);cellchat<-aggregateNet(cellchat)
cellchat<-netAnalysis_computeCentrality(cellchat,slot.name='netP')
saveRDS(cellchat,file.path(OUT,'cellchat_fresh_normalized_nboot100.rds'))
lr<-subsetCommunication(cellchat);write_result(lr,'cellchat_fresh_all_LR')
tumour<-c('nmSC','myeSC') # original source annotation is nmSC, not nSMC
focus<-lr[(lr$source%in%tumour&lr$target%in%c('C3high_Myeloid','C3low_Myeloid'))|(lr$target%in%tumour&lr$source%in%c('C3high_Myeloid','C3low_Myeloid')),]
write_result(focus,'cellchat_fresh_myeloid_tumour_LR')
net<-data.frame(group=rownames(cellchat@net$weight),outgoing=rowSums(cellchat@net$weight),incoming=colSums(cellchat@net$weight));net$total<-net$outgoing+net$incoming;net<-net[order(-net$total),];write_result(net,'cellchat_fresh_group_centrality_summary');print(net)
pw<-aggregate(focus$prob,list(pathway=focus$pathway_name,source=focus$source,target=focus$target),sum);names(pw)[4]<-'sum_probability';write_result(pw,'cellchat_fresh_focus_pathway_summary')

# Patient-level expression support is descriptive, not a second CellChat p-value
# nor proof of functional signalling. Complexes require every measured subunit.
resolve<-function(g){if(g%in%rownames(mat))return(g);if(g%in%rownames(cellchat@DB$complex)){z<-as.character(unlist(cellchat@DB$complex[g,,drop=FALSE]));return(z[!is.na(z)&nzchar(z)])};character()}
linear<-mat;linear@x<-expm1(linear@x)
pm<-interaction(md$patient,md$cellchat_label,drop=TRUE,sep='|');Z<-sparseMatrix(i=seq_along(pm),j=as.integer(pm),x=1,dims=c(length(pm),nlevels(pm)))
tot<-as.numeric(table(pm));avg<-as.matrix(linear%*%Z);avg<-sweep(avg,2,tot,'/');colnames(avg)<-levels(pm)
detect<-as.matrix((linear>0)%*%Z);detect<-sweep(detect,2,tot,'/');colnames(detect)<-levels(pm)
selected<-focus[focus$pathway_name%in%c('MIF','SPP1','APP'),];patient_rows<-list()
for(i in seq_len(nrow(selected))){r<-selected[i,];lig<-resolve(r$ligand);rec<-resolve(r$receptor);available<-length(lig)>0&&length(rec)>0&&all(c(lig,rec)%in%rownames(avg))
 for(p in unique(md$patient)){
  a<-paste(p,r$source,sep='|');b<-paste(p,r$target,sep='|');valid<-available&&all(c(a,b)%in%colnames(avg))
  if(valid)valid<-tot[match(a,colnames(avg))]>=10&&tot[match(b,colnames(avg))]>=10
  support<-if(valid)all(detect[lig,a]>=.1)&all(detect[rec,b]>=.1)else NA
  val<-if(valid)exp(mean(log(avg[lig,a]+1e-8))+mean(log(avg[rec,b]+1e-8)))else NA_real_
  patient_rows[[length(patient_rows)+1]]<-data.frame(patient=p,source=r$source,target=r$target,pathway=r$pathway_name,interaction=r$interaction_name,eligible=valid,expression_support_10pct=support,ligand_receptor_expression_product=val)
 }
}
if(length(patient_rows))write_result(do.call(rbind,patient_rows),'cellchat_fresh_patient_expression_support')
write_result(data.frame(n_patients=11,n_cells=ncol(mat),nboot=100,reused_cached_pvalues=FALSE,normalized=TRUE,patient_inference='pooled-cell permutations; not patient-level significance',tumour_label_correction='nmSC rather than nSMC'),'cellchat_fresh_run_contract')
end_log()
