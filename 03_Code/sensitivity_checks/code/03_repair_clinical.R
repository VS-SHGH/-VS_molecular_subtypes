source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'00_config.R'))
start_log('03_repair_clinical')
suppressPackageStartupMessages({library(readxl);library(MASS);library(logistf);library(sandwich);library(lmtest);library(pROC)})
raw <- as.data.frame(read_excel(file.path(DATA,'cl_bulk_with_Subtype.xlsx')))
units <- raw[!raw$Subtype%in%c('C1','C2','C3'),,drop=FALSE]
cl <- raw[raw$Subtype%in%c('C1','C2','C3'),,drop=FALSE]
su <- as.data.frame(read_excel(file.path(DATA,'VS 手术分析.xlsx')))
stopifnot(nrow(cl)==38,!anyDuplicated(cl$SampleID),!anyDuplicated(su$SampleID),setequal(cl$SampleID,su$SampleID))
su <- su[match(cl$SampleID,su$SampleID),,drop=FALSE]
stopifnot(identical(su$SampleID,cl$SampleID))
parse_size <- function(s){v<-as.numeric(strsplit(gsub('×','*',s,fixed=TRUE),'\\*')[[1]]);stopifnot(length(v)==3,!anyNA(v));max(v)}
size_mm <- vapply(cl$TumorSize_raw,parse_size,numeric(1));size_mm[cl$SampleID=='P32']<-22
df <- data.frame(SampleID=cl$SampleID,Subtype=factor(cl$Subtype,levels=c('C1','C2','C3')),Age=as.numeric(cl$Age),size_raw=cl$TumorSize_raw,size_mm=size_mm,size_cm=size_mm/10,size_provenance=ifelse(cl$SampleID=='P32','author_confirmed_20260906','parsed_original_mm'),texture=ordered(unname(c('软'='soft','中等'='medium','硬'='hard')[su$肿瘤质地]),levels=c('soft','medium','hard')),adhesion=unname(c('有'=1,'无'=0)[su$脑干面粘连]),NK=as.numeric(cl$自然杀伤细胞),size_large=as.integer(cl$SizeGrade=='large'))
stopifnot(!anyNA(df),sum(df$adhesion)==10)
write_result(df,'clinical_analysis_deidentified');print(table(df$Subtype,df$texture));print(table(df$Subtype,df$adhesion))
unit_out<-data.frame(variable=names(cl)[15:91],unit=as.character(units[1,15:91]))
write_result(unit_out,'clinical_laboratory_units');print(unit_out[unit_out$variable%in%c('自然杀伤细胞','总T淋巴细胞'),])

results<-list();joint<-list();modelcache<-list()
add <- function(m,label,outcome,kind){
 if(kind=='firth'){
  b<-coef(m);se<-sqrt(diag(m$var));lo<-m$ci.lower;hi<-m$ci.upper;p<-m$prob;ci_method<-'penalized_profile_likelihood';v<-m$var
 }else{
  b<-coef(m);se<-sqrt(diag(vcov(m)))[names(b)];v<-vcov(m)[names(b),names(b),drop=FALSE]
  if(kind=='ordinal'){
   ci<-tryCatch(suppressMessages(confint(m)),error=function(e)NULL)
   if(is.null(ci)||any(!is.finite(ci))){lo<-b-1.96*se;hi<-b+1.96*se;ci_method<-'Wald_profile_failed'}else{lo<-ci[names(b),1];hi<-ci[names(b),2];ci_method<-'profile_likelihood'}
   p<-2*pnorm(-abs(b/se))
  }else{ci<-confint(m);lo<-ci[names(b),1];hi<-ci[names(b),2];p<-summary(m)$coefficients[names(b),4];ci_method<-'t_interval'}
 }
 k<-grep('^Subtype',names(b));w<-as.numeric(t(b[k])%*%solve(v[k,k,drop=FALSE],b[k]));jp<-pchisq(w,length(k),lower.tail=FALSE)
 joint[[label]]<<-data.frame(model=label,outcome=outcome,method=kind,n=nrow(df),joint_subtype_Wald=w,df=length(k),p=jp)
 results[[label]]<<-data.frame(model=label,outcome=outcome,method=kind,n=nrow(df),term=names(b),beta=unname(b),SE=unname(se),estimate=if(kind=='linear')unname(b)else exp(unname(b)),CI_low=if(kind=='linear')unname(lo)else exp(unname(lo)),CI_high=if(kind=='linear')unname(hi)else exp(unname(hi)),p=unname(p),CI_method=ci_method)
 modelcache[[label]]<<-m
}
for(label in c('unadjusted','size_adjusted','size_age_sensitivity')){
 rhs<-switch(label,unadjusted='Subtype',size_adjusted='Subtype + size_cm',size_age_sensitivity='Subtype + size_cm + Age')
 fm<-logistf(as.formula(paste('adhesion ~',rhs)),data=df,pl=TRUE,control=logistf.control(maxit=1000),plcontrol=logistpl.control(maxit=1000));add(fm,paste0('adhesion_',label),'brainstem_adhesion','firth')
 om<-polr(as.formula(paste('texture ~',rhs)),data=df,Hess=TRUE,method='logistic');add(om,paste0('texture_',label),'texture_soft_medium_hard','ordinal')
}
for(label in c('age_adjusted','age_size_binary_legacy','age_size_continuous')){
 rhs<-switch(label,age_adjusted='Subtype + Age',age_size_binary_legacy='Subtype + Age + size_large',age_size_continuous='Subtype + Age + size_cm')
 lm0<-lm(as.formula(paste('NK ~',rhs)),df);add(lm0,paste0('NK_',label),'NK_recorded_laboratory_measure','linear')
}
rr<-do.call(rbind,results);jj<-do.call(rbind,joint)
jj$q_primary_two_surgical_endpoints<-NA_real_;i<-jj$model%in%c('adhesion_size_adjusted','texture_size_adjusted');jj$q_primary_two_surgical_endpoints[i]<-p.adjust(jj$p[i],'BH')
write_result(rr,'clinical_model_coefficients');write_result(jj,'clinical_omnibus_tests');saveRDS(modelcache,file.path(OUT,'clinical_model_fits.rds'))
print(rr[grepl('size_adjusted|age_size_continuous',rr$model)&grepl('Subtype',rr$term),]);print(jj)
# Proportional-odds diagnostic: compare subtype log ORs at both thresholds.
th<-list();for(t in c(1,2)){z<-df;z$harder<-as.integer(as.integer(z$texture)>t);m<-logistf(harder~Subtype+size_cm,data=z,pl=TRUE);th[[t]]<-data.frame(threshold=t,term=names(coef(m)),OR=exp(coef(m)),low=exp(m$ci.lower),high=exp(m$ci.upper),p=m$prob)}
write_result(do.call(rbind,th),'texture_threshold_sensitivity')
nkfit<-modelcache$NK_age_size_continuous;hc<-coeftest(nkfit,vcov.=vcovHC(nkfit,type='HC3'))
write_result(data.frame(term=rownames(hc),beta=hc[,1],SE_HC3=hc[,2],t=hc[,3],p=hc[,4]),'NK_continuous_size_HC3_sensitivity')
write_result(data.frame(shapiro_p=shapiro.test(residuals(nkfit))$p.value,BP_p=bptest(nkfit)$p.value,max_cooks=max(cooks.distance(nkfit)),model_R2=summary(nkfit)$r.squared),'NK_linear_diagnostics')

# Internal marker-selection sensitivity. The historical, complete sequence of
# exploratory decisions is not reconstructable; this is not a complete correction
# of every researcher choice and is NOT external validation of NK or subtypes.
laboratory<-names(cl)[15:91]
num <- function(x){s<-trimws(as.character(x));s[grepl('^[<>≤≥]',s)]<-NA_character_;suppressWarnings(as.numeric(s))}
X<-as.data.frame(lapply(cl[,laboratory,drop=FALSE],num),check.names=FALSE)
X[['NK_to_total_T_ratio']]<-X[['自然杀伤细胞']]/X[['总T淋巴细胞']]
X<-as.matrix(X);X[!is.finite(X)]<-NA_real_;y<-as.integer(df$Subtype=='C3')
auc0<-function(y,s){if(length(unique(y))!=2)return(NA_real_);r<-rank(s,ties.method='average');n1<-sum(y==1);n0<-sum(y==0);(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)}
train_select<-function(idx,only_nk=FALSE){
 xx<-X[idx,,drop=FALSE];eligible<-which(colMeans(is.finite(xx))>=.8&apply(xx,2,function(z)length(unique(z[is.finite(z)]))>=5))
 if(only_nk)eligible<-match('自然杀伤细胞',colnames(X))
 med<-apply(xx[,eligible,drop=FALSE],2,median,na.rm=TRUE)
 aa<-vapply(seq_along(eligible),function(i){z<-xx[,eligible[i]];z[is.na(z)]<-med[i];auc0(y[idx],z)},numeric(1));strength<-pmax(aa,1-aa)
 w<-which.max(strength);j<-eligible[w];z<-X[idx,j];z[is.na(z)]<-med[w];sgn<-ifelse(aa[w]>=.5,1,-1)
 list(j=j,median=med[w],direction=sgn,train_auc=strength[w],fit=logistf(y~x,data=data.frame(y=y[idx],x=z),pl=FALSE))
}
predict_train<-function(m,ids){z<-X[ids,m$j];z[is.na(z)]<-m$median;as.numeric(predict(m$fit,newdata=data.frame(x=z),type='response'))}
nkroc<-pROC::roc(y,df$NK,levels=c(0,1),direction='<',quiet=TRUE);ci<-as.numeric(ci.auc(nkroc));full<-train_select(seq_len(nrow(df)))
write_result(data.frame(metric=c('NK_apparent_AUC','NK_apparent_DeLong_low','NK_apparent_DeLong_high','selected_marker_apparent_AUC'),value=c(as.numeric(auc(nkroc)),ci[1],ci[3],full$train_auc),note=c('same 38 patients, no selection adjustment','not selection-adjusted','not selection-adjusted',colnames(X)[full$j])),'NK_apparent_AUC')
set.seed(20260906);reps<-100L;cv<-list();choices<-list()
for(b in seq_len(reps)){
 fold<-integer(length(y));for(g in 0:1){ids<-which(y==g);fold[ids]<-sample(rep(1:5,length.out=length(ids)))}
 pp<-pn<-rep(NA_real_,length(y))
 for(k in 1:5){tr<-which(fold!=k);te<-which(fold==k);m<-train_select(tr);n<-train_select(tr,TRUE);pp[te]<-predict_train(m,te);pn[te]<-predict_train(n,te);choices[[length(choices)+1]]<-data.frame(repetition=b,fold=k,marker=colnames(X)[m$j])}
 cv[[b]]<-data.frame(repetition=b,selected_marker_OOF_AUC=auc0(y,pp),fixed_NK_OOF_AUC=auc0(y,pn))
 if(b%%20==0)cat('CV',b,'/',reps,'\n')
}
cv<-do.call(rbind,cv);write_result(cv,'NK_selection_aware_repeated_CV');write_result(do.call(rbind,choices),'NK_training_fold_selected_markers')
write_result(data.frame(method=c('selection_within_training_fold','fixed_NK_only'),mean_AUC=c(mean(cv$selected_marker_OOF_AUC),mean(cv$fixed_NK_OOF_AUC)),split_q025=c(quantile(cv$selected_marker_OOF_AUC,.025),quantile(cv$fixed_NK_OOF_AUC,.025)),split_q975=c(quantile(cv$selected_marker_OOF_AUC,.975),quantile(cv$fixed_NK_OOF_AUC,.975)),interval_type='distribution across repeated splits, NOT a population 95% confidence interval',patients=38,folds=5,repetitions=100),'NK_CV_summary')
end_log()
