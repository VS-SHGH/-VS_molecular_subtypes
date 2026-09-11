# Prepare the submission table from computed, unrounded cross-method results.
# No correlation values are inferred from the earlier rounded submission table.
a <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(sub('^--file=','',a[1])),'..','..'))
out <- file.path(root,'outputs','cibersort_mcp')
d <- read.csv(file.path(out,'cibersort_mcp_correlation.csv'),check.names=FALSE)
d$Cell_Population[d$Cell_Population=='B cells (naive+memory)'] <- 'B cells naive'
d$Cell_Population[d$Cell_Population=='NK cells (resting+activated)'] <- 'NK cells resting'
fields <- c('Cell_Population','MCP_CellType','CIBERSORT_CellType','N','Spearman_Rho','P_value','Mean_MCP','Mean_CIBERSORT')
stopifnot(nrow(d)==13L, all(d$N==38L), all(fields %in% names(d)))
write.csv(d[,fields],file.path(root,'02_Submission_Materials','03_Supplementary_Tables','Supplementary_Table3_CrossValidation.csv'),row.names=FALSE)

# Recover the two P-value conventions for the displayed CD8 correlation.
# For n=38, R's AS89 route is an approximation, not an exact permutation test.
i <- which(d$Cell_Population=='CD8 T cells')
n <- d$N[i]; rho <- d$Spearman_Rho[i]; S <- (n^3-n)*(1-rho)/6
stopifnot(abs(S-round(S))<1e-6)
p_as89 <- min(2*.Call(stats:::C_pRho, round(S)+2,as.integer(n),TRUE),1)
write.csv(data.frame(N=n,Spearman_rho=rho,S=S,P_AS89=p_as89,P_asymptotic=d$P_value[i]),
          file.path(out,'cd8_pvalue_conventions.csv'),row.names=FALSE)
print(data.frame(N=n,Spearman_rho=rho,P_AS89=p_as89,P_asymptotic=d$P_value[i]))
