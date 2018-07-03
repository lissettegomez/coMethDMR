# Input:
# regionsBeta - a matrix of beta values for one region, with rownames = cpg ids, column names = sample ids
# pheno_table - table with phenotype and covariates, the phenotype variable should be called "pheno" in this table. All the covariats will be included in the model.


mxmodel<-function(regionsBeta,pheno_table)  {
  regionsBeta$ProbeID<-row.names(regionsBeta)
  regionsBeta_t<-reshape(regionsBeta, varying = colnames(regionsBeta[-dim(regionsBeta)[2]]), v.names="beta", direction="long", time=colnames(regionsBeta[-dim(regionsBeta)[2]]),timevar="Sample")
  temp<-merge(regionsBeta_t,pheno_table,by="Sample")
  temp$Mvalue<-log2(temp$beta/(1-temp$beta))
  ps=999
  tryCatch({
    f<-lmer(Mvalue ~pheno + (pheno|ProbeID) +(1|Sample) + . , temp) #### fix the covariates
    ps<-coef(summary(f))[2,5]
    slope<-coef(summary(f))[2,1]
  }, error=function(e){})
  correlation<-cor(t(regionsBeta[-dim(regionsBeta)[2]]),method="spearman", use="pairwise.complete.obs")
  median.corr<-median(correlation[lower.tri(correlation)])
  l<-list(slope,ps,median.corr)
  names(l)<-c("slope","Pval.mixed.model","median.corr")
  return(l)
}
