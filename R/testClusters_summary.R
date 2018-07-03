# input
# list_regions - list of regions of cpgs
# betaMatrix - a matrix of beta values, with rownames = cpg ids, column names = sample ids
# pheno_table - table with phenotype and covariates, the phenotype variable should be called "pheno" in this table. All the covariats will be included in the model.

# output = data frame with region, chr, start, end, estimate, pvalue, gene_name

testClusters_summary <-function(list_regions, betaMatrix, pheno_table){

  regionsBeta<-lapply(list_regions, function (item) betaMatrix[which(rownames(betaMatrix)%in%item),])

  library(lmerTest)

  mxmodel_pvalue<-lapply(regionsBeta, mxmodel, pheno_table, pheno, covariates)


  split.name<-lapply(names(mxmodel_pvalue), strsplit, split="[.]")
  split.name.unlist<-lapply(split.name, unlist)
  region.id<-lapply(split.name.unlist, function(x) x[[1]])
  #close.subregion.id<-lapply(pfc.closeCpGs.regions.final,region.name,location.file)
  contiguous.subregion.id<-lapply(split.name.unlist, function(x) x[[2]])
  results<-data.frame("region.id"=unlist(region.id),
                      "contiguous.subregion.id"=unlist(contiguous.subregion.id),
                      "slope"=unlist(unname(lapply(mxmodel_pvalue, function(item) item$`slope`))),
                      "Pval.mixed.model"=unlist(unname(lapply(mxmodel_pvalue, function(item) item$`Pval.mixed.model`))),
                      "median.corr"=unlist(unname(lapply(mxmodel_pvalue, function(item) item$`median.corr`))))

  return(results)





}
