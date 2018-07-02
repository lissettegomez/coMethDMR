
# Input
# list_regions - list of predefined regions, e.g. cpg island
# location.file - a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg
# maxGap - an integer, genomic locations within maxGap from each other are placed into the same cluster
# min.cpgs - an integer, minimum nubmer of cpgs for resulting regions
# betaMatrix - a matrix of beta values, with rownames = cpg ids, column names = sample ids
# min.region.size - min number of cpgs in the sub-regions
# threshold.r - min correlation between a cpg in the region with the rest of the CpGs

# Output = list of contiguous co-methylated regions



CoMethClusters <- function (list_regions, location.file, maxGap, min.cpgs, betaMatrix, min.region.size, threshold.r){
  colnames(location.file)[1] <- "cpg"
  closeCpGs.regions<-lapply(list_regions,closeCpGs.regions,location.file,maxGap=maxGap,min.cpgs=min.cpgs)
  closeCpGs.regions.final<-unlist(closeCpGs.regions,recursive=F)
  names(closeCpGs.regions.final)<-lapply(closeCpGs.regions.final,region.name,location.file)

  #rownames(betaMatrix)<-betaMatrix[,1]
  #betaMatrix<-betaMatrix[,-1]
  matrix.cgi<-lapply(closeCpGs.regions.final, function (item) betaMatrix[which(rownames(betaMatrix)%in%item),])

  contiguous.regions.out<-lapply(matrix.cgi,contiguous.regions, min.region.size=min.region.size, location.file=location.file, threshold.r=threshold.r)
  contiguous.filtered<-lapply(contiguous.regions.out, function(item) item$`n.sub.regions` > 0)
  contiguous.filtered.vec<-unlist(contiguous.filtered, recursive = FALSE)
  contiguous.filtered.regions<-contiguous.regions.out[contiguous.filtered.vec]
  contiguous.filtered.cpgs<-lapply(contiguous.filtered.regions, function(x) x[[3]])
  contiguous.filtered.cpgs.vec<-unlist(contiguous.filtered.cpgs,recursive = FALSE)

  split.name<-lapply(names(contiguous.filtered.cpgs.vec), strsplit, split="[.]")
  split.name.unlist<-lapply(split.name, unlist)
  region.id<-lapply(split.name.unlist, function(x) x[[1]])
  names(contiguous.filtered.cpgs.vec)<-region.id

  return(contiguous.filtered.cpgs.vec)
}

#test<-CoMethClusters(list_regions=cgi[1:1000], location.file=location.file, maxGap=200, min.cpgs=3, betaMatrix=betaMatrix, min.region.size=3, threshold.r=0.5)
