
##########################################################################################################

#### Input
# CpGs = a list of CpG ids
# location.file = a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg
# maxGap = an integer, genomic locations within maxGap from each other are placed into the same cluster
# min.cpgs = an integer, minimum nubmer of cpgs for resulting regions

#### Output
# a list, each item is a list of CpG ids in the contiguous region (i.e. same cluster)

##############################################################################################################


#' Extract clusters of close CpGs
#'
#' @param CpGs a list of CpG ids
#' @param location.file a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg
#' @param maxGap an integer, genomic locations within maxGap from each other are placed into the same cluster
#' @param min.cpgs an integer, minimum nubmer of cpgs for resulting regions
#'
#' @return a list, each item is a list of CpG ids in the contiguous region (i.e. same cluster)
#'
#' @importFrom bumphunter clusterMaker
#'
#' @export
#'
#' @examples
#'
closeCpGs.regions <- function(CpGs, location.file, maxGap, min.cpgs){

  #library(bumphunter)

  tryCatch({

      CpGs<-unlist(CpGs,recursive = F)

      CpGs.location <-location.file[which(location.file$cpg%in%CpGs),]

      CpGs.location <- CpGs.location[order(as.numeric(as.character(CpGs.location$MAPINFO))),]

      chr<-as.numeric(as.character(CpGs.location$CHR))

      pos<-as.numeric(as.character(CpGs.location$MAPINFO))

      CpGs.location$cluster <- clusterMaker(chr, pos, maxGap = maxGap)

      cluster.num <- unique(CpGs.location$cluster)

      cpgs.region.list <- list()

      for (i in cluster.num [1:length(cluster.num)]){

        region <- subset(CpGs.location, cluster== i);

        cpgs.region.list[[i]] <- factor(region$cpg)}  ## need to factor

  }, error=function(e){})

  cpgs.region.list3 <- cpgs.region.list[lapply(cpgs.region.list, length) >= min.cpgs]

  if (length(cpgs.region.list3) >0 ){
    return(cpgs.region.list3)
  }

  return ()

}
