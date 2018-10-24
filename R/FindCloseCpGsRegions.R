
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
#'    CpGsChr22_char<-c("cg02953382", "cg12419862", "cg24565820", "cg04234412",
#'       "cg04824771", "cg09033563", "cg10150615", "cg18538332", "cg20007245",
#'       "cg23131131", "cg25703541")
#'
closeCpGs.regions <- function(CpGs_char, arrayType=c("450k","EPIC"), maxGap=200, min.cpgs=3){

  arrayType <- match.arg(arrayType)

  switch(arrayType,
         "450k" = {
           CpGlocations_df = IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
         },
         "EPIC" = {
           CpGlocations_df = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations
         }
  )

    tryCatch({

      ### Subset the location Data Frame  ###
      CpGlocations_df <- as.data.frame(CpGlocations_df)
      CpGlocations_df$cpg <- row.names(CpGlocations_df)
      row.names(CpGlocations_df) <- NULL
      CpGs_df <- CpGlocations_df[which(CpGlocations_df$cpg%in%CpGs_char),]

      ###  Re-order this Subset  ###
      order_idx <- order(CpGs_df$chr, CpGs_df$pos)
      CpGsOrdered_df <- CpGs_df[order_idx, ]

      chr<-CpGsOrdered_df$chr

      pos<-CpGsOrdered_df$pos

      CpGsOrdered_df$cluster <- clusterMaker(chr, pos, maxGap)

      clusterNum <- unique(CpGsOrdered_df$cluster)

      CpGsRegion_ls <- list()

      for (i in clusterNum [1:length(clusterNum)]){

        region <- subset(CpGsOrdered_df, cluster== i);

        CpGsRegion_ls[[i]] <- factor(region$cpg)}  ## need to factor

  }, error=function(e){})

  cpgs.region.list3 <- CpGsRegion_ls[lapply(CpGsRegion_ls, length) >= min.cpgs]

  if (length(cpgs.region.list3) >0 ){
    return(cpgs.region.list3)
  }

  return ()

}
