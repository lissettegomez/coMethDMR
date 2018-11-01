
#' Extract clusters of close by CpGs from a genomic region
#'
#' @param CpGs_char a list of CpG ids
#' @param arrayType Type of array, 450k or EPIC
#' @param maxGap an integer, genomic locations within maxGap from each other
#'    are placed into the same cluster
#' @param minCpGs an integer, minimum number of CpGs for resulting regions
#'
#' @return a list, each item is a list of CpG ids in the contiguous region
#'    (i.e. same cluster)
#'
#' @importFrom bumphunter clusterMaker
#'
#' @export
#'
#' @examples
#'    CpGs_char <- c("cg02505293", "cg03618257", "cg04421269", "cg17885402",
#'        "cg19890033", "cg20566587", "cg27505880")
#'    CloseBySingleRegion(CpGs_char, arrayType="450k", maxGap=100, minCpGs=3)
#'
CloseBySingleRegion <- function(CpGs_char, arrayType=c("450k","EPIC"), maxGap = 200, minCpGs = 3){

  CpGsOrdered_df <- OrderCpGsByLocation(CpGs_char, arrayType, output = "dataframe")

  ### Find close by clusters ###
  chr <- CpGsOrdered_df$chr
  pos <- CpGsOrdered_df$pos
  CpGsOrdered_df$cluster <- clusterMaker(chr, pos, maxGap = maxGap)

  ### Create list of vectors of CpGs in each cluster ###
  CpGsRegion_ls <- split(CpGsOrdered_df$cpg, CpGsOrdered_df$cluster)

  ### Filter for clusters with number of CpGs >= minCpGs ###
  CpGsRegionMinCpGs_ls <- CpGsRegion_ls[lengths(CpGsRegion_ls) >= minCpGs]

  if (length(CpGsRegionMinCpGs_ls) > 0 ){
    CpGsRegionMinCpGs_ls
  } else {
    NULL
  }

}
