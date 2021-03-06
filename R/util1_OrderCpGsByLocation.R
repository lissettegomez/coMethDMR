#' Order CpGs by genomic location
#'
#' @param CpGs_char vector of CpGs
#' @param arrayType Type of array, 450k or EPIC
#' @param output vector of CpGs or dataframe with CpGs, CHR, MAPINFO
#' @param genome Human genome of reference hg19 or hg38
#'
#' @return vector of CpGs ordered by location or dataframe with CpGs ordered by location (cpg),
#' chromosome (chr), position (pos)
#' @export
#'
#' @importFrom  sesameData sesameDataGet
#'
#' @examples
#'    CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
#'    OrderCpGsByLocation(
#'      CpGs_char,
#'      genome = "hg38",
#'      arrayType = "EPIC",
#'      output = "dataframe"
#'    )
#'
#'    CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
#'    OrderCpGsByLocation(
#'      CpGs_char,
#'      genome = "hg19",
#'      arrayType = "450k",
#'      output = "vector"
#'    )
#'
OrderCpGsByLocation <- function(
  CpGs_char,
  genome = c("hg19","hg38"),
  arrayType = c("450k","EPIC"),
  output = c("vector", "dataframe")
){

  arrayType <- match.arg(arrayType)
  genome <- match.arg(genome)

  # Available manifest files are
  # "EPIC.hg19.manifest"  "EPIC.hg38.manifest"
  # "HM27.hg19.manifest"  "HM27.hg38.manifest"
  # "HM450.hg19.manifest" "HM450.hg38.manifest"
  manifest <- paste(
    ifelse(arrayType == "450k","HM450","EPIC"),
    genome, "manifest",
    sep = "."
  )

  CpGlocations.gr <- sesameDataGet(manifest)

  CpGs.gr <- CpGlocations.gr[as.character(CpGs_char)] %>% sort

  ### Select and return output ###
  if (output == "dataframe") {
    CpGsOrdered_df <- as.data.frame(CpGs.gr)[,c("seqnames","start")]
    CpGsOrdered_df$cpg <- names(CpGs.gr)
    colnames(CpGsOrdered_df) <- c("chr","pos","cpg")
    CpGsOrdered_df
  } else {
    as.character(names(CpGs.gr))
  }
}
