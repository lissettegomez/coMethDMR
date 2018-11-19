
#' Extract contiguous comethylated CpG regions
#'
#' @param betaMatrix matrix of beta values, with row names = CpG ids,
#'    column names = sample ids
#' @param file name of input file with close by regions
#' @param fileType extension of input file: "gmt" or "RDS"
#' @param arrayType Type of array, 450k or EPIC
#' @param returnAllCpGs indicates if outputting all the CpGs in the region
#'    when there is not a contiguous comethylated region or
#'    only the CpGs in the contiguous comethylated regions
#' @param ...
#'
#' @return a list of two items:
#'    1. Contiguous_Regions - a data frame with
#'    CpG = CpG name,
#'    Chr = chromosome number,
#'    MAPINFO = genomic position,
#'    r_drop = correlation between the CpG with rest of the CpGs,
#'    keep = indicator for co-methylated CpG,
#'    keep_contiguous = cotiguous comethylated subregion number
#'    for all the regions
#'    2. CpGs_subregions - lists of CpGs in each contiguous co-methylated
#'    subregion from all the regions
#' @export
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    CoMethAllRegions(
#'       betaMatrix = betaMatrixChr22_df,
#'       file = system.file(
#'          "extdata",
#'          "CpGislandsChr22_ex.RDS",
#'          package = 'coMethDMR',
#'          mustWork = TRUE),
#'       fileType = "RDS",
#'       arrayType = "450k",
#'       returnAllCpGs = FALSE
#'       )
#'
CoMethAllRegions <- function(betaMatrix,
                             file, fileType = c("gmt","RDS"),
                             arrayType = c("450k","EPIC"),
                             returnAllCpGs = FALSE,
                             ...){

  arrayType <- match.arg(arrayType)
  fileType <- match.arg(fileType)

  ### Read file of close by CpGs ###
  switch(
    fileType,
    "RDS" = {
      closeByGenomicRegion_ls <- readRDS(file)
    },
    "gmt" = {
      closeByGenomicRegion_ls <- read_gmt(file)
    }
  )


  ### Extract contiguous comethylated region(s) from each close by region ###
  # A pathway is a set of CpGs in a region
  coMethCpGsAllREgions_ls <- lapply(
    unname(closeByGenomicRegion_ls$pathways),
    FUN = CoMethSingleRegion,
    betaMatrix, arrayType, returnAllCpGs
  )


  ### Return list of cotiguous comethylated CpGs by Regions ###
  out_ContigRegions <- (lapply(coMethCpGsAllREgions_ls, `[[`, 1))
  out_coMethCpGsAll <- unlist(
    lapply(coMethCpGsAllREgions_ls, `[[`, 2),
    recursive = FALSE
  )


  ### return output ###
  list(
    CpGsSubregions    = out_coMethCpGsAll,
    contiguousRegions = out_ContigRegions
  )

}
