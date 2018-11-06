
#' Extract contiguous comethylated CpG regions
#'
#' @param closeByGenomicRegionFile a file with the close by regions
#' @param betaMatrix matrix of beta values, with row names = CpG ids,
#'    column names = sample ids
#' @param arrayType Type of array, 450k or EPIC
#' @param fileType the output files can be saved as .gmt or .RDS,
#'    .gmt files can be open as flat text file, .RDS files are 50% the size of
#'    .gmt files, but cannot be open
#' @param dataDir link to the directory where the filtered files will be saved
#' @param returnAllCpGs indicates if outputting all the CpGs in the region
#'    when there is not a contiguous comethylated region or
#'    only the CpGs in the contiguous comethylated regions
#' @param ...
#'
#' @return a text file with a list of two items for each region:
#'    1. Contiguous_Regions - a data frame with CpG = CpG name,
#'    Chr = chromosome number,
#'    MAPINFO = genomic position,
#'    r_drop = correlation between the CpG with rest of the CpGs,
#'    keep = indicator for co-methylated CpG,
#'    keep_contiguous = cotiguous comethylated subregion number
#'    2. CpGs_subregions - lists of CpGs in each contiguous co-methylated
#'    subregion
#'    a file with the contiguous comethylated CpG regions
#' @export
#'
#' @examples
CoMethAllRegions <- function(closeByGenomicRegionFile, betaMatrix, arrayType = c("450k","EPIC"), returnAllCpGs = TRUE, fileType = c("gmt","RDS"), dataDir, ...){

  arrayType <- match.arg(arrayType)
  fileType <- match.arg(fileType)

  ### Read file of close by CpGs ###
  if (fileType == "RDS") {
    closeByGenomicRegion_ls <- readRDS(closeByGenomicRegionFile)
  } else {
    closeByGenomicRegion_ls <- read_gmt(closeByGenomicRegionFile)
  }

  ### Extract contiguous comethylated region(s) from each close by region ###
  coMethCpGsAllREgions_ls <- lapply(
    unname(closeByGenomicRegion_ls$pathways), CoMethSingleRegion, betaMatrix, arrayType
  )



  ### Return list of cotiguous comethylated CpGs by Regions ###
  out_coMethCpGsAll <- unlist(
    lapply(coMethCpGsAllREgions_ls, function(x) x[[2]]), recursive = F
  )

  ### Select and return output ###
  path_char <- paste(dataDir, "coMethAllRegions", sep = "/")
  fileName <- paste(path_char, fileType, sep = ".")
  if (fileType == "RDS") {
    saveRDS(out_coMethCpGsAll, file = fileName)
  } else {
    write_gmt(out_coMethCpGsAll, file = fileName)
  }


  out_coMethCpGsAll


}
