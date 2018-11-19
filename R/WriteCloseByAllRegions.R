

#' Extract clusters of close by CpGs from all genomic regions
#'
#' @param genomicRegionType Type of region
#' @param arrayType Type of array, 450k or EPIC
#' @param maxGap an integer, genomic locations within maxGap from each other
#'    are placed into the same cluster
#' @param minCpGs an integer, minimum number of CpGs for resulting regions
#' @param fileType the output files can be saved as .gmt or .RDS,
#'    .gmt files can be open as flat text file, .RDS files are 50% the size of
#'    .gmt files, but they can only be read in R enviroment.
#' @param file file where the close by regions will be saved
#' @param ...
#'
#' @return a file with the close by regions
#'
#' @details For \code{maxGap} = 200 and \code{minCpGs} = 3, we have calculated
#'    the close-by regions output files already. They are in /inst/extdata/
#'    File extension should not be supplied. We add it automatically through
#'    the fileType argument.
#'    Creating and writing this file will take around 25 minutes for
#'    genomicRegionType = "ISLAND"
#'
#'
#' @importFrom pathwayPCA read_gmt write_gmt CreatePathwayCollection
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   CloseByAllRegions(
#'      genomicRegionType = "ISLAND", arrayType = "450k", maxGap = 50,
#'      minCpGs = 3, fileType = "gmt", file = "closeByRegions"
#'   )
#' }
#'
WriteCloseByAllRegions <- function(
  file,
  genomicRegionType = c("ISLAND", "TSS1500", "TSS200", "UTR5", "EXON1",
                        "GENEBODY", "UTR3", "NSHORE", "SSHORE", "NSHELF", "SSHELF"),
  arrayType = c("450k","EPIC"),
  maxGap = 200,
  minCpGs = 3,
  fileType = c("gmt","RDS"),
  ...
  ){

  if(maxGap == 200 & minCpGs == 3){

    warning(
      paste("A file of close by CpGs for maxGap = 200 and minCpGs = 3
            already exist at /inst/extdata/",
            genomicRegionType, "3_200.rda", sep = "")
      )

  } else {

    genomicRegionType <- match.arg(genomicRegionType)
    arrayType <- match.arg(arrayType)

    ### Read gmtFile ###
    fileName <- paste(genomicRegionType, "Ind.gmt", sep = "")
    data_path <- system.file(
      "extdata",fileName, package = 'coMethDMR', mustWork = TRUE
    )
    gmtFile <- read_gmt(data_path)

    ### Extract regions with three or more CpGs ###
    gmtFile3CpGs_ls <- gmtFile$pathway [lapply(gmtFile$pathway, length) >= 3]

    ### Find close by clusters in all the regions ###
    ### 45.92571 secs for 1000 regions
    closeByRegions_ls <- lapply(
      gmtFile3CpGs_ls, CloseBySingleRegion, arrayType, maxGap, minCpGs
    )

    ### Remove 'NULL' elements of the list ###
    closeByRegionsNoNull_ls <- unlist(closeByRegions_ls, recursive = FALSE)

    ### Order CpGs in each cluster by location to name the cluster ###
    ### 8.202557 secs for 162 regions, after unlisting 1000 regions from previous step
    closeByRegionsOrderedDF_ls <- lapply(
      closeByRegionsNoNull_ls,
      FUN = OrderCpGsByLocation,
      arrayType,
      output = "dataframe"
    )

    ### Name each cluster with genomic region ###
    closeByRegionsNames_ls <- lapply(
      closeByRegionsOrderedDF_ls, FUN = NameRegion
    )

    ### Order CpGs in each cluster by location ###
    closeByRegionsOrdered_ls <- lapply(
      closeByRegionsOrderedDF_ls, function(x) x[,"cpg"]
    )

    ### Create gmt file ###
    out_CloseByRegions <- CreatePathwayCollection(
      pathways = closeByRegionsOrdered_ls,
      TERMS = unlist(closeByRegionsNames_ls)
    )

    ### Select and return output ###
    fileType <- match.arg(fileType)
    fileName <- paste(file, fileType, sep = ".")
    message(
      paste("Writing to file", fileName)
    )
    if (fileType == "RDS") {
      saveRDS(out_CloseByRegions, file = fileName)
    } else {
      write_gmt(out_CloseByRegions, file = fileName)
    }

    # END if else
  }

}
