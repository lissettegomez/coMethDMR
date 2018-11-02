

#' Extract clusters of close by CpGs from all genomic regions
#'
#' @param genomicRegionType Type of region
#' @param arrayType Type of array, 450k or EPIC
#' @param maxGap an integer, genomic locations within maxGap from each other
#'    are placed into the same cluster
#' @param minCpGs an integer, minimum number of CpGs for resulting regions
#' @param fileType the output files can be saved as .gmt or .RDS,
#'    .gmt files can be open as flat text file, .RDS files are 50% the size of
#'    .gmt files, but cannot be open
#' @param dataDir link to the directory where the filtered files will be saved
#' @param ...
#'
#' @return a file with the close by regions
#'
#' @importFrom pathwayPCA read_gmt write_gmt CreatePathwayCollection
#'
#' @export
#'
#' @examples
#'    CloseByAllRegions(
#'      genomicRegionType = "ISLAND", arrayType = "450k", maxGap = 50,
#'      minCpGs = 3, fileType = "gmt", dataDir = "inst"
#'    )
#'
CloseByAllRegions <- function(
  genomicRegionType = c("ISLAND", "TSS1500", "TSS200", "UTR5", "EXON1", "GENEBODY", "UTR3", "NSHORE", "SSHORE", "NSHELF", "SSHELF"),
  arrayType = c("450k","EPIC"), maxGap = 200, minCpGs = 3, fileType = c("gmt","RDS"), dataDir, ...
  ){

  if(maxGap == 200 & minCpGs == 3){
    print(
      paste("A file of close by CpGs for this region already exists: /inst/extdata/", genomicRegionType,"3_200.rda",sep=""), quote = F
    )
  } else {

    genomicRegionType <- match.arg(genomicRegionType)
    arrayType <- match.arg(arrayType)

    ### Read gmtFile ###
    gmtFile <- read_gmt(paste("inst/extdata/", genomicRegionType, "Ind.gmt", sep = ""))

    ### Extract regions with three or more CpGs ###
    gmtFile3CpGs_ls <- gmtFile$pathway [lapply(gmtFile$pathway, length) >= 3]

    ### Find close by clusters in all the regions ###
    ### 45.92571 secs for 1000 regions
    closeByRegions_ls <- lapply(
      gmtFile3CpGs_ls, CloseBySingleRegion, arrayType, maxGap, minCpGs
    )

    ### Remove 'NULL' elements of the list ###
    closeByRegionsNoNull_ls <- unlist(closeByRegions_ls, recursive = F)

    ### Order CpGs in each cluster by location to name the cluster ###
    ### 8.202557 secs for 162 regions, after unlisting 1000 regions from previous step
    closeByRegionsOrderedDF_ls <- lapply(
      closeByRegionsNoNull_ls, OrderCpGsByLocation, arrayType, output = "dataframe"
    )

    ### Name each cluster with genomic region ###
    closeByRegionsNames_ls <- lapply(closeByRegionsOrderedDF_ls, NameRegion)

    ### Order CpGs in each cluster by location ###
    closeByRegionsOrdered_ls <- lapply(closeByRegionsOrderedDF_ls, function(x) x[,"cpg"])

    ### Create gmt file ###
    out_CloseByRegions <- CreatePathwayCollection(
      pathways = closeByRegionsOrdered_ls, TERMS = unlist(closeByRegionsNames_ls)
    )

    ### Select and return output ###
    fileType <- match.arg(fileType)
    path_char <- paste(dataDir, genomicRegionType, sep = "/")
    fileName <- paste(path_char, fileType, sep = ".")
    if (fileType == "RDS") {
      saveRDS(out_CloseByRegions, file = fileName)
    } else {
      write_gmt(out_CloseByRegions, file = fileName)
    }
  }

}
