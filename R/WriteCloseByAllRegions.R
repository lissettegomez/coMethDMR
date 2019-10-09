

#' Extract clusters of CpG probes located closely
#'
#' @param fileName Name of the RDS file where the output genomic regions will be saved.
#' @param regions GRanges of input genomic regions
#' @param arrayType Type of array, can be "450k" or "EPIC"
#' @param maxGap an integer, genomic locations within maxGap from each other
#'    are placed into the same cluster
#' @param minCpGs an integer, minimum number of CpGs for each resulting region
#' @param ... Dots for internal arguments. Currently unused.
#'
#' @return a file with the genomic regions containing CpGs located closely within each
#'  inputing pre-defined genomic region
#'
#' @details For \code{maxGap} = 200 and \code{minCpGs} = 3, we already calculated
#'    the clusters of CpGs. They are saved in folder \code{/inst/extdata/}.
#'
#' @export
#'
#' @examples
#'
#' data(regions)
#'
#' WriteCloseByAllRegions(
#'    regions = regions, arrayType = "EPIC", maxGap = 50,
#'    minCpGs = 3, fileName = "closeByRegions.rds"
#' )
#'
WriteCloseByAllRegions <- function(
  fileName,
  regions,
  arrayType = c("450k","EPIC"),
  maxGap = 200,
  minCpGs = 3,
  ...
  ){

  if(maxGap == 200 & minCpGs == 3){

    warning(
      paste("A file of close by CpGs for maxGap = 200 and minCpGs = 3
            for genic and intergenic regions already exist at /inst/extdata/")
      )

  } else {

    arrayType <- match.arg(arrayType)

    switch(
      arrayType,
      "450" = {
        loc_df <- as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations)
      },
      "EPIC" = {
        loc_df <- as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations)
      }
    )

    ### Convert input from GRanges to list of vectors of CpGs ###
    regions_df <- as.data.frame(regions)
    regions_ls <- split(regions_df, seq(nrow(regions_df)))

    regionCpGs_ls <- lapply(
      regions_ls,
      function(x)
        rownames(loc_df[which(loc_df$chr %in% as.character(x[ ,"seqnames"])
                              & loc_df$pos >= x [, "start"]
                              & loc_df$pos <= x[ , "end"]),])
    )

    ### Extract regions with three or more CpGs ###
    region3CpGs_ls <- regionCpGs_ls [lapply(regionCpGs_ls, length) >= minCpGs]

    ### Find close by clusters in all the regions ###
    ### 45.92571 secs for 1000 regions
    closeByRegions_ls <- lapply(
      region3CpGs_ls, CloseBySingleRegion, arrayType, maxGap, minCpGs
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

    names(closeByRegionsOrdered_ls) <- closeByRegionsNames_ls


    ### Select and return output ###
    message(
      paste("Writing to file", fileName)
    )

    saveRDS(closeByRegionsOrdered_ls, fileName)

  }

}
