#' Annotate \code{coMethDMR} Pipeline Results
#'
#' @description Given a data frame with regions of the genome, add UCSC
#'    Reference Gene information, probe IDs, and relationship to known CpG
#'    islands.
#'
#' @param lmmRes_df A data frame returned by \code{\link{lmmTest}} or a parallel
#'    implementation thereof. This data frame must contain the following
#'    columns: \code{chrom} = the chromosome the region is on; \code{start} =
#'    the region start point; \code{end} = the region end point; and
#'    \code{regionType} = a character string marking from which region type the
#'    region came from. See \strong{Details} for more information.
#' @param arrayType Type of array: 450k or EPIC
#'
#' @return a data frame with location of the genomic region's chromosome
#'    (\code{chrom}), start (\code{start}), and end (\code{end}); the number of
#'    CpGs in the region (\code{nCpGs}); results for testing association of
#'    methylation in individual CpGs with continuous phenotype (\code{Estimate},
#'    \code{StdErr}, \code{Stat}, and \code{pValue}); the type of region
#'    (\code{regionType}) and its relationship to known CpG islands
#'    (\code{Relation_to_UCSC_CpG_Island}); UCSC annotation information
#'    (\code{UCSC_RefGene_Group}, \code{UCSC_RefGene_Accession}, and
#'    \code{UCSC_RefGene_Name}); and a list of all of the probes in that region
#'    (\code{probes}).
#'
#' @details The region types include \code{"NSHORE"}, \code{"NSHELF"},
#'    \code{"SSHORE"}, \code{"SSHELF"}, \code{"TSS1500"}, \code{"TSS200"},
#'    \code{"UTR5"}, \code{"EXON1"}, \code{"GENEBODY"}, \code{"UTR3"}, and
#'    \code{"ISLAND"}.
#'
#' @export
#'
#' @examples
#'    lmmResults_df <- data.frame(
#'      chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'      start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'      end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'      regionType = c("TSS1500", "EXON1", "ISLAND", "TSS200", "ISLAND"),
#'      stringsAsFactors = FALSE
#'    )
#'
#'    AnnotateResults(
#'      lmmRes_df = lmmResults_df,
#'      arrayType = "450k"
#'    )
#'
AnnotateResults <- function(lmmRes_df, arrayType = c("450k","EPIC")){
  # browser()

  ###  Check Inputs  ###
  stopifnot(
    "data.frame" %in% class(lmmRes_df),
    all(c("chrom", "start", "end", "regionType") %in% colnames(lmmRes_df))
  )
  arrayType <- match.arg(arrayType)

  lmmRes_df$start <- as.integer(lmmRes_df$start)
  lmmRes_df$end   <- as.integer(lmmRes_df$end)


  ###  Pull Database  ###
  switch(
    arrayType,
    "450k" = {

      locations_df <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
      UCSCinfo_df  <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other

    },
    "EPIC" = {

      locations_df <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations
      UCSCinfo_df  <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other

    }
  )

  # Locations
  locations_df <- as.data.frame(locations_df)
  locations_df$cpg <- row.names(locations_df)
  rownames(locations_df) <- NULL

  # UCSC Gene Info
  UCSCinfo_df <- as.data.frame(UCSCinfo_df)
  interestingColumns_char <- c(
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Accession",
    "UCSC_RefGene_Group"
  )
  UCSCinfo_df <- UCSCinfo_df[, interestingColumns_char]



  ###  Define Wrapper Function  ###
  AnnotateRow <- function(row_df, loc_df, info_df){
    # browser()

    ###  Filter Data Frames  ###
    # Extract Row Region
    chr   <- row_df$chrom
    start <- row_df$start
    end   <- row_df$end

    # Find Probes in that Region
    chr_df  <- loc_df[which(loc_df$chr == chr), ]
    inRegion_idx <- which(chr_df$pos >= start & chr_df$pos <= end)
    out_df <- chr_df[inRegion_idx, ]
    probes_char <- out_df$cpg

    # Find UCSC Annotation Information for those Probes
    infoOut_df <- info_df[probes_char, ]


    ###  Wrangle UCSC Annotation  ###
    refGeneGroup_char <- unlist(
      sapply(
        infoOut_df$UCSC_RefGene_Group,
        strsplit, ";",
        USE.NAMES = FALSE
      )
    )
    refGeneGroup_char <- sort(unique(refGeneGroup_char))

    refGeneAcc_char <- unlist(
      sapply(
        infoOut_df$UCSC_RefGene_Accession,
        strsplit, ";",
        USE.NAMES = FALSE
      )
    )
    refGeneAcc_char <- sort(unique(refGeneAcc_char))

    refGeneName_char <- unlist(
      sapply(
        infoOut_df$UCSC_RefGene_Name,
        strsplit, ";",
        USE.NAMES = FALSE
      )
    )
    refGeneName_char <- sort(unique(refGeneName_char))


    ###  Return Annotated 1-Row Data Frame  ###
    row_df$Relation_to_UCSC_CpG_Island <- ifelse(
      test = row_df$regionType %in%
        c("NSHELF", "NSHORE", "ISLAND", "SSHORE", "SSHELF"),
      yes  = row_df$regionType,
      no   = ""
    )
    row_df$UCSC_RefGene_Group <-
      paste0(unique(refGeneGroup_char), collapse = ";")
    row_df$UCSC_RefGene_Accession <-
      paste0(unique(refGeneAcc_char), collapse = ";")
    row_df$UCSC_RefGene_Name <-
      paste0(unique(refGeneName_char), collapse = ";")
    row_df$probes <-
      paste0(unique(probes_char), collapse = ";")

    row_df

  }

  resultsAnno_ls <- lapply(seq_len(nrow(lmmRes_df)), function(row){

    AnnotateRow(
      row_df = lmmRes_df[row, ],
      loc_df = locations_df,
      info_df = UCSCinfo_df
    )

  })

  do.call(rbind, resultsAnno_ls)

}
