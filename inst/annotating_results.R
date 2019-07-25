# Annotate Results Function
# Gabriel Odom
# 2019-07-23


###  Linear Mixed Model Output  ###
# Given the lmmTest() results list in coMethDMR/vignettes/, extract and bind
#   the results into a data frame

resultsAll_ls <- readRDS("vignettes/BiocParallel_AllRegions_Out.RDS")

lmmResults_ls <- lapply(resultsAll_ls, function(regionType){

  out_ls <- lapply(regionType, `[[`, "modelFits_df")
  do.call(rbind, out_ls)

})

lmmResults_df <- do.call(rbind, lmmResults_ls)
row.names(lmmResults_df) <- NULL


###  Illumina Annotation and "Other" Info  ###
# 450k
anno_df <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
otherInfo_df <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
# EPIC
# anno_df <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations
# otherInfo_df <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other

anno_df <- as.data.frame(anno_df)
anno_df$cpg <- row.names(anno_df)

###  Add Location Rows to Ranges from Results  ###
# For each row in lmmResults_df, we want to find the location of the probes in
#   that region. Also, we want to add the UCSC information.
AnnotateRow <- function(row_df, loc_df, other_df){
  # browser()

  ###  Filter Data Frames  ###
  # Extract Row Region
  chr   <- row_df$chrom
  start <- as.integer(row_df$start)
  end   <- as.integer(row_df$end)

  # Find Probes in that Region
  chr_df  <- loc_df[which(loc_df$chr == chr), ]
  rownames(chr_df) <- NULL
  inRegion_idx <- which(chr_df$pos >= start & chr_df$pos <= end)
  out_df <- chr_df[inRegion_idx, ]

  # EDIT: this is hella wrong. The position and probe IDs are totally different
  #   Position is where on the genome we measure, but the probe ID is just
  #   gobbledygook. It's a coincidence that they are both 8 characters.
  # # Transform Probe IDs
  # # Add leading 0s until the number is 8 characters wide; append leading "cg"
  # probes_char <- sprintf("%08d", sort(out_df$pos))
  probes_int  <- as.integer(gsub("cg", "", out_df$cpg))
  probes_char <- sprintf("%08d", sort(probes_int))
  probes_char <- paste0("cg", probes_char)
  # EDIT 2: Lily said that we don't need to order the probe IDs


  # Find UCSC Annotation Information for those Probes
  interestingColumns_char <- c(
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Accession",
    "UCSC_RefGene_Group"
  )
  otherOut_df <- other_df[probes_char, interestingColumns_char]


  ###  Wrangle UCSC Annotation  ###
  refGeneGroup_char <- unlist(
    sapply(
      otherOut_df$UCSC_RefGene_Group,
      strsplit, ";",
      USE.NAMES = FALSE
    )
  )
  refGeneGroup_char <- sort(unique(refGeneGroup_char))

  refGeneName_char <- unlist(
    sapply(
      otherOut_df$UCSC_RefGene_Name,
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
  row_df$UCSC_RefGene_Group <- paste0(unique(refGeneGroup_char), collapse = ";")
  row_df$UCSC_RefGene_Name <- paste0(unique(refGeneName_char), collapse = ";")
  row_df$probes <- paste0(unique(probes_char), collapse = ";")

  row_df

}

# Test
AnnotateRow(
  row_df = lmmResults_df[1, ],
  loc_df = anno_df,
  other_df = otherInfo_df
)

CpGsInfoOneRegion(
  regionName_char = "chr22:17082770-17082787",
  betas_df = betaMatrixChr22_df,
  pheno_df, contPheno_char = "stage",
  covariates_char = c("age.brain", "sex"),
  arrayType = "450k"
)


###  Apply over Results  ###
resultsAnno_ls <- lapply(seq_len(nrow(lmmResults_df)), function(row){

  AnnotateRow(
    row_df = lmmResults_df[row, ],
    loc_df = anno_df,
    other_df = otherInfo_df
  )

})

resultsAnno_df <- do.call(rbind, resultsAnno_ls)

res_dir <- "~/Dropbox (BBSR)/GabrielOdom/coMethDMR/vignette_parallel_computing/"
write.csv(
  resultsAnno_df,
  paste0(res_dir, "test_annotated_results_20190724.csv")
)
