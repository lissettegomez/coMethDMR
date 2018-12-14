# Test CoMethAllRegions()
# Gabriel Odom
# 2018-10-31

library(coMethDMR)

CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
CpGsOrdered_df <- OrderCpGsByLocation(CpGs_char, arrayType=c("EPIC"), output = "dataframe")
out_char <- NameRegion(CpGsOrdered_df)
stopifnot(out_char == "chr10:100028236-100028499")
