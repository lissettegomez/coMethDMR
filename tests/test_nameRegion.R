# Test CoMethAllRegions()
# Gabriel Odom
# 2018-10-31

library(coMethDMR)

data(CpGsOrdered_df)
out_char <- NameRegion(CpGsOrdered_df)
stopifnot(out_char == "chr10:100028236-100028499")
