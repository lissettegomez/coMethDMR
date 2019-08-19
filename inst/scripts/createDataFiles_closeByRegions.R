
library(coMethDMR)

## https://github.com/rforge/ima/blob/master/www/fullannotInd.rda
annotData <- load("fullannotInd.rda")

for (g in 2:length(annotData)){

  regionName <- annotData[g]
  regionData <- get(regionName)
  regionData_SID <- regionData$SID
  region3 <- regionData_SID [lapply(regionData_SID, length) >=3]
  region3_200 <- lapply(region3,
                       CloseBySingleRegion,
                       arrayType = "450k",
                       maxGap = 200,
                       minCpGs = 3)
  region3_200 <- unlist(region3_200, recursive=F)
  region3_200_df <- lapply(region3_200,
                            OrderCpGsByLocation,
                            arrayType = c("450k"),
                            output = "dataframe")
  names(region3_200) <- lapply(region3_200_df, NameRegion)
  saveRDS(region3_200, paste0(gsub("Ind","",regionName),"_3_200.rds"))

}


