## Input
#cpgs - vector of CpGs
#location.file - a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg

## Output
#genome location of the CpGs region

#######################################

#' Title
#'
#' @param cpgs vector of CpGs
#' @param location.file a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg
#'
#' @return genome location of the CpGs region
#'
#' @export
#'
#' @examples
#'
#'
region.name<-function(cpgs,location.file){

cpgs.location<-location.file[which(location.file$cpg%in%cpgs),]
      cpgs.location.ordered<-cpgs.location[order(cpgs.location$CHR,as.numeric(as.character(cpgs.location$MAPINFO))),]
	cpgs.location.name<-paste("Chr",cpgs.location.ordered$CHR[1],":",cpgs.location.ordered$MAPINFO[1],"-",cpgs.location.ordered$MAPINFO[length(cpgs)], sep="")
	return(cpgs.location.name)
}




