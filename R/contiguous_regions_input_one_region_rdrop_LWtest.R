# Title - extract contiguous co-methylated regions

## INPUT
# beta.matrix - a matrix of beta values, with rownames = cpg ids, column names = sample ids
# min.region.size - min number of cpgs in the sub-regions
# location.file - a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg
# threshold.r - min correlation between a cpg in the region with the rest of the CpGs


## OUTPUT
# A list of three items:
# 1. n.sub-regions - number of contiguous co-methylated sub-regions
# 2. contiguous.regions - a data frame of probe location (ProbeID, CHR, MAPINFO),
#    index (ind), r.drop (correlation between the cpg with rest of the CpGs),
#    indicator for contiguous co-methylated CpGs (keep, keep.contiguous)
# 3. CpGs.subregions - lists of CpGs in each contiguous co-methylated regions

############################################

#' Extract contiguous co-methylated regions
#'
#' @param beta.matrix a matrix of beta values, with rownames = cpg ids, column names = sample ids
#' @param min.region.size min number of cpgs in the sub-regions
#' @param location.file a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg
#' @param threshold.r min correlation between a cpg in the region with the rest of the CpGs
#'
#' @return A list of three items:
#' 1. n.sub-regions - number of contiguous co-methylated sub-regions
#' 2. contiguous.regions - a data frame of probe location (ProbeID, CHR, MAPINFO),
#'    index (ind), r.drop (correlation between the cpg with rest of the CpGs),
#'    indicator for contiguous co-methylated CpGs (keep, keep.contiguous)
#' 3. CpGs.subregions - lists of CpGs in each contiguous co-methylated regions
#'
#' @importFrom psych alpha
#' @importFrom bumphunter getSegments
#'
#'
#' @export
#'
#' @examples
#'
#'
#'
contiguous.regions <-function (beta.matrix, min.region.size, location.file, threshold.r) {

  #library(psych)
  #library(bumphunter)

    tryCatch({
      ## order the CpGs in the beta.matrix by location
      cluster.location<-location.file[which(location.file$cpg%in%rownames(beta.matrix)),]
      cluster.location.ordered<-cluster.location[order(cluster.location$CHR,as.numeric(as.character(cluster.location$MAPINFO))),]
      beta.matrix.ordered<-beta.matrix[as.character(cluster.location.ordered$cpg),]
      beta.matrix.cluster<-t(beta.matrix.ordered)

      ## calculate alpha
      alphaaa<-alpha(beta.matrix.cluster, warnings = FALSE)

      # drop CpGs with r.drop < threshold.r
      alpha.cpg<-row.names(subset(alphaaa$item.stats, alphaaa$item.stats$r.drop < threshold.r))  ###drop these cpgs
      cpgs<-as.data.frame(rownames(alphaaa$alpha.drop))
      colnames(cpgs)<-"cpg"
      cpgs$alpha<-ifelse(row.names(alphaaa$alpha.drop)%in%alpha.cpg,0,1)   ##(drop=0, keep=1)
      cpgs$ind<-1:dim(beta.matrix.cluster)[2]

      ## get contiguous regions of CpGs
      contiguous.region<-getSegments(cpgs$alpha,cutoff=1)

      l<-data.frame(matrix(ncol=1,nrow=0))
      if (length(contiguous.region$upIndex)>0){
      for (j in 1:length(contiguous.region$upIndex)){l<-rbind(l,length(contiguous.region$upIndex[[j]]))}
      up.3<-as.numeric(rownames(subset(l,l[,1]>=min.region.size)))
      if (length(up.3)>0){
        contiguous.region.cpgs<-data.frame(matrix(ncol=2,nrow=0))
        for (u in 1:length(up.3))
        {contiguous.region.cpgs<-rbind(contiguous.region.cpgs,cbind(as.data.frame(subset(cpgs,ind%in%contiguous.region$upIndex[[up.3[u]]],select="cpg")),rep(u,length(contiguous.region$upIndex[[up.3[u]]]))))}
      } else {contiguous.region.cpgs<-cbind(as.data.frame(cpgs$cpg),rep(0,length(cpgs$cpg)))}
      } else {contiguous.region.cpgs<-cbind(as.data.frame(cpgs$cpg),rep(0,length(cpgs$cpg)))}
      colnames(contiguous.region.cpgs)<-c("ProbeID","subisland")

      ## output l1: data frame of location, alpha values and contiguous CpGs
      alpha.info<-cbind(cluster.location.ordered[,c("CHR","MAPINFO")],ind=1:dim(cluster.location.ordered)[1], r.drop=alphaaa$item.stats[,"r.drop"],keep=cpgs$alpha)
      alpha.info.keep.contiguous<-merge(alpha.info,contiguous.region.cpgs,by.x="row.names",by.y="ProbeID",all.x=T,sort=F)
      alpha.info.keep.contiguous<-alpha.info.keep.contiguous[order(alpha.info.keep.contiguous$ind),]
      colnames(alpha.info.keep.contiguous)[c(1,7)]<-c("ProbeID","keep.contiguous")
      alpha.info.keep.contiguous$keep.contiguous[is.na(alpha.info.keep.contiguous$keep.contiguous)]<-0
      alpha.info.keep.contiguous<-alpha.info.keep.contiguous[order(as.numeric(as.character(alpha.info.keep.contiguous$MAPINFO))),]
      #l1<-list(alpha.info.keep.contiguous)

      ## output l2 and l3: extract contiguous CpGs and location
      si<-sort(unique(contiguous.region.cpgs$subisland))
      l2<-list()
      l3<-list()
      # for (s in si[1:length(si)]){if (s==0) {l2[[1]]<-contiguous.region.cpgs$ProbeID; l3[[1]]<-paste("Chr",alpha.info.keep.contiguous$CHR[1],":",alpha.info.keep.contiguous$MAPINFO[1],"-",alpha.info.keep.contiguous$MAPINFO[length(alpha.info.keep.contiguous$MAPINFO)])}
      #   else {l2[[s]]<-subset(contiguous.region.cpgs,subisland==s)[,"ProbeID"];
      #         l3[[s]]<-paste("Chr",alpha.info.keep.contiguous$CHR[1],":",subset(alpha.info.keep.contiguous,keep.contiguous==s,select="MAPINFO")[1,1],"-",subset(alpha.info.keep.contiguous, keep.contiguous==s,select="MAPINFO")[length(l2[[s]]),1],sep="")}
      #   }

      for (s in si[1:length(si)]){if (s==0) {l2[[1]]<-contiguous.region.cpgs$ProbeID; names(l2) <- paste("Chr",alpha.info.keep.contiguous$CHR[1],":",alpha.info.keep.contiguous$MAPINFO[1],"-",alpha.info.keep.contiguous$MAPINFO[length(alpha.info.keep.contiguous$MAPINFO)])}
        else {l2[[s]]<- factor(subset(contiguous.region.cpgs,subisland==s)[,"ProbeID"]);
              l3[[s]]<-paste("Chr",alpha.info.keep.contiguous$CHR[1],":",subset(alpha.info.keep.contiguous,keep.contiguous==s,select="MAPINFO")[1,1],"-",subset(alpha.info.keep.contiguous, keep.contiguous==s,select="MAPINFO")[length(l2[[s]]),1],sep="")}
        }

      names (l2) <- unlist(l3, recursive = FALSE)


    }, error=function(e){})

    ## LW 3/21 - number of regions
    n.regions <- max(si)


    l4<-list(n.regions, alpha.info.keep.contiguous,l2)
    names (l4) <- c("n.sub.regions", "contiguous.regions", "CpGs.subregions")
    return(l4)
}


######## USAGE
#
# setwd("M:/Consulting_Core/Pathway_methylation/work/cpg/AD/PFCanalysis_updatedJuly5/contiguous_subregions")
#
# blood.betas.bmiq<-read.csv("blood_betas_BMIQnormalized.csv")
# rownames(blood.betas.bmiq)<-blood.betas.bmiq[,1]
# blood.betas.bmiq<-blood.betas.bmiq[,-1]
# location.file <- readRDS ("C:/Users/lxw391/Box Sync/METHODS-SHARED/METHOD-METHYL-GENE-TEST/AD/DATA/cpg.locations.RDS")
# colnames(location.file)[1] <- "cpg"
# cgis <- readRDS ("C:/Users/lxw391/Box Sync/METHODS-SHARED/METHOD-METHYL-GENE-TEST/AD/DATA/ISLANDInd.RDS")
# cgi <- cgis$SID
# cgi.3 <- cgi [lapply(cgi, length) >=30]
# cgi.3<-unname(cgi.3)
# blood.matrix.cgi<-lapply(cgi.3[1:10], function (item) blood.betas.bmiq[item,])
# beta.matrix<-blood.matrix.cgi[[2]]
# contiguous.regions (beta.matrix, min.region.size=3, location.file, threshold.r=0.5)
#
# library(corrplot)
# corr.values <- cor (beta.matrix.cluster, method = "spearman")
# corrplot(corr.values, method="number", number.cex = 0.5, tl.cex = 0.7)
