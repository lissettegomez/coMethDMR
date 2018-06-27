setwd("M:/Consulting_Core/Pathway_methylation/work/cpg/AD/PFCanalysis_updatedJuly5/contiguous_subregions")

# pfc<-readRDS("C:/Users/lxg255/Box Sync/METHOD-METHYL-GENE-TEST/AD/DATA/pfc.raw.rds")
# 
# load("M:/Consulting_Core/Pathway_methylation/work/cpg/AD/PFCanalysis_updatedJuly5/mixed_model_beta_filtered_alpha/fullannotInd.rda")
# annot<-as.data.frame(fullannot)
# type<-annot[c("NAME","INFINIUM_DESIGN_TYPE")]

##Normalization
# set.seed(5000)
# source("BMIQ_1.3.R")
# pfc$ProbeID<-rownames(pfc)
# PFCbetas.type<-merge(type,pfc,by.x="NAME",by.y="ProbeID")
# type12<-ifelse(PFCbetas.type$INFINIUM_DESIGN_TYPE=='I',1,2)
# PFCbetas.bmiq<-as.data.frame(matrix(ncol=0,nrow=485577))
# for (i in 3:112) {betas.sample<-PFCbetas.type[,i]; PFCbetas.bmiq.sample<-BMIQ(beta.v=betas.sample,design.v=type12,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=i-2); PFCbetas.bmiq<-cbind(PFCbetas.bmiq,PFCbetas.bmiq.sample$nbeta)}
# colnames(PFCbetas.bmiq)<-colnames(PFCbetas.type)[3:112]
# rownames(PFCbetas.bmiq)<-PFCbetas.type$NAME
# write.csv(PFCbetas.bmiq,"PFCbetas_BMIQnormalized.csv")
PFCbetas.bmiq<-read.csv("PFCbetas_BMIQnormalized.csv")
rownames(PFCbetas.bmiq)<-PFCbetas.bmiq[,1]
PFCbetas.bmiq<-PFCbetas.bmiq[,-1]


source("M:/Consulting_Core/Pathway_methylation/work/cpg/functions/contiguous_regions_input_one_region_rdrop_LWtest.R")
#source("M:/Consulting_Core/Pathway_methylation/work/cpg/functions/subislands.R")
#source("M:/Consulting_Core/Pathway_methylation/work/cpg/functions/common_cpgs.R")
source("M:/Consulting_Core/Pathway_methylation/work/cpg/functions/closeCpGs_regions_LWedit.R")
#source("M:/Consulting_Core/Pathway_methylation/work/cpg/functions/alpha_iqr_range.R")
source("M:/Consulting_Core/Pathway_methylation/work/cpg/functions/cpgs_location_name.R")


location.file <- readRDS ("C:/Users/lxg255/Box Sync/METHOD-METHYL-GENE-TEST/AD/DATA/cpg.locations.RDS")
colnames(location.file)[1] <- "cpg"

cgis <- readRDS ("C:/Users/lxg255/Box Sync/METHOD-METHYL-GENE-TEST/AD/DATA/ISLANDInd.RDS")
cgi <- cgis$SID
cgi.3 <- cgi [lapply(cgi, length) >=3]
#cgi.3<-unname(cgi.3)

pfc.closeCpGs.regions<-lapply(cgi.3,closeCpGs.regions,location.file,maxGap=200,min.cpgs=3)
pfc.closeCpGs.regions.final<-unlist(pfc.closeCpGs.regions,recursive=F)
names(pfc.closeCpGs.regions.final)<-lapply(pfc.closeCpGs.regions.final,region.name,location.file)
#pfc.matrix.cgi<-lapply(pfc.closeCpGs.regions.final, function (item) PFCbetas.bmiq[item,])
pfc.matrix.cgi<-lapply(pfc.closeCpGs.regions.final, function (item) PFCbetas.bmiq[which(rownames(PFCbetas.bmiq)%in%item),])

pfc.contiguous.regions<-lapply(pfc.matrix.cgi,contiguous.regions, min.region.size=3, location.file=location.file, threshold.r=0.5)
pfc.contiguous.filtered<-lapply(pfc.contiguous.regions, function(item) item$`n.sub.regions` > 0)
pfc.contiguous.filtered.vec<-unlist(pfc.contiguous.filtered, recursive = FALSE)
pfc.contiguous.filtered.regions<-pfc.contiguous.regions[pfc.contiguous.filtered.vec]
pfc.contiguous.filtered.cpgs<-lapply(pfc.contiguous.filtered.regions, function(x) x[[3]])
pfc.contiguous.filtered.cpgs.vec<-unlist(pfc.contiguous.filtered.cpgs,recursive = FALSE)
#final.subregions.beta<-lapply(pfc.contiguous.filtered.cpgs.vec, function (item) PFCbetas.bmiq[item,])
final.subregions.beta<-lapply(pfc.contiguous.filtered.cpgs.vec, function (item) PFCbetas.bmiq[which(rownames(PFCbetas.bmiq)%in%item),])

library(lmerTest)
prop.neuron<-readRDS("M:/Consulting_Core/Pathway_methylation/work/cpg/AD/LWtest/CETS/pfc.prop.neurons.RDS")
pheno<-readRDS("C:/Users/lxg255/Box Sync/METHOD-METHYL-GENE-TEST/AD/DATA/pheno.pfc.RDS")
pheno.prop.neuron<-merge(pheno,prop.neuron,by.x="sample.id",by.y="row.names")


mxmodel<-function(beta.matrix2,pheno.prop.neuron)  {
  beta.matrix2$ProbeID<-row.names(beta.matrix2)
  beta.matrix.t<-reshape(beta.matrix2, varying = colnames(beta.matrix2[-dim(beta.matrix2)[2]]), v.names="beta", direction="long", time=colnames(beta.matrix2[-dim(beta.matrix2)[2]]),timevar="Sample")
  temp<-merge(beta.matrix.t,pheno.prop.neuron,by.x="Sample",by.y="sample.id")
  temp$Mvalue<-log2(temp$beta/(1-temp$beta))
  ps=999
  tryCatch({
    f<-lmer(Mvalue ~stage + (stage|ProbeID) +(1|Sample) + age.brain + sex + as.factor(Mplate)+ prop.neuron, temp)
    ps<-coef(summary(f))[2,5]
    slope<-coef(summary(f))[2,1]
  }, error=function(e){})
  correlation<-cor(t(beta.matrix2[-dim(beta.matrix2)[2]]),method="spearman", use="pairwise.complete.obs")
  median.corr<-median(correlation[lower.tri(correlation)])
  l<-list(slope,ps,median.corr)
  names(l)<-c("slope","Pval.mixed.model","median.corr")
  return(l)
}


mxmodel.pvalue<-lapply(final.subregions.beta, mxmodel, pheno.prop.neuron=pheno.prop.neuron)


split.name<-lapply(names(mxmodel.pvalue), strsplit, split="[.]")
split.name.unlist<-lapply(split.name, unlist)
region.id<-lapply(split.name.unlist, function(x) x[[1]])
#close.subregion.id<-lapply(pfc.closeCpGs.regions.final,region.name,location.file)
contiguous.subregion.id<-lapply(split.name.unlist, function(x) x[[2]])
results<-data.frame("region.id"=unlist(region.id),
                    "contiguous.subregion.id"=unlist(contiguous.subregion.id),
                    "slope"=unlist(unname(lapply(mxmodel.pvalue, function(item) item$`slope`))),
                    "Pval.mixed.model"=unlist(unname(lapply(mxmodel.pvalue, function(item) item$`Pval.mixed.model`))),
                    "median.corr"=unlist(unname(lapply(mxmodel.pvalue, function(item) item$`median.corr`))))

write.csv(results,"ADresultsRstudio_spearman_corr_slope.csv")




##analysis by sex

##contiguous regions extracted from female and male samples, then mixed model on female and male separately

pheno.prop.neuron.female<-subset(pheno.prop.neuron,sex=="Sex: FEMALE")
pheno.prop.neuron.male<-subset(pheno.prop.neuron,sex=="Sex: MALE")


mxmodel.bysex<-function(beta.matrix2,pheno.prop.neuron)  {
  beta.matrix2$ProbeID<-row.names(beta.matrix2)
  beta.matrix.t<-reshape(beta.matrix2, varying = colnames(beta.matrix2[-dim(beta.matrix2)[2]]), v.names="beta", direction="long", time=colnames(beta.matrix2[-dim(beta.matrix2)[2]]),timevar="Sample")
  temp<-merge(beta.matrix.t,pheno.prop.neuron,by.x="Sample",by.y="sample.id")
  temp$Mvalue<-log2(temp$beta/(1-temp$beta))
  ps=999
  tryCatch({
    f<-lmer(Mvalue ~stage + (stage|ProbeID) +(1|Sample) + age.brain + as.factor(Mplate)+ prop.neuron, temp)
    ps<-coef(summary(f))[2,5]
    slope<-coef(summary(f))[2,1]
  }, error=function(e){})
  beta.matrix2.bysex<-
  correlation<-cor(t(beta.matrix2[,as.character(pheno.prop.neuron$sample.id)]),method="spearman", use="pairwise.complete.obs")
  median.corr<-median(correlation[lower.tri(correlation)])
  l<-list(slope,ps,median.corr)
  names(l)<-c("slope","Pval.mixed.model","median.corr")
  return(l)
}



mxmodel.pvalue.female<-lapply(final.subregions.beta, mxmodel.bysex, pheno.prop.neuron=pheno.prop.neuron.female)
mxmodel.pvalue.male<-lapply(final.subregions.beta, mxmodel.bysex, pheno.prop.neuron=pheno.prop.neuron.male)




split.name<-lapply(names(mxmodel.pvalue.female), strsplit, split="[.]")
split.name.unlist<-lapply(split.name, unlist)
region.id<-lapply(split.name.unlist, function(x) x[[1]])
subregion.id<-lapply(split.name.unlist, function(x) x[[2]])
results.female<-data.frame("region.id"=unlist(region.id),
                    "subregion.id"=unlist(subregion.id),
                    "slope"=unlist(unname(lapply(mxmodel.pvalue.female, function(item) item$`slope`))),
                    "Pval.mixed.model"=unlist(unname(lapply(mxmodel.pvalue.female, function(item) item$`Pval.mixed.model`))),
                    "median.corr"=unlist(unname(lapply(mxmodel.pvalue.female, function(item) item$`median.corr`))))

write.csv(results.female,"ADresults.female.slope.csv")



split.name<-lapply(names(mxmodel.pvalue.male), strsplit, split="[.]")
split.name.unlist<-lapply(split.name, unlist)
region.id<-lapply(split.name.unlist, function(x) x[[1]])
subregion.id<-lapply(split.name.unlist, function(x) x[[2]])
results.male<-data.frame("region.id"=unlist(region.id),
                    "subregion.id"=unlist(subregion.id),
                    "slope"=unlist(unname(lapply(mxmodel.pvalue.male, function(item) item$`slope`))),
                    "Pval.mixed.model"=unlist(unname(lapply(mxmodel.pvalue.male, function(item) item$`Pval.mixed.model`))),
                    "median.corr"=unlist(unname(lapply(mxmodel.pvalue.male, function(item) item$`median.corr`))))

write.csv(results.male,"ADresults.male.slope.csv")



##extract contiguous regions on female samples only 

female.samples<-as.character(pheno.prop.neuron.female$sample.id)
PFCbetas.bmiq.female<-PFCbetas.bmiq[,female.samples]

pfc.matrix.cgi.female<-lapply(pfc.closeCpGs.regions.final, function (item) PFCbetas.bmiq.female[which(rownames(PFCbetas.bmiq.female)%in%item),])

pfc.contiguous.regions<-lapply(pfc.matrix.cgi.female,contiguous.regions, min.region.size=3, location.file=location.file, threshold.r=0.5)
pfc.contiguous.filtered<-lapply(pfc.contiguous.regions, function(item) item$`n.sub.regions` > 0)
pfc.contiguous.filtered.vec<-unlist(pfc.contiguous.filtered, recursive = FALSE)
pfc.contiguous.filtered.regions<-pfc.contiguous.regions[pfc.contiguous.filtered.vec]
pfc.contiguous.filtered.cpgs<-lapply(pfc.contiguous.filtered.regions, function(x) x[[3]])
pfc.contiguous.filtered.cpgs.vec<-unlist(pfc.contiguous.filtered.cpgs,recursive = FALSE)
#final.subregions.beta<-lapply(pfc.contiguous.filtered.cpgs.vec, function (item) PFCbetas.bmiq[item,])
final.subregions.beta.female<-lapply(pfc.contiguous.filtered.cpgs.vec, function (item) PFCbetas.bmiq[which(rownames(PFCbetas.bmiq)%in%item),])

mxmodel.pvalue.female<-lapply(final.subregions.beta.female, mxmodel.bysex, pheno.prop.neuron=pheno.prop.neuron.female)

split.name<-lapply(names(mxmodel.pvalue.female), strsplit, split="[.]")
split.name.unlist<-lapply(split.name, unlist)
region.id<-lapply(split.name.unlist, function(x) x[[1]])
subregion.id<-lapply(split.name.unlist, function(x) x[[2]])
results<-data.frame("region.id"=unlist(region.id),
                    "subregion.id"=unlist(subregion.id),
                    "Pval.mixed.model"=unlist(unname(lapply(mxmodel.pvalue.female, function(item) item$`Pval.mixed.model`))),
                    "median.corr"=unlist(unname(lapply(mxmodel.pvalue.female, function(item) item$`median.corr`))))

write.csv(results,"ADresults.female.contiguous.regions.csv")


##extract contiguous regions on male samples only

male.samples<-as.character(pheno.prop.neuron.male$sample.id)
PFCbetas.bmiq.male<-PFCbetas.bmiq[,male.samples]

pfc.matrix.cgi.male<-lapply(pfc.closeCpGs.regions.final, function (item) PFCbetas.bmiq.male[which(rownames(PFCbetas.bmiq.male)%in%item),])

pfc.contiguous.regions<-lapply(pfc.matrix.cgi.female,contiguous.regions, min.region.size=3, location.file=location.file, threshold.r=0.5)
pfc.contiguous.filtered<-lapply(pfc.contiguous.regions, function(item) item$`n.sub.regions` > 0)
pfc.contiguous.filtered.vec<-unlist(pfc.contiguous.filtered, recursive = FALSE)
pfc.contiguous.filtered.regions<-pfc.contiguous.regions[pfc.contiguous.filtered.vec]
pfc.contiguous.filtered.cpgs<-lapply(pfc.contiguous.filtered.regions, function(x) x[[3]])
pfc.contiguous.filtered.cpgs.vec<-unlist(pfc.contiguous.filtered.cpgs,recursive = FALSE)
#final.subregions.beta<-lapply(pfc.contiguous.filtered.cpgs.vec, function (item) PFCbetas.bmiq[item,])
final.subregions.beta.male<-lapply(pfc.contiguous.filtered.cpgs.vec, function (item) PFCbetas.bmiq[which(rownames(PFCbetas.bmiq)%in%item),])

mxmodel.pvalue.male<-lapply(final.subregions.beta.male, mxmodel.bysex, pheno.prop.neuron=pheno.prop.neuron.male)

split.name<-lapply(names(mxmodel.pvalue.female), strsplit, split="[.]")
split.name.unlist<-lapply(split.name, unlist)
region.id<-lapply(split.name.unlist, function(x) x[[1]])
subregion.id<-lapply(split.name.unlist, function(x) x[[2]])
results<-data.frame("region.id"=unlist(region.id),
                    "subregion.id"=unlist(subregion.id),
                    "Pval.mixed.model"=unlist(unname(lapply(mxmodel.pvalue.male, function(item) item$`Pval.mixed.model`))),
                    "median.corr"=unlist(unname(lapply(mxmodel.pvalue.male, function(item) item$`median.corr`))))

write.csv(results,"ADresults.male.contiguous.regions.csv")





###########################################################################

### individual CpGs analysis for most significant regions
# region1<-final.subregions.beta$`Chr22:19709548-19709755.Chr22:19709548-19709755`
# region2<-final.subregions.beta$`Chr2:241721922-241722113.Chr2:241721922-241722113`
# region3<-final.subregions.beta$`Chr7:27146237-27146445.Chr7:27146237-27146445`
# region4<-final.subregions.beta$`Chr21:47855893-47856140.Chr21:47855893-47856020`
# region5<-final.subregions.beta$`Chr16:29675846-29676071.Chr16:29675846-29676071`
# region6<-final.subregions.beta$`Chr17:56355299-56355431.Chr17:56355299-56355431`

top10f.regions<-read.table("top10regions_females.txt")
top10f.regions$V3<-paste(top10f.regions$V1,".",top10f.regions$V2,sep="")

# pdf("linear_correlation_plots_6top_regions.pdf")
# region<-list(region1,region2,region3,region4,region5,region6)

pdf("linear_correlation_plots_10topf_regions.pdf")

n=1
# for (i in region[1:6]){

for (i in top10f.regions[1:10,3]){
  
  
## individual CpGs pvalues  
res<-data.frame(matrix(ncol=2,nrow=0))
res.annot<-data.frame(matrix(ncol=8,nrow=0))
# subregion<-as.data.frame(i)
subregion<-as.data.frame(final.subregions.beta[[i]])
subregion<-subregion[,female.samples]
correlation<-cor(t(subregion),method="spearman", use="pairwise.complete.obs")
mean.corr<-mean(correlation[lower.tri(correlation)])
subregion$ProbeID<-row.names(subregion)
subregion.t<-reshape(subregion, varying = colnames(subregion[-dim(subregion)[2]]), v.names="beta", direction="long", time=colnames(subregion[-dim(subregion)[2]]),timevar="Sample")
temp<-merge(subregion.t,pheno.prop.neuron.female,by.x="Sample",by.y="sample.id")
temp$Mvalue<-log2(temp$beta/(1-temp$beta))
temp.loc<-merge(temp,annot[c("NAME","CHR","MAPINFO")],by.x="ProbeID",by.y="NAME")
temp.loc.ordered<-temp.loc[order(as.numeric(as.character(temp.loc$MAPINFO))),]
probe<-unique(temp.loc.ordered$ProbeID)
for (j in 1:length(probe)) { tmp<-coef(summary(lm(Mvalue ~ stage + age.brain + as.factor(Mplate)+ prop.neuron, data=subset(temp.loc.ordered, ProbeID==probe[j]))));
                                       gene<-subset(annot,NAME==probe[j],select=c("UCSC_REFGENE_NAME","UCSC_REFGENE_GROUP","MAPINFO"));
                                       gene$region<-i
                                       gene$mean.correlation<-mean.corr
                                       res[1,]= tmp[2, c(1,4)];
                                       rownames(res)<-probe[j];
                                       res.annot<-rbind(res.annot,merge(res,gene,by="row.names"))
                                       
                                       }
colnames(res.annot) <- c("ProbeID","slope.estimate", "slope.pval", "UCSC_REFGENE_NAME","UCSC_REFGENE_GROUP","MAPINFO","region","mean.correlation")
write.csv(res.annot, paste("individual_CpGs_pvalues_top_female_region",n,".csv",sep=""))

## linear plot
#pdf(paste("linear_plot_region",n,".pdf",sep=""))
library(ggplot2)
means.cpg <-aggregate(Mvalue ~ stage + ProbeID, data=temp.loc.ordered, mean)
means.cpg$av.beta.value<-means.cpg$beta
# av <-aggregate (Mvalue ~ ProbeID, data=means.cpg, mean) 
# av$mean <-av$Mvalue
# av <-subset(av, select=-Mvalue)
# means.cpg.center <-merge (means.cpg, av, by="ProbeID")
# means.cpg.center$Mvalue.centered <-means.cpg.center$Mvalue -means.cpg.center$mean
plot<-ggplot(means.cpg, aes(x=stage, y=av.beta.value, color=ProbeID) ) + geom_point(size=3, shape=21) + geom_line() + ggtitle(paste("region", i))
print(plot)


## subregion correlation plot
#pdf(paste("correlation_plot_region",n,".pdf",sep=""))
require(corrplot)
subregion<-subregion[probe,]
temp2<-subregion[-dim(beta.matrix2)[2]]
temp3<-t(temp2)
corr.spearman <-cor (temp3, method="spearman", use="pairwise.complete.obs")
print (corrplot(corr.spearman, method="number", number.cex = 1))
n=n+1

}
dev.off()


#####
### individual CpGs analysis for most significant regions in females
library(limma)


region1f<-final.subregions.beta$`Chr18:74799495-74799572.Chr18:74799495-74799572`
region2f<-final.subregions.beta$`Chr19:46270244-46270412.Chr19:46270244-46270412`


female.samples<-as.character(pheno.prop.neuron.female$sample.id)
region1f.female<-region1f[,female.samples]
region2f.female<-region2f[,female.samples]

male.samples<-as.character(pheno.prop.neuron.male$sample.id)
region1f.male<-region1f[,male.samples]
region2f.male<-region2f[,male.samples]

mod.f = model.matrix(~stage, data = pheno.prop.neuron.female)
mod.m = model.matrix(~stage, data = pheno.prop.neuron.male)
region1f.female.adj<-removeBatchEffect(region1f.female,batch=as.factor(pheno.prop.neuron.female$Mplate), covariates=pheno.prop.neuron.female[,c("age.brain","prop.neuron")], design= mod.f)
region1f.male.adj<-removeBatchEffect(region1f.male,batch=as.factor(pheno.prop.neuron.male$Mplate), covariates=pheno.prop.neuron.male[,c("age.brain","prop.neuron")], design= mod.m)

region2f.female.adj<-removeBatchEffect(region2f.female,batch=as.factor(pheno.prop.neuron.female$Mplate), covariates=pheno.prop.neuron.female[,c("age.brain","prop.neuron")], design= mod.f)
region2f.male.adj<-removeBatchEffect(region2f.male,batch=as.factor(pheno.prop.neuron.male$Mplate), covariates=pheno.prop.neuron.male[,c("age.brain","prop.neuron")], design= mod.m)




subregion.f<-as.data.frame(region1f.female.adj)
subregion.f<-as.data.frame(region2f.female.adj)

subregion.f$ProbeID<-row.names(subregion.f)
subregion.f.t<-reshape(subregion.f, varying = colnames(subregion.f[-dim(subregion.f)[2]]), v.names="beta", direction="long", time=colnames(subregion.f[-dim(subregion.f)[2]]),timevar="Sample")

subregion.m<-as.data.frame(region1f.male.adj)
subregion.m<-as.data.frame(region2f.male.adj)

subregion.m$ProbeID<-row.names(subregion.m)
subregion.m.t<-reshape(subregion.m, varying = colnames(subregion.m[-dim(subregion.m)[2]]), v.names="beta", direction="long", time=colnames(subregion.m[-dim(subregion.m)[2]]),timevar="Sample")


temp.f<-merge(subregion.f.t,pheno.prop.neuron.female,by.x="Sample",by.y="sample.id")
temp.m<-merge(subregion.m.t,pheno.prop.neuron.male,by.x="Sample",by.y="sample.id")
#temp$Mvalue<-log2(temp$beta/(1-temp$beta))

temp.f.loc<-merge(temp.f,annot[c("NAME","CHR","MAPINFO")],by.x="ProbeID",by.y="NAME")
temp.m.loc<-merge(temp.m,annot[c("NAME","CHR","MAPINFO")],by.x="ProbeID",by.y="NAME")

temp.f.loc.ordered<-temp.f.loc[order(as.numeric(as.character(temp.f.loc$MAPINFO))),]
temp.m.loc.ordered<-temp.m.loc[order(as.numeric(as.character(temp.m.loc$MAPINFO))),]

means.cpg.f <-aggregate(beta ~ stage + ProbeID , data=temp.f.loc.ordered, mean)
means.cpg.f$sex<-"female"
means.cpg.m <-aggregate(beta ~ stage + ProbeID , data=temp.m.loc.ordered, mean)
means.cpg.m$sex<-"male"
#av <-aggregate (beta ~ ProbeID, data=means.cpg, mean) 
#av$mean <-av$beta
#av <-subset(av, select=-beta)
#means.cpg.center <-merge (means.cpg, av, by="ProbeID")
#means.cpg.center$beta.centered <-means.cpg.center$beta -means.cpg.center$mean
means.cpg<-rbind(means.cpg.f,means.cpg.m)
means.cpg$group<-paste(means.cpg$ProbeID,"_",means.cpg$sex,sep="")
means.cpg$av.beta.value<-means.cpg$beta
plot<-ggplot(means.cpg, aes(x=stage, y=av.beta.value, group=group) ) + geom_smooth (se = FALSE) + ggtitle("Chr18:74799495-74799572") + aes(color=sex)
print(plot)


plot<-ggplot(means.cpg, aes(x=stage, y=av.beta.value, group=group) ) + geom_point() + stat_smooth (method = "lm", formula = y ~ x, size = 1) + ggtitle("Chr18:74799495-74799572") + aes(color=sex)+ theme_bw() + theme(plot.title=element_text(size=10)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
