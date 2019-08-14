
library(pathwayPCA)

genomic_region <- c("ISLAND_3_200", "NSHELF_3_200", "SSHELF_3_200", "NSHORE_3_200",
                  "SSHORE_3_200", "TSS200_3_200", "TSS1500_3_200",
                  "GENEBODY_3_200", "UTR3_3_200", "UTR5_3_200", "EXON1_3_200")

for (g in genomic_region[1:length(genomic_region)]){
  region_file <- readRDS(paste0(g,".RDS"))
  out_CloseByRegions <- CreatePathwayCollection(
     sets_ls = unname(region_file),
     TERMS = names(region_file),
     setType = "regions"
  )
  write_gmt(out_CloseByRegions, file = paste0(g, ".gmt"), setType = "regions")
}
