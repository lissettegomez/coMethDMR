## coMethDMR: Accurate identification of co-methylated and differentially methylated regions in epigenome-wide association studies 

## Description
coMethDMR is an R package that identifies genomic regions associated with continuous phenotypes by optimally leverages covariations 
among CpGs within predefined genomic regions. Instead of testing all CpGs within a genomic region, coMethDMR carries out an additional 
step that selects comethylated sub-regions first without using any outcome information. Next, coMethDMR tests association between 
methylation within the sub-region and continuous phenotype using a random coefficient mixed effects model, which models both variations 
between CpG sites within the region and differential methylation simultaneously.

## Installation

The latest version can be installed by

```{r eval=FALSE, message=FALSE, warning=FALSE, results='hide'}
library(devtools)
install_github("lissettegomez/coMethDMR")
```
After installation, the coMethDMR package can be loaded into R using:

```{r eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
library(coMethDMR)
```

## Manual

The reference manual for coMethDMR can be downloaded from /docs/coMethDMR_0.0.0.9001.pdf. Two vignettes are available in the same directory: "1_Introduction_coMethDMR_10-9-2019.pdf" and "2_BiocParallel_for_coMethDMR_geneBasedPipeline.pdf"

## References
Gomez L, Odom GJ, Young JI, Martin ER, Liu L, Chen X, Griswold AJ, Gao Z, Zhang L, Wang L (2019) coMethDMR: Accurate identification of co-methylated and differentially methylated regions in epigenome-wide association studies with continuous phenotypes. Nucleic Acids Research, gkz590, https://doi.org/10.1093/nar/gkz590
