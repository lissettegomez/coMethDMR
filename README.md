# coMethDMR
An unsupervised approach for identifying differentially methylated regions in Illumina arrays

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

The reference manual for coMethDMR can be downloaded from /docs/coMethDMR.pdf.

## References

A manuscript on details of the coMethDMR methodology is coming soon.
