## coMethDMR: Accurate identification of co-methylated and differentially methylated regions in epigenome-wide association studies 

by Lissette Gomez, Gabriel J. Odom, Juan I. Young, Eden R. Martin, Lizhong Liu, Xi Chen, Anthony J. Griswold, Zhen Gao, Lanyu Zhang and Lily Wang

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
Manuscirpt is available at https://www.biorxiv.org/content/10.1101/615427v1
