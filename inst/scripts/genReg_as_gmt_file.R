# Test Using .gmt Files for Genomic Region Lists
# Gabriel Odom
# 20181003

# Load the Genomic Region list
load("inst/extdata/NSHELF3_200.rda")

# Convert from factor to character
nShelf_GenReg <- lapply(NSHELF3_200, as.character)

# Load the pathwayPCA package for .gmt utility functions
library(devtools)
install_github("gabrielodom/pathwayPCA")

# Create pathway collection object as Genomic Region container
out_PC <- CreatePathwayCollection(
  pathways = nShelf_GenReg, TERMS = names(nShelf_GenReg)
)

# Write to file
write_gmt(out_PC, file = "inst/extdata/test_NShelf3_200.gmt")

# Read the Genomic Region list
read_gmt("inst/extdata/test_NShelf3_200.gmt")
