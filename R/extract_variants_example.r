library(dplyr)
library(readr)

## Old version - successfully extracts rs72725854_G AND rs72725854_T

remotes::install_github("icallumwebb/ukbrapR_Callum", force = TRUE)
varlist <- data.frame(rsid=c("rs72725854"), chr=c(8))
imputed_genotypes <- ukbrapR:::extract_variants(varlist, overwrite=TRUE, 
                                                progress=TRUE, verbose=TRUE, very_verbose=TRUE)