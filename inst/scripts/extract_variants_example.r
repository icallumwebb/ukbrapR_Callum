library(dplyr)
library(readr)

.rs.restartR() # run this line first to restart R before installing old version of the package
remove.packages("ukbrapR")
remotes::install_github("icallumwebb/ukbrapR_Callum", force = TRUE)
varlist <- data.frame(rsid=c("rs72725854"), chr=c(8))
imputed_genotypes <- ukbrapR:::extract_variants(varlist, overwrite=TRUE, 
                                                progress=TRUE, verbose=TRUE, very_verbose=TRUE)

## Old version - successfully extracts rs72725854_G AND rs72725854_T

.rs.restartR() # run this line first to restart R before installing old version of the package
remove.packages("ukbrapR")
remotes::install_github("lcpilling/ukbrapR@v0.3.9", force = TRUE)
varlist <- data.frame(rsid=c("rs72725854"), chr=c(8))
imputed_genotypes_old <- ukbrapR:::extract_variants(varlist, overwrite=TRUE, 
                                                progress=TRUE, verbose=TRUE, very_verbose=TRUE)
