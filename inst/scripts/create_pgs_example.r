### YOU WILL NEED TO RUN THESE TWO LINES SEPARETELY (AND INDIVIDUALLY), THEN THE REST OF THE SCRIPT
remotes::install_github("icallumwebb/ukbrapR_Callum", force = TRUE)

library(dplyr)

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

library(readxl)
library(stringr)
library(tidyr)

## GRSs trained to predict general prostate cancer diagnosis

system('dx download Callum/ContiGWAS/Conti2021supplementarytables.xlsx') # This is the supplementary table file from Conti et al. (2021), downloadable at: https://www.nature.com/articles/s41588-020-00748-0
system('dx download Callum/WangGWAS/Wang2023supplementarytables.xlsx') # This is the supplementary table file from Wang et al. (2023), downloadable at: https://pmc.ncbi.nlm.nih.gov/articles/PMC10841479/
system('dx download Callum/SchumacherGWAS/Schumacher.txt') # This is the list of SNPs and weights from Schumacher et al. (2018), downloadable at: https://www.pgscatalog.org/publication/PGP000019/ 
system('dx download Callum/SchumacherGWAS/BARCODE1.txt') # This is the list of SNPs and weights from BARCODE1 (2021), downloadable at: https://www.pgscatalog.org/publication/PGP000726/ 

## GRSs trained to predict aggressive prostate cancer diagnosis

system('dx download Callum/SeibertGWAS/Seibert.txt') # This is the list of SNPs and weights from Seibert et al. (2018), downloadable at: https://www.pgscatalog.org/publication/PGP000047/ 
system('dx download Callum/SeibertGWAS/Pagadala.txt') # This is the list of SNPs and weights from Pagadala et al. (2022), downloadable at: https://www.pgscatalog.org/publication/PGP000400/

############################
# Define withdrawal filter #
############################

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite') ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
  df2 = df %>% left_join(
    read_csv("withdrawn_20260310.csv", col_names = FALSE, show_col_types = FALSE) %>%
      mutate(w=1) %>%
      rename(eid = X1),
    by='eid'
  ) %>%
    filter(is.na(w))
  return(df2)
}

#####################################
# Option #1 - Conti Multiethnic GRS #
#####################################


# The below block of code reads the conti xlsx file, and puts it into a form that ukbrapR can use to make a GRS (the outputted tsv)

#If you change the effect_weight column it should be easy to run a different GRS

Conti <- read_excel("Conti2021supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
Conti=Conti[1:269,]
Conti=arrange(Conti,Chromosome,Position)
Conti=Conti%>%rename(
  rsID=`rs*`,
  CHR=Chromosome,
  POS=Position,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`Multiethnic Analysis`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Conti2=Conti%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Conti2$CHR[is.na(Conti2$CHR)]="X"
write.table(Conti2,'Conti.tsv',quote=FALSE,sep='\t',row.names = FALSE)

# Check a few sample RSIDs against the MFI to verify format
head(Conti2$rsID, 10)  # verify these look like "rs1234567"
head(Conti2$POS, 10)   # verify positions are numeric, not scientific notation

# DIAGNOSTIC: Check CHR distribution and types before writing TSV
cat("CHR column class:", class(Conti2$CHR), "\n")
cat("CHR distribution:\n")
print(table(Conti2$CHR, useNA="always"))
cat("Sample of CHR values:", paste(head(unique(Conti2$CHR), 10), collapse=", "), "\n")
cat("Any NA in CHR?", sum(is.na(Conti2$CHR)), "\n")
cat("Any NA in POS?", sum(is.na(Conti2$POS)), "\n")
cat("Any NA in rsID?", sum(is.na(Conti2$rsID)), "\n")
# Check for hidden whitespace/special chars in RSIDs
cat("RSID nchar range:", range(nchar(Conti2$rsID, type="bytes"), na.rm=TRUE), "\n")
cat("All RSIDs start with 'rs'?", all(grepl("^rs", Conti2$rsID, perl=TRUE), na.rm=TRUE), "\n")

# ALSO: read the TSV back and check what readr sees (this is what create_pgs will read)
test_read <- readr::read_tsv('Conti.tsv', show_col_types=TRUE)
cat("After re-reading TSV - CHR class:", class(test_read$CHR), "\n")
cat("After re-reading TSV - CHR distribution:\n")
print(table(test_read$CHR, useNA="always"))

conti_out=ukbrapR:::create_pgs(
  in_file='Conti.tsv',
  out_file='Conti_multi_ethnic.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Conti_multi_ethnic.pgs.diagnostic_extracted.txt', header=TRUE)  # in BGEN workflow, check the diagnostic file instead of .bim
#outscore=read.table('Conti_multi_ethnic.pgs.profile',header=T)

GRS <- read.table("Conti_multi_ethnic.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "Conti_multi_ethnic.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Conti_multi_ethnic.pgs.tsv"))