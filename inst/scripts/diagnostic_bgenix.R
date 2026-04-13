# DIAGNOSTIC: Manual bgenix test for chr 1 (works) vs chr 2 (fails?)
# Run this AFTER the Conti data prep block (so Conti2 and bgenix exist)

library(readr)
library(dplyr)

Conti2 <- read_tsv('Conti.tsv', show_col_types=FALSE)

# ---- CHR 1 (known working) ----
chr1 <- Conti2 %>% filter(CHR == "1")
cat("=== CHR 1: ", nrow(chr1), "variants ===\n")
cat("First 3 RSIDs:", paste(head(chr1$rsID, 3), collapse=", "), "\n")
cat("First 3 positions:", paste(head(chr1$POS, 3), collapse=", "), "\n")

# ---- CHR 2 (likely failing) ----
chr2 <- Conti2 %>% filter(CHR == "2")
cat("\n=== CHR 2: ", nrow(chr2), "variants ===\n")
cat("First 3 RSIDs:", paste(head(chr2$rsID, 3), collapse=", "), "\n")
cat("First 3 positions:", paste(head(chr2$POS, 3), collapse=", "), "\n")

# ---- Test 1: Does .bgi index exist for chr 1 and chr 2? ----
cat("\n=== TEST 1: .bgi index files ===\n")
system("ls -la '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.bgen.bgi' 2>&1")
system("ls -la '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.bgen.bgi' 2>&1")

# ---- Test 2: Manual bgenix RSID extraction for chr 2 ----
writeLines(chr2$rsID, "_diag_rsids.txt")
cat("\n=== TEST 2: bgenix -incl-rsids for chr 2 ===\n")
system("~/_ukbrapr_tools/bgenix -g '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.bgen' -incl-rsids _diag_rsids.txt -list 2>&1 | head -20")
cat("(line count):\n")
system("~/_ukbrapr_tools/bgenix -g '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.bgen' -incl-rsids _diag_rsids.txt -list 2>&1 | wc -l")

# ---- Test 3: Manual bgenix range extraction for chr 2 ----
ranges <- paste0("02:", chr2$POS, "-", chr2$POS)
writeLines(ranges, "_diag_ranges.txt")
cat("\n=== TEST 3: bgenix -incl-range for chr 2 ===\n")
cat("Range file first 3 lines:\n")
system("head -3 _diag_ranges.txt")
system("~/_ukbrapr_tools/bgenix -g '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.bgen' -incl-range _diag_ranges.txt -list 2>&1 | head -20")
cat("(line count):\n")
system("~/_ukbrapr_tools/bgenix -g '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.bgen' -incl-range _diag_ranges.txt -list 2>&1 | wc -l")

# ---- Test 4: Check MFI for a chr 2 RSID and position ----
cat("\n=== TEST 4: MFI lookup for first chr 2 variant ===\n")
cat("Looking for RSID:", chr2$rsID[1], "\n")
system(paste0("grep -w '", chr2$rsID[1], "' '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.mfi.txt' 2>&1 | head -3"))
cat("Looking for position:", chr2$POS[1], "\n")
system(paste0("grep -w '", chr2$POS[1], "' '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.mfi.txt' 2>&1 | head -3"))

# ---- Test 5: For comparison, same tests on chr 1 (known working) ----
writeLines(chr1$rsID, "_diag_rsids_chr1.txt")
cat("\n=== TEST 5: bgenix -incl-rsids for chr 1 (comparison) ===\n")
system("~/_ukbrapr_tools/bgenix -g '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.bgen' -incl-rsids _diag_rsids_chr1.txt -list 2>&1 | wc -l")
cat("MFI for first chr 1 RSID:", chr1$rsID[1], "\n")
system(paste0("grep -w '", chr1$rsID[1], "' '/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.mfi.txt' 2>&1 | head -3"))

# ---- Test 6: Check what files exist in imputed directory ----
cat("\n=== TEST 6: List BGEN + BGI files ===\n")
system("ls '/mnt/project/Bulk/Imputation/UKB imputation from genotype/' | grep -E '\\.(bgen|bgi)$' | head -30")

# ---- Cleanup ----
system("rm -f _diag_rsids.txt _diag_ranges.txt _diag_rsids_chr1.txt")

# ---- Test 7: Also check the output files from the last create_pgs run ----
cat("\n=== TEST 7: Scoring output files ===\n")
if (file.exists("Conti_multi_ethnic.pgs.sscore.vars")) {
  vars <- readLines("Conti_multi_ethnic.pgs.sscore.vars")
  cat("Variants used for scoring (", length(vars), "lines):\n")
  cat(head(vars, 25), sep="\n")
} else {
  cat("sscore.vars file not found\n")
}

cat("\n=== TEST 7b: First 10 lines of varlist sent to plink2 ===\n")
if (file.exists("Conti_multi_ethnic.pgs.varlist.txt")) {
  vl <- read_tsv("Conti_multi_ethnic.pgs.varlist.txt", show_col_types=FALSE)
  cat("Varlist rows:", nrow(vl), " columns:", paste(names(vl), collapse=", "), "\n")
  print(head(vl, 10))
} else {
  cat("varlist.txt file not found\n")
}
