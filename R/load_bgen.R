
#' Load BGEN file into memory
#'
#' @description Use plink2 to export a BGEN file to RAW format and load as a data frame
#'
#' @return A data frame
#'
#' @author Callum Webb, adapted from load_bed.R by Luke Pilling
#'
#' @name load_bgen
#'
#' @param in_bgen A string. BGEN file prefix (without .bgen extension). A matching .sample file must also be present.
#' @param verbose Logical. Be verbose (show individual steps),
#'        \code{default=FALSE}
#' @param very_verbose Logical. Be very verbose (show individual steps & show terminal output from plink etc),
#'        \code{default=FALSE}
#'
#' @examples
#'
#' liver_variants <- load_bgen(in_bgen="liver_cirrhosis.imputed.variants")
#'
#' @export
#'
load_bgen <- function(
	in_bgen,
	verbose=FALSE,
	very_verbose=FALSE
)  {

	if (class(in_bgen)[1] == "character")  {
		if (length(in_bgen) > 1)  cli::cli_abort("Input file path needs to be length 1")
		if (! file.exists(stringr::str_c(in_bgen, ".bgen")))    cli::cli_abort("Input BGEN file not found")
		if (! file.exists(stringr::str_c(in_bgen, ".sample")))  cli::cli_abort("Input sample file not found")
	}

	#
	# get plink2
	ukbrapR:::prep_tools(get_plink2=TRUE, verbose=verbose, very_verbose=very_verbose)

	# Step 1: BGEN -> PGEN
	# Converting to plink2's native PGEN format first ensures that multi-allelic loci stored as a
	# single tri/multi-allelic BGEN record are preserved intact before being split in step 2.
	# Unique CHR:POS:REF:ALT IDs are assigned here so that each allele can be distinguished.
	if (verbose) cli::cli_alert("Use plink2 to convert BGEN to PGEN")
	c1 <- paste0("~/_ukbrapr_tools/plink2 --bgen ", in_bgen, ".bgen ref-first --sample ", in_bgen,
	             ".sample --set-all-var-ids @:#:$r:$a --new-id-max-allele-len 200 missing --make-pgen --out _ukbrapr_tmp_pgen")
	if (very_verbose)  {
		system(c1)
	} else {
		system(stringr::str_c(c1, " >/dev/null"))
	}
	if (! file.exists("_ukbrapr_tmp_pgen.pgen"))  cli::cli_abort("plink2 failed to convert BGEN to PGEN. Try with `very_verbose=TRUE` to see terminal output.")

	# Step 2: PGEN -> RAW
	# --split-multiallelics turns each multi-allelic record into separate bi-allelic records so
	# that --export A produces one dosage column per alt allele (e.g., both T and G for a
	# tri-allelic SNP), matching the behaviour of the old BGEN->BED->RAW pipeline.
	if (verbose) cli::cli_alert("Use plink2 to export PGEN to RAW text file (splitting multi-allelic records)")
	c2 <- paste0("~/_ukbrapr_tools/plink2 --pfile _ukbrapr_tmp_pgen --split-multiallelics",
	             " --set-all-var-ids @:#:$r:$a --new-id-max-allele-len 200 missing --export A --out _ukbrapr_tmp_raw")
	if (very_verbose)  {
		system(c2)
	} else {
		system(stringr::str_c(c2, " >/dev/null"))
	}
	if (! file.exists("_ukbrapr_tmp_raw.raw"))  cli::cli_abort("plink2 failed to export RAW from PGEN. Try with `very_verbose=TRUE` to see terminal output.")

	# load genotype data to format
	if (verbose) cli::cli_alert("Read into memory and format")
	# PLINK RAW files are whitespace-delimited (tabs or spaces)
	geno_df <- utils::read.table("_ukbrapr_tmp_raw.raw", header=TRUE, sep="", check.names=FALSE, stringsAsFactors=FALSE)
	names(geno_df)[1] <- "eid"   # plink2 writes "#FID" as first column name
	geno_df <- geno_df |>
		dplyr::select(-dplyr::any_of(c("IID", "PAT", "MAT", "SEX", "PHENO")))

	# remove tmp files
	system("rm -f _ukbrapr_tmp_pgen.pgen _ukbrapr_tmp_pgen.pvar _ukbrapr_tmp_pgen.psam _ukbrapr_tmp_pgen.log _ukbrapr_tmp_raw.raw _ukbrapr_tmp_raw.log")

	# return
	return(geno_df)

}
