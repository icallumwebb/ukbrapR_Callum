
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

	# use plink2 to export BGEN to RAW text file
	# For multi-allelic records, plink2 --export A produces one dosage column per alt allele.
	# --set-all-var-ids @_$a names each variant as variant_ALT_ALLELE (e.g., rs72725854_T, rs72725854_G)
	# so column names are unique and correspond to the alt allele being exported.
	if (verbose) cli::cli_alert("Use plink2 to export BGEN to RAW text file")
	c1 <- paste0("~/_ukbrapr_tools/plink2 --bgen ", in_bgen, ".bgen ref-first --sample ", in_bgen,
	             ".sample --set-all-var-ids @_$a --export A --out _ukbrapr_tmp")
	if (very_verbose)  {
		system(c1)
	} else {
		system(stringr::str_c(c1, " >/dev/null"))
	}
	if (! file.exists("_ukbrapr_tmp.raw"))  cli::cli_abort("plink2 failed to export RAW from BGEN. Try with `very_verbose=TRUE` to see terminal output.")

	# load genotype data to format
	if (verbose) cli::cli_alert("Read into memory and format")
	# PLINK RAW files are whitespace-delimited (tabs or spaces)
	geno_df <- utils::read.table("_ukbrapr_tmp.raw", header=TRUE, sep="", check.names=FALSE, stringsAsFactors=FALSE)
	names(geno_df)[1] <- "eid"   # plink2 writes "#FID" as first column name
	
	# For multi-allelic variants, plink2 may generate duplicate column names (e.g., both alt alleles
	# named "rs72725854_A"). Repair these before dplyr::select() which requires unique names.
	# Only apply to genotype columns (skip metadata FID/IID/PAT/MAT/SEX/PHENO columns).
	metadata_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO", "PHENOTYPE", "eid")
	geno_col_idx <- which(!(names(geno_df) %in% metadata_cols))
	if (length(geno_col_idx) > 0 && anyDuplicated(names(geno_df)[geno_col_idx]) > 0)  {
		new_names <- names(geno_df)
		new_names[geno_col_idx] <- make.unique(names(geno_df)[geno_col_idx], sep="_dup")
		names(geno_df) <- new_names
		if (verbose) cli::cli_alert_info("Repaired duplicate genotype column names")
	}
	
	geno_df <- geno_df |>
		dplyr::select(-dplyr::any_of(c("IID", "PAT", "MAT", "SEX", "PHENO")))

	# remove tmp files
	system("rm -f _ukbrapr_tmp.raw _ukbrapr_tmp.log")

	# return
	return(geno_df)

}
