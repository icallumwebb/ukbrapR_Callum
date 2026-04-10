
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
	if (verbose) cli::cli_alert("Use plink2 to export BGEN to RAW text file")
	c1 <- paste0("~/_ukbrapr_tools/plink2 --bgen ", in_bgen, ".bgen ref-first --sample ", in_bgen, ".sample --export A --out _ukbrapr_tmp")
	if (very_verbose)  {
		system(c1)
	} else {
		system(stringr::str_c(c1, " >/dev/null"))
	}

	# load genotype data to format
	if (verbose) cli::cli_alert("Read into memory and format")
	# PLINK RAW files are whitespace-delimited (tabs or spaces), so parse with read_table2
	geno_df <- readr::read_table2("_ukbrapr_tmp.raw", progress=FALSE, show_col_types=FALSE)
	names(geno_df)[1] <- "eid"   # plink2 writes "#FID" as first column name
	geno_df <- geno_df |>
		dplyr::select(-dplyr::any_of(c("IID", "PAT", "MAT", "SEX", "PHENO")))

	# remove tmp files
	system("rm _ukbrapr_tmp*")

	# return
	return(geno_df)

}
