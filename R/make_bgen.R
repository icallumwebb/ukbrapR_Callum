#' Extract variants from imputed genotype file(s) into a single BGEN file
#'
#' @description For a given set of identifiers (RSID or "build37 position") extract from the UK Biobank imputed genotypes (BGEN files, field 22828) into a single BGEN file. Uses bgenix to subset each per-chromosome BGEN, then cat-bgen to concatenate them.
#'
#' @return A single BGEN file (and .sample file)
#'
#' @author Callum Webb, adapted from make_bed.R by Luke Pilling
#'
#' @name make_imputed_bgen
#'
#' @param in_file A data frame or file path. Contains at least two columns: (`rsID` and `CHR`) OR (`CHR` and `POS`). Other columns are ignored.
#' @param out_bgen A string. Output BGEN file prefix.
#' @param use_pos Logical. Use genomic position (CHR and POS) instead of RSID,
#'        \code{default=FALSE}
#' @param progress Logical. Show progress through each individual file,
#'        \code{default=TRUE}
#' @param verbose Logical. Be verbose (show individual steps),
#'        \code{default=FALSE}
#' @param very_verbose Logical. Be very verbose (show individual steps & show terminal output from bgenix/plink etc),
#'        \code{default=FALSE}
#'
#' @examples
#'
#' varlist <- data.frame(rsid=c("rs1800562","rs429358"), chr=c(6,19))
#' make_imputed_bgen(in_file=varlist, out_bgen="example_variants_by_rsid")
#'
#' @export
#'
make_imputed_bgen <- function(
    in_file,
    out_bgen,
    use_pos=FALSE,
    progress=TRUE,
    verbose=FALSE,
    very_verbose=FALSE
)  {

  # start
  start_time <- Sys.time()

  #
  #
  # check inputs
  if (very_verbose)  verbose <- TRUE
  if (verbose) cli::cli_alert("Checking inputs")

  # load user-provided varlist file
  varlist <- NULL

  if (class(in_file)[1] == "character")  {
    if (length(in_file)>1)  cli::cli_abort("Input file path needs to be length 1")
    if (! file.exists(in_file))  cli::cli_abort("Input file not found")
    varlist <- readr::read_tsv(in_file, progress=FALSE, show_col_types=FALSE)
  } else if (! any(class(in_file) %in% c("data.frame","tbl","tbl_df")))  {
    cli::cli_abort(c(
      "{.var in_file} must be a data.frame (or tibble), or a character string",
      "x" = "You've supplied a {.cls {class(in_file)}} vector."
    ))
  } else {
    varlist <- in_file
  }

  # check varlist formatting
  if (use_pos)  {
    varlist$rsid <- ""
    varlist <- ukbrapR:::prep_varlist(varlist, doing_pgs=FALSE, need_pos=TRUE, verbose=verbose)
  } else {
    varlist <- ukbrapR:::prep_varlist(varlist, doing_pgs=FALSE, need_pos=FALSE, verbose=verbose)
  }

  # check output format
  if (! class(out_bgen)=="character")  cli::cli_abort("Output file prefix needs to be a character string")
  if (length(out_bgen)>1)  cli::cli_abort("Output file prefix needs to be length 1")

  #
  # get bgenix and plink2
  ukbrapR:::prep_tools(get_plink=FALSE, get_plink2=TRUE, get_bgen=TRUE, verbose=verbose, very_verbose=very_verbose)

  #
  #
  # for each CHR
  chrs <- unique(varlist$chr)
  n_chrs <- length(chrs)

  # show progress
  cli::cli_alert("Extracting {nrow(varlist)} variant{?s} from {n_chrs} imputed BGEN file{?s}")

  bgen_files      <- character(0)
  sample_file_path <- NULL

  # loop over files...
  for (ii in 1:n_chrs)  {

    # this CHR
    chr      <- chrs[ii]
    chr_time <- Sys.time()

    # get variants list for this file
    if (use_pos)  {
      this_chr <- chr
      if (chr %in% c(1:9))  this_chr <- stringr::str_c("0", chr)   # imputed BGENs have 0 prefix for chrs <10
      varlist_sub <- varlist |> dplyr::filter(chr==!!chr) |> dplyr::mutate(variant_range=stringr::str_c(this_chr, ":", pos, "-", pos))
      readr::write_tsv(dplyr::select(varlist_sub, variant_range), "_ukbrapr_tmp_range.txt", col_names=FALSE, progress=FALSE)
    } else {
      varlist_sub <- varlist |> dplyr::filter(chr==!!chr)
      readr::write_tsv(dplyr::select(varlist_sub, rsid), "_ukbrapr_tmp_rsids.txt", col_names=FALSE, progress=FALSE)
    }

    # path to BGEN
    bgen_path <- stringr::str_c("/mnt/project/Bulk/Imputation/UKB\\ imputation\\ from\\ genotype/ukb22828_c", chr, "_b0_v3.bgen")

    # use bgenix to extract subset of BGEN into a per-chr tmp file
    if (verbose) cli::cli_alert(stringr::str_c("Using bgenix to extract the positions from chr", chr))
    tmp_bgen <- stringr::str_c("_ukbrapr_tmp_chr", chr, ".bgen")
    if (use_pos)  {
      c1 <- stringr::str_c("~/_ukbrapr_tools/bgenix -g ", bgen_path, " -incl-range _ukbrapr_tmp_range.txt > ", tmp_bgen)
    } else {
      c1 <- stringr::str_c("~/_ukbrapr_tools/bgenix -g ", bgen_path, " -incl-rsids _ukbrapr_tmp_rsids.txt > ", tmp_bgen)
    }
    if ( very_verbose)  system(c1)
    if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))

    # did it work?
    if (! file.exists(tmp_bgen))  cli::cli_abort("BGENIX failed to extract from the UKB imputed BGEN. Try with `very_verbose=TRUE` to see terminal output.")

    # does the BGEN actually contain variants? -- create index and check list file length
    c1 <- stringr::str_c("~/_ukbrapr_tools/bgenix -g ", tmp_bgen, " -index")
    if ( very_verbose)  system(c1)
    if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))
    c1 <- stringr::str_c("~/_ukbrapr_tools/bgenix -g ", tmp_bgen, " -list > _ukbrapr_tmp_list.txt")
    if ( very_verbose)  system(c1)
    if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))
    n_rows <- as.integer(system("wc -l < _ukbrapr_tmp_list.txt", intern=TRUE))
    system("rm _ukbrapr_tmp_list.txt")

    # if no variants in the BGEN (nrow of list file <=3) then skip this CHR
    if (n_rows > 3)  {
      bgen_files <- c(bgen_files, tmp_bgen)
      if (is.null(sample_file_path))  {
        sample_file_path <- stringr::str_replace_all(
          stringr::str_c("/mnt/project/Bulk/Imputation/UKB\\ imputation\\ from\\ genotype/ukb22828_c", chr, "_b0_v3.sample"),
          stringr::fixed("\\"), ""
        )
      }
    } else {
      cli::cli_warn(stringr::str_c("Variants on CHR ", chr, " are in the input varlist but are missing from imputed BGEN"))
      system(stringr::str_c("rm -f ", tmp_bgen, " ", tmp_bgen, ".bgi"))
    }

    system("rm -f _ukbrapr_tmp_range.txt _ukbrapr_tmp_rsids.txt")

    # give update
    if (progress)  cli::cli_alert_info(stringr::str_c("Extracted from imputed BGEN ", ii, " of ", n_chrs, " [", prettyunits::pretty_sec(as.numeric(difftime(Sys.time(), chr_time, units="secs"))), "]"))

  }

  # no variants found in any CHR?
  if (length(bgen_files) == 0)  cli::cli_abort("No variants were extracted from any imputed BGEN.")

  # concatenate all per-chr BGENs into a single BGEN
  if (length(bgen_files) == 1)  {
    if (verbose) cli::cli_alert("Renaming single-chr BGEN")
    system(stringr::str_c("mv ", bgen_files[1], " ", out_bgen, ".bgen"))
    system(stringr::str_c("mv ", bgen_files[1], ".bgi ", out_bgen, ".bgen.bgi"))
  } else {
    if (verbose) cli::cli_alert("Concatenating per-chr BGENs with cat-bgen")
    bgen_args <- stringr::str_c("-g ", bgen_files, collapse=" ")
    c1 <- stringr::str_c("~/_ukbrapr_tools/cat-bgen ", bgen_args, " -og ", out_bgen, ".bgen")
    if ( very_verbose)  system(c1)
    if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))
    for (f in bgen_files)  system(stringr::str_c("rm -f ", f, " ", f, ".bgi"))
  }

  # check output
  if (! file.exists(stringr::str_c(out_bgen, ".bgen")))  cli::cli_abort("Failed to create the BGEN. Try with `very_verbose=TRUE` to see terminal output.")

  # copy sample file alongside the BGEN
  out_sample <- stringr::str_c(out_bgen, ".sample")
  if (file.exists(out_sample))  unlink(out_sample, force=TRUE)
  ok_copy <- file.copy(sample_file_path, out_sample, overwrite=TRUE)
  if (!ok_copy)  cli::cli_abort("Failed to copy .sample file. Check write permissions in your working directory or choose a different output prefix.")

  # finished
  if (n_chrs>1)  {
    cli::cli_progress_done()
    options(cli.progress_show_after = 2)
  }

  cli::cli_alert_success(stringr::str_c("Imputed BGEN made in ", prettyunits::pretty_sec(as.numeric(difftime(Sys.time(), start_time, units="secs")))))

}


#' Extract variants from DRAGEN BGEN file(s) into a single BGEN file
#'
#' @description For a given set of genomic coordinates extract the UK Biobank WGS DRAGEN variant calls (from the BGEN format, field 24309) into a single BGEN file. Uses bgenix to subset each per-chromosome BGEN, then cat-bgen to concatenate them.
#'
#' This assumes your project has access to the WGS BGEN files released April 2025. If not, run `ukbrapR:::make_dragen_bed_from_pvcfs()`.
#'
#' @return A single BGEN file (and .sample file)
#'
#' @author Luke Pilling
#'
#' @name make_dragen_bgen
#'
#' @param in_file A data frame or file path. Contains at least two columns: `chr` and `pos` (in build 38). Other columns are ignored.
#' @param out_bgen A string. Output BGEN file prefix.
#' @param progress Logical. Show progress through each individual file,
#'        \code{default=TRUE}
#' @param verbose Logical. Be verbose (show individual steps),
#'        \code{default=FALSE}
#' @param very_verbose Logical. Be very verbose (show individual steps & show terminal output from bgenix/plink etc),
#'        \code{default=FALSE}
#'
#' @examples
#'
#' make_dragen_bgen(in_file=system.file("files", "pgs_liver_cirrhosis.txt", package="ukbrapR"), out_bgen="liver_cirrhosis.dragen.variants")
#'
#' @export
#'
make_dragen_bgen <- function(
	in_file,
	out_bgen,
	progress=TRUE,
	verbose=FALSE,
	very_verbose=FALSE
)  {

	# start
	start_time <- Sys.time()

	#
	#
	# check inputs
	if (very_verbose)  verbose <- TRUE
	if (verbose) cli::cli_alert("Checking inputs")

	# load user-provided varlist file
	varlist <- NULL

	if (class(in_file)[1] == "character")  {
		if (length(in_file)>1)  cli::cli_abort("Input file path needs to be length 1")
		if (! file.exists(in_file))  cli::cli_abort("Input file not found")
		varlist <- readr::read_tsv(in_file, progress=FALSE, show_col_types=FALSE)
	} else if (! any(class(in_file) %in% c("data.frame","tbl","tbl_df")))  {
		cli::cli_abort(c(
			"{.var in_file} must be a data.frame (or tibble), or a character string",
			"x" = "You've supplied a {.cls {class(in_file)}} vector."
		))
	} else {
		varlist <- in_file
	}

	varlist$rsid <- ""
	varlist <- ukbrapR:::prep_varlist(varlist, doing_pgs=FALSE, verbose=verbose)

	if (! class(out_bgen)=="character")  cli::cli_abort("Output file prefix needs to be a character string")
	if (length(out_bgen)>1)  cli::cli_abort("Output file prefix needs to be length 1")

	#
	# get bgenix and plink2
	ukbrapR:::prep_tools(get_plink=FALSE, get_plink2=TRUE, get_bgen=TRUE, verbose=verbose, very_verbose=very_verbose)

	#
	#
	# for each CHR
	chrs   <- unique(varlist$chr)
	n_chrs <- length(chrs)

	# show progress
	cli::cli_alert("Extracting {nrow(varlist)} variant{?s} from {n_chrs} DRAGEN BGEN file{?s}")

	bgen_files       <- character(0)
	sample_file_path <- NULL

	# loop over files...
	for (ii in 1:n_chrs)  {

		# this CHR
		chr      <- chrs[ii]
		chr_time <- Sys.time()

		# get variants list for this file
    varlist_sub <- varlist |> dplyr::filter(chr==!!chr) |> dplyr::mutate(variant_range=stringr::str_c(chr, ":", pos, "-", pos))
    readr::write_tsv(dplyr::select(varlist_sub, variant_range), "_ukbrapr_tmp_range.txt", col_names=FALSE, progress=FALSE)

		# path to BGEN
		bgen_path <- stringr::str_c("/mnt/project/Bulk/DRAGEN\\ WGS/DRAGEN\\ population\\ level\\ WGS\\ variants\\,\\ BGEN\\ format\\ \\[500k\\ release\\]/ukb24309_c", chr, "_b0_v1.bgen")

		# check it exists - exit if not
		if (! file.exists(stringr::str_replace_all(bgen_path, stringr::fixed("\\"), "")) )  {
			cli::cli_abort(c(
				stringr::str_c("DRAGEN BGEN file not found: ukb24309_c", chr, "_b0_v1.bgen"),
				"Has your Project been updated since April 2025? If not, you probably don't have the new BGENs.",
				"Consider using `ukbrapR:::make_dragen_bed_from_pvcfs()`"
			))
		}

		# use bgenix to extract subset of BGEN into a per-chr tmp file
		if (verbose) cli::cli_alert(stringr::str_c("Using bgenix to extract the positions from chr", chr))
		tmp_bgen <- stringr::str_c("_ukbrapr_tmp_chr", chr, ".bgen")
		c1 <- stringr::str_c("~/_ukbrapr_tools/bgenix -g ", bgen_path, " -incl-range _ukbrapr_tmp_range.txt > ", tmp_bgen)
		if ( very_verbose)  system(c1)
		if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))

		# did it work?
		if (! file.exists(tmp_bgen))  cli::cli_abort("BGENIX failed to extract from the UKB BGEN. Try with `very_verbose=TRUE` to see terminal output.")

		# does the BGEN actually contain variants? -- create index and check list file length
		c1 <- stringr::str_c("~/_ukbrapr_tools/bgenix -g ", tmp_bgen, " -index")
		if ( very_verbose)  system(c1)
		if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))
		c1 <- stringr::str_c("~/_ukbrapr_tools/bgenix -g ", tmp_bgen, " -list > _ukbrapr_tmp_list.txt")
		if ( very_verbose)  system(c1)
		if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))
		n_rows <- as.integer(system("wc -l < _ukbrapr_tmp_list.txt", intern=TRUE))
		system("rm _ukbrapr_tmp_list.txt")

		# if no variants in the BGEN (nrow of list file <=3) then skip this CHR
		if (n_rows > 3)  {
			bgen_files <- c(bgen_files, tmp_bgen)
			if (is.null(sample_file_path))  {
				sample_file_path <- stringr::str_replace_all(
					stringr::str_c("/mnt/project/Bulk/DRAGEN\\ WGS/DRAGEN\\ population\\ level\\ WGS\\ variants\\,\\ BGEN\\ format\\ \\[500k\\ release\\]/ukb24309_c", chr, "_b0_v1.sample"),
					stringr::fixed("\\"), ""
				)
			}
		} else {
			cli::cli_warn(stringr::str_c("Variants on CHR ", chr, " are in the input varlist but are missing from DRAGEN BGEN"))
			system(stringr::str_c("rm -f ", tmp_bgen, " ", tmp_bgen, ".bgi"))
		}

		system("rm -f _ukbrapr_tmp_range.txt")

		# give update
		if (progress)  cli::cli_alert_info(stringr::str_c("Extracted from DRAGEN BGEN chr", chr, " (", ii, " of ", n_chrs, ") [", prettyunits::pretty_sec(as.numeric(difftime(Sys.time(), chr_time, units="secs"))), "]"))

	}

	# no variants found in any CHR?
	if (length(bgen_files) == 0)  cli::cli_abort("No variants were extracted from any DRAGEN BGEN.")

	# concatenate all per-chr BGENs into a single BGEN
	if (length(bgen_files) == 1)  {
		if (verbose) cli::cli_alert("Renaming single-chr BGEN")
		system(stringr::str_c("mv ", bgen_files[1], " ", out_bgen, ".bgen"))
		system(stringr::str_c("mv ", bgen_files[1], ".bgi ", out_bgen, ".bgen.bgi"))
	} else {
		if (verbose) cli::cli_alert("Concatenating per-chr BGENs with cat-bgen")
		bgen_args <- stringr::str_c("-g ", bgen_files, collapse=" ")
		c1 <- stringr::str_c("~/_ukbrapr_tools/cat-bgen ", bgen_args, " -og ", out_bgen, ".bgen")
		if ( very_verbose)  system(c1)
		if (!very_verbose)  system(stringr::str_c(c1, " 2>/dev/null"))
		for (f in bgen_files)  system(stringr::str_c("rm -f ", f, " ", f, ".bgi"))
	}

	# check output
	if (! file.exists(stringr::str_c(out_bgen, ".bgen")))  cli::cli_abort("Failed to create the BGEN. Try with `very_verbose=TRUE` to see terminal output.")

	# copy sample file alongside the BGEN
  out_sample <- stringr::str_c(out_bgen, ".sample")
  if (file.exists(out_sample))  unlink(out_sample, force=TRUE)
  ok_copy <- file.copy(sample_file_path, out_sample, overwrite=TRUE)
  if (!ok_copy)  cli::cli_abort("Failed to copy .sample file. Check write permissions in your working directory or choose a different output prefix.")

	# finished
	if (n_chrs>1)  {
		cli::cli_progress_done()
		options(cli.progress_show_after = 2)
	}

	cli::cli_alert_success(stringr::str_c("DRAGEN BGEN made in ", prettyunits::pretty_sec(as.numeric(difftime(Sys.time(), start_time, units="secs")))))

}
