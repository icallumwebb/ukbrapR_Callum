#' Create a polygenic score
#'
#' @description Use a user-provided list of genetic variants with weights for a trait to create a polygenic score. Uses the imputed BGEN files (field 22828) or WGS DRAGEN BGEN files (field 24309).
#'
#' Uses plink2 to create the score (https://www.cog-genomics.org/plink/2.0/score). The returned score is the sum of the effect alleles weighted by the provided beta coefficients, divided by the number of non-missing alleles (i.e. the average score per allele).
#'
#' If selecting the DRAGEN data as the source, this assumes your project has access to the WGS BGEN files released April 2025. If not, run `ukbrapR:::make_dragen_bed_from_pvcfs()` to use [tabix] and [plink] to subset the [DRAGEN WGS pVCF files].
#'
#' @return A data frame
#'
#' @author Luke Pilling
#'
#' @name create_pgs
#'
#' @param in_file A data frame or file path. Must contain rsid, chr, pos, effect_allele, other_allele, beta. For imputed genos pos is build 37. For DRAGEN pos is build 38. Other columns are ignored.
#' @param out_file A string. Prefix for output files (optional)
#'        \code{default="tmp"}
#' @param pgs_name A string. Variable name for created PGS (optional)
#'        \code{default="pgs"}
#' @param source A string. Either "imputed" or "dragen" - indicating whether the variants should be from "UKB imputation from genotype" (field 22828) or "DRAGEN population level WGS variants, BGEN format [500k release]" (field 24309). Can instead be a path to a local BED file, if `is_bed=TRUE`.
#'        \code{default="imputed"}
#' @param use_imp_pos Logical. If source imputed, use position instead of rsID to extract variants?,
#'        \code{default=FALSE}
#' @param is_bed Logical. If you already have a BED file containing the required variants set this to TRUE and provide a path to the BED file in the `source` option,
#'        \code{default=FALSE}
#' @param overwrite Logical. Overwrite output genotype files? (If out_file is left as 'tmp' overwrite is set to TRUE). When `is_bed=FALSE`, this controls overwriting output BGEN files.
#'        \code{default=FALSE}
#' @param progress Logical. Show progress through each individual file,
#'        \code{default=FALSE}
#' @param verbose Logical. Be verbose (show individual steps),
#'        \code{default=FALSE}
#' @param very_verbose Logical. Be very verbose (show individual steps & show terminal output from Plink etc),
#'        \code{default=FALSE}
#'
#' @examples
#'
#' # example variant list and weights from GWAS of liver cirrhosis
#' #  - Innes 2020 Gastroenterology doi:10.1053/j.gastro.2020.06.014
#' #  - Position in build 38
#' varlist <- system.file("files", "pgs_liver_cirrhosis.txt", package="ukbrapR")
#'
#' # Create PGS from imputed data using RSID
#' liver_pgs <- create_pgs(in_file=varlist, out_file="liver_cirrhosis.imputed.pgs", pgs_name="liver_cirrhosis_imputed_pgs")
#'
#' # Create PGS from DRAGEN WGS data using CHR and POS
#' liver_pgs <- create_pgs(in_file=varlist, out_file="liver_cirrhosis.dragen.pgs", pgs_name="liver_cirrhosis_dragen_pgs", source="dragen")
#'
#' # For these allele weights, we has position in build 37 and will get imputed data using this not RSIDs
#' #  - Bladder Cancer GWAS, Graff 2021 (https://doi.org/10.1038/s41467-021-21288-z)
#' varlist2 <- readr::read_tsv("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000071/ScoringFiles/Harmonized/PGS000071_hmPOS_GRCh37.txt.gz", comment="#")
#' varlist2 <- varlist2 |> dplyr::rename(chr=chr_name, pos=chr_position)
#' bladder_cancer_pgs  <- create_pgs(in_file=varlist2, out_file="bladder_cancer.imputed.pgs", pgs_name="bladder_cancer_imputed_pgs", use_imp_pos=TRUE)
#'
#' @export
#'
create_pgs <- function(
    in_file,
    out_file="tmp",
    pgs_name="pgs",
    source="imputed",
    use_imp_pos=FALSE,
    is_bed=FALSE,
    overwrite=FALSE,
    progress=FALSE,
    verbose=FALSE,
    very_verbose=FALSE
)  {

	# start up messages
  pkg_version <- utils::packageVersion("ukbrapR")
  cli::cli_alert_info("{.pkg ukbrapR} v{pkg_version}")
  .ukbrapr_startup_notice()

  start_time <- Sys.time()

  #
  #
  # check inputs
  if (very_verbose)  verbose <- TRUE
  if (verbose) cli::cli_alert("Checking inputs")

  # local BED?
  geno_path <- NULL
  if (is_bed)  {

    geno_path <- source

    # have .bed suffix? If so, remove
    if (stringr::str_ends(geno_path, ".bed"))  geno_path <- stringr::str_sub(geno_path, end = -5)

    # check exists:
    if (! file.exists(stringr::str_c(geno_path, ".bed")))  cli::cli_abort("Local BED file not found")
    if (verbose) cli::cli_alert("Local BED file found")

    # is it DRAGEN format? see if .bim only contains "chr" IDs (no rsIDs)
    source <- "imputed"
    bim <- readr::read_tsv(stringr::str_c(geno_path, ".bim"), col_names=c("chr","id","null","pos","a1","a2"), progress=FALSE, show_col_types=FALSE)
    if ( all( stringr::str_detect(bim$id, "chr") ) )  source <- "dragen"

  }  else  {

    # imputed or dragen?
    if (! source %in% c("imputed","dragen"))  cli::cli_abort("{.var source} must be either \"imputed\" or \"dragen\"")

    geno_path <- out_file
  }

  # load user-provided varlist file (only first two TSV cols are used: must be chr, bp)
  varlist <- NULL

  # if it's a character string, assume user has provided a file path
  if (class(in_file)[1] == "character")  {

    if (length(in_file)>1)  cli::cli_abort("Input file path needs to be length 1")

    # does input file exist?
    if (! file.exists(in_file))  cli::cli_abort("Input file not found")
    varlist <- readr::read_tsv(in_file, progress=FALSE, show_col_types=FALSE)

  } else if (! any(class(in_file) %in% c("data.frame","tbl","tbl_df")))  {

    cli::cli_abort(c(
      "{.var in_file} must be a data.frame (or tibble), or a character string",
      "x" = "You've supplied a {.cls {class(in_file)}}."
    ))

  } else {  # user has passed a data frame
    varlist <- in_file
  }

  # check varlist formatting and save
  readr::write_tsv(varlist, stringr::str_c(out_file, ".input.txt"))
  varlist <- ukbrapR:::prep_varlist(varlist, doing_pgs=TRUE, verbose=verbose)
  varlist <- varlist |>
    dplyr::mutate(
      effect_allele=stringr::str_to_upper(effect_allele),
      other_allele=stringr::str_to_upper(other_allele)
    )
  out_file_varlist <- stringr::str_c(out_file, ".varlist.txt")

  # check output format
  if (! class(out_file)=="character")  cli::cli_abort("Output file prefix needs to be a character string")
  if (length(out_file)>1)  cli::cli_abort("Output file prefix needs to be length 1")
  if (out_file=="tmp")  overwrite <- TRUE
  if (!is_bed & file.exists(paste0(out_file,".bgen")) & !overwrite)  cli::cli_abort("Output BGEN already exists. To overwrite, set option `overwrite=TRUE`")

  #
  #
  # prepare genotype files -- if local BED provided then just check plink2 is available
    if (!is_bed & source == "imputed")  ukbrapR::make_imputed_bgen(in_file=varlist, out_bgen=geno_path, use_pos=use_imp_pos, progress=progress, verbose=verbose, very_verbose=very_verbose)
    if (!is_bed & source == "dragen")   ukbrapR::make_dragen_bgen(in_file=varlist, out_bgen=geno_path, progress=progress, verbose=verbose, very_verbose=very_verbose)
  if (is_bed)  ukbrapR:::prep_tools(get_plink=FALSE, get_plink2=TRUE, get_bgen=FALSE, verbose=verbose, very_verbose=very_verbose)


  # did it work?
  if ( is_bed & ! file.exists(stringr::str_c(geno_path, ".bed")))   cli::cli_abort("Local BED file not found. Try with `very_verbose=TRUE` to see terminal output.")
  if (!is_bed & ! file.exists(stringr::str_c(geno_path, ".bgen")))  cli::cli_abort("Failed to make the BGEN. Try with `very_verbose=TRUE` to see terminal output.")

  #
  #
  # create PGS

    # Replace user-provided `rsid` with internal IDs when needed.
    # - DRAGEN BGENs generally require CHR:POS:REF:ALT-like IDs.
    # - Imputed BGEN with use_imp_pos=TRUE also benefits from ID remapping, since many
    #   extracted variants do not retain user RSIDs as plink2 scoring IDs.

    varlist$rsid_old <- varlist$rsid
    need_id_remap <- identical(source, "dragen") || (!is_bed && identical(source, "imputed") && isTRUE(use_imp_pos))

    if (need_id_remap && is_bed)  {

      # BED path: read the BIM file and match by CHR, POS, and alleles
      varinfo <- readr::read_tsv(stringr::str_c(geno_path, ".bim"), col_names=c("chr","id","null","pos","a1","a2"), progress=FALSE, show_col_types=FALSE)

    } else if (need_id_remap) {

      # BGEN path: use plink2 to write a PGEN (which includes a .pvar variant info file)
      # and match by CHR, POS, and alleles. We use --make-pgen for compatibility with
      # plink2 builds where --make-pvar is not available.
      if (verbose) cli::cli_alert("Getting variant info from BGEN")
      c1 <- paste0("~/_ukbrapr_tools/plink2 --bgen ", geno_path, ".bgen ref-first --sample ", geno_path, ".sample --make-pgen --out _ukbrapr_tmp_pvar")
      if (very_verbose)  {
        system(c1)
      } else {
        system(stringr::str_c(c1, " >/dev/null"))
      }
      if (! file.exists("_ukbrapr_tmp_pvar.pvar"))  cli::cli_abort("Plink2 failed to write variant info (.pvar) from BGEN. Try with `very_verbose=TRUE` to see terminal output.")
      # PVAR header: ##fileformat=PARvX.X then #CHROM POS ID REF ALT
      # read_tsv with comment="##" skips the ##fileformat line; #CHROM becomes the first column name
      varinfo_raw <- readr::read_tsv("_ukbrapr_tmp_pvar.pvar", comment="##", progress=FALSE, show_col_types=FALSE)
      varinfo <- varinfo_raw |>
        dplyr::rename(id=ID, pos=POS, a1=REF, a2=ALT) |>
        dplyr::mutate(
          chr=as.integer(stringr::str_remove(as.character(.data[[names(varinfo_raw)[1]]]), "^chr")),
          a1=stringr::str_to_upper(a1),
          a2=stringr::str_to_upper(a2)
        ) |>
        tidyr::separate_rows(a2, sep=",") |>
        dplyr::mutate(a2=stringr::str_trim(a2)) |>
        dplyr::select(chr, id, pos, a1, a2)
      system("rm -f _ukbrapr_tmp_pvar*")

    }

    if (need_id_remap && verbose)  {
      n_mapped <- sum(varlist$rsid != varlist$rsid_old, na.rm=TRUE)
      cli::cli_alert_info("Mapped {n_mapped} of {nrow(varlist)} score IDs to genotype IDs")
    }

    # create ID for each row of the varlist - make sure alleles match
    if (need_id_remap) for (ii in 1:nrow(varlist))  {

      # keep rows where CHR and POS match
      r <- varinfo[ varinfo$chr==varlist$chr[ii] & varinfo$pos==varlist$pos[ii] , ]

      # any matched?
      if (nrow(r)>0)  {

        # keep rows where both alleles are genotyped
        r <- r[ r$a1 %in% c(varlist$effect_allele[ii],varlist$other_allele[ii]) & r$a2 %in% c(varlist$effect_allele[ii],varlist$other_allele[ii]) , ]

        # any matched?
        if (nrow(r)>0)  {

          # update rsid to match plink2's internal variant ID
          # if multiple (shouldn't be!) use first
          varlist$rsid[ii] <- r$id[1]

        }

      }

    }

  # save the varlist for plink
  readr::write_tsv(varlist, out_file_varlist, progress=FALSE)

  # Plink
  if (verbose) cli::cli_alert("Make PGS")
    if (is_bed)  {
      c1 <- paste0("~/_ukbrapr_tools/plink2 --bfile ", geno_path, " --score ", out_file_varlist, " 1 4 6 header cols=+scoresums,+scoreavgs --out ", out_file)
    } else {
      c1 <- paste0("~/_ukbrapr_tools/plink2 --bgen ", geno_path, ".bgen ref-first --sample ", geno_path, ".sample --score ", out_file_varlist, " 1 4 6 header cols=+scoresums,+scoreavgs --out ", out_file)
    }
    if (very_verbose)  {
      system(c1)
    } else {
      system(stringr::str_c(c1, " >/dev/null"))
    }

  # did it work?
  if (! file.exists(stringr::str_c(out_file, ".sscore")))  cli::cli_abort("Plink failed to make the allele score. Try with `very_verbose=TRUE` to see terminal output.")

  # just extract EID and SCORE to a .tsv file -- remove participants with invalid EIDs < 0
  system(stringr::str_c("echo \"eid\t", pgs_name, "\" > ", out_file, ".tsv"))
  system(stringr::str_c("awk 'NR > 1 && $1 > 0 { print $1\"\t\"$5 }' ", out_file, ".sscore >> ", out_file, ".tsv"))

  # load
  pgs <- readr::read_tsv(stringr::str_c(out_file, ".tsv"), progress=FALSE, show_col_types=FALSE)

  #
  #
  # finished
  cli::cli_alert_success(stringr::str_c("PGS created! See file {.file ", out_file, ".tsv}"))
  if (verbose) cli::cli_alert_info(c("Time taken: ", "{prettyunits::pretty_sec(as.numeric(difftime(Sys.time(), start_time, units=\"secs\")))}."))

  return(pgs)

}
