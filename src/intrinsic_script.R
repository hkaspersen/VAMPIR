# ------------- AMR gene analysis - Intrinsic genes ---------------

## This track of the ARIBA analysis script analyses AMR gene reports
## from ARIBA and generates result files based on user input of
## selected genes of interest.

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
in_report_loc <- args[1]
output_loc <- args[2]
genes <- args[3]
ending <- as.character(args[4])
gyr_par_fix <- args[5]

# adjust parameters for filtering
if (grepl("all", genes, ignore.case = TRUE) == TRUE) {
    genes <- "ALL"
  } else {
    genes <- unlist(strsplit(genes, ",", fixed = TRUE))
  }

# ------------------------ Load libraries -------------------------

packages <-
  c(
    "dplyr",
    "tidyr",
    "purrr",
    "stringr",
    "impoRt",
    "vampfunc",
    "funtools"
  )

suppressPackageStartupMessages(
  invisible(lapply(packages, function(x)
  library(
    x,
    character.only = T,
    quietly = T,
    warn.conflicts = FALSE
  ))))

# -------------------------- Analysis ----------------------------

# Create output directory
dir.create(paste0(output_loc, "/amr_in/"), showWarnings = FALSE)
amr_output <- paste0(output_loc, "/amr_in/")

## Intrinsic genes

in_data <- get_data(in_report_loc,
                    ending,
                    convert = TRUE) %>%
  fix_gene_names(ending, db = "res")

in_flags <- check_flags(in_data)

write.table(in_flags,
            paste0(amr_output, "intrinsic_flag_report.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

if (all(in_flags$flag_result == 0) == TRUE) {
  print("No flags accepted, please check the flag report")
  stop()
}

in_table <- create_table(in_data, acquired = FALSE)

if (exists("gyr_par_fix") == TRUE) {
  in_table <- fix_gyr_par_results(in_table)
} else {
  in_table <- in_table
}

if ("ALL" %in% genes) {
  in_table_filtered <- in_table
} else {
  in_table_filtered <- filter_table(in_table, genes)
  in_flags <- filter_table(in_flags, genes)
}

in_report <- create_report(in_table_filtered, mut = FALSE)
in_mut_report <- create_report(in_table_filtered, mut = TRUE)

in_stats <- calc_stats(in_table_filtered)

## Write results to file

write.table(in_report,
            paste0(amr_output, "intrinsic_gene_report.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(in_mut_report,
            paste0(amr_output, "intrinsic_mut_report.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(in_stats,
            paste0(amr_output, "intrinsic_gene_stats.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)