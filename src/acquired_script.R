# --------------- AMR gene analysis - Acquired genes --------------

## This track of the ARIBA analysis script analyses AMR gene 
## reports from ARIBA and generates result files based on user
## input of selected genes of interest.

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
ac_report_loc <- args[1]
output_loc <- args[2]
genes <- args[3]
ending <- as.character(args[4])

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
    "stringr",
    "knitr",
    "impoRt",
    "vampfunc",
    "funtools"
  )

invisible(lapply(packages, function(x)
  library(
    x,
    character.only = T,
    quietly = T,
    warn.conflicts = FALSE
  )))

# -------------------------- Analysis ----------------------------

# Create output directory
dir.create(paste0(output_loc, "/amr_ac/"), showWarnings = FALSE)
amr_output <- paste0(output_loc, "/amr_ac/")

ac_data <- get_data(ac_report_loc,
                    ending,
                    convert = TRUE) %>%
  fix_gene_names(ending, db = "res")

ac_flags <- check_flags(ac_data)

write.table(
  ac_flags,
  paste0(amr_output, "acquired_flag_report.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (all(ac_flags$flag_result == 0) == TRUE) {
  print("No flags accepted, please check the flag report")
  stop()
}

ac_table <- create_table(ac_data)

if ("ALL" %in% genes) {
  ac_table_filtered <- ac_table
} else {
  ac_table_filtered <- filter_table(ac_table, genes)
  ac_flags <- filter_table(ac_flags, genes)
}

ac_report <- create_report(ac_table_filtered)
ac_stats <- calc_stats(ac_table_filtered)

## Write results to file

write.table(
  ac_report,
  paste0(amr_output, "acquired_gene_report.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  ac_stats,
  paste0(amr_output, "acquired_gene_stats.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)