# ------------------------ Plasmid analysis ------------------------

## This track of the ariba analysis script analyses plasmid reports
## from ARIBA and gives summary statistic reports on the inc-types
## reported

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
report_loc <- args[1]
output_loc <- args[2]
ending <- as.character(args[3])

# ------------------------ Load libraries -------------------------

packages <-
  c(
    "dplyr",
    "tidyr",
    "purrr",
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

# Create output folder
dir.create(paste0(output_loc, "/plasmid/"), showWarnings = FALSE)
plasmid_output <- paste0(output_loc, "/plasmid/")

plasmid_data <- get_data(report_loc,
                         ending,
                         convert = TRUE) %>%
  fix_gene_names(ending)

plasmid_flags <- check_flags(plasmid_data)
plasmid_table <- create_table(plasmid_data)
plasmid_report <- create_report(plasmid_table)
plasmid_quant <- calc_stats(plasmid_table)

# Write to output folder
write.table(plasmid_report,
            paste0(plasmid_output, "plasmid_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(plasmid_quant,
            paste0(plasmid_output, "plasmid_stats.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(plasmid_flags,
            paste0(plasmid_output, "plasmid_flags.tsv"),
            sep = "\t",
            row.names = FALSE)