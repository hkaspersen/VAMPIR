# ------------------------ Plasmid analysis ------------------------

## This track analyses plasmid reports from Mob-Suite

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
report_loc <- args[1]
output_loc <- args[2]

# ------------------------ Load libraries -------------------------

packages <-
  c(
    "ggplot2",
    "dplyr",
    "tidyr",
    "purrr",
    "impoRt"
  )

invisible(lapply(packages, function(x)
  library(
    x,
    character.only = T,
    quietly = T,
    warn.conflicts = FALSE
  )))

# -------------------------- Functions ----------------------------

# Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]),
                                collapse = ", ")

# -------------------------- Analysis ----------------------------

# Create output folder
dir.create(paste0(output_loc, "/plasmid/"), showWarnings = FALSE)
plasmid_output <- paste0(output_loc, "/plasmid/")

plasmid_data <- get_data(report_loc,
                         "mobtyper_aggregate_report.txt",
                         convert = TRUE)

clean_plasmid_data <- plasmid_data %>%
  mutate(ref = sub("/mobtyper_aggregate_report.txt", "", ref))

# Write to output folder
write.table(clean_plasmid_data,
            paste0(plasmid_output, "plasmid_report.tsv"),
            sep = "\t",
            row.names = FALSE)