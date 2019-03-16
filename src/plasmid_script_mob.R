# ------------------------ Plasmid analysis ------------------------

## This track analyses plasmid reports from Mob-Suite

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
report_loc <- args[1]
output_loc <- args[2]

# ------------------------ Load libraries -------------------------

suppressPackageStartupMessages(if (!require("pacman")) 
  install.packages("pacman"))
suppressPackageStartupMessages(
  pacman::p_load(
    ggplot2, 
    dplyr,
    tidyr,
    purrr
  )
)

# -------------------------- Functions ----------------------------

# Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]),
                                collapse = ", ")

# Identifies filenames in input folder
file_names_plasmid <- function(filepath) {
  files <- list.files(path = filepath,
                      pattern = "mobtyper_aggregate_report.txt",
                      recursive = T)
  return(files)
}

# Import plasmid data from report.tsv
get_plasmid_data <- function(filepath) {
  files <- file_names_plasmid(filepath)
  
  data_list <- lapply(files,
                      FUN = function(file) {
                        read.delim(
                          paste0(filepath, "/", file),
                          stringsAsFactors = F,
                          header = TRUE,
                          sep = "\t"
                        )
                      })
  
  names(data_list) <- files
  data <- bind_rows(lapply(
    data_list, function(x) map(x, as.character)
  ), .id = "ref")
  return(data)
}

# -------------------------- Analysis ----------------------------

# Create output folder
dir.create(paste0(output_loc, "/plasmid/"), showWarnings = FALSE)
plasmid_output <- paste0(output_loc, "/plasmid/")

plasmid_data <- get_plasmid_data(report_loc)

clean_plasmid_data <- plasmid_data %>%
  mutate(ref = sub("/mobtyper_aggregate_report.txt", "", ref))

# Write to output folder
write.table(clean_plasmid_data,
            paste0(plasmid_output, "plasmid_report.tsv"),
            sep = "\t",
            row.names = FALSE)