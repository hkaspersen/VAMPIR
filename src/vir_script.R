# ----------------------- Virulence analysis ----------------------

## This track of the ariba analysis script analyses virulence
## reports from ARIBA and gives summary statistic reports for 
## selected genes of interest.

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
report_loc <- args[1]
vir_database <- args[2]
vir_genes <- args[3]
output_loc <- args[4]
ending <- as.character(args[5])

# adjust parameters for filtering
if (grepl("all", vir_genes, ignore.case = TRUE) == TRUE) {
  vir_genes <- "ALL"
} else {
  vir_genes <- unlist(strsplit(vir_genes, ",", fixed = TRUE))
}

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

# --------------------------- Analysis ---------------------------

# Create output directory
dir.create(paste0(output_loc, "/vir/"), showWarnings = FALSE)
vir_output <- paste0(output_loc, "/vir/")

# import data
if (vir_database == "virfinder") {
  vir_data <- get_data(report_loc,
                       ending,
                       convert = TRUE) %>%
    fix_gene_names(ending, db = "virfinder") 
}

if (vir_database == "vfdb") {
  vir_data <- get_data(report_loc,
                       ending,
                       convert = TRUE) %>%
    fix_gene_names(ending, db = "vfdb")
}

vir_flags <- check_flags(vir_data)

write.table(vir_flags,
            paste0(vir_output, "virulence_flags.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

if (all(vir_flags$flag_result == 0) == TRUE) {
  print("No flags accepted, please check the flag report")
  stop()
}

vir_table <- create_table(vir_data, acquired = TRUE)
  
if ("ALL" %in% vir_genes) {
  vir_table_filtered <- vir_table
} else {
  vir_table_filtered <- filter_table(vir_table, vir_genes)
}
  
vir_report <- create_report(vir_table_filtered)
vir_quant <- calc_stats(vir_table_filtered)
  
write.table(vir_report,
            paste0(vir_output, "virulence_report_detailed.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(vir_quant,
            paste0(vir_output, "virulence_stats.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
