# ------------------------ Plasmid analysis ------------------------

## This track of the ariba analysis script analyses plasmid reports
## from ARIBA and gives summary statistic reports for selected genes 
## of interest.

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

# Confidence interval function
get_binCI <- function(x, n) as.numeric(
  setNames(binom.test(x,n)$conf.int*100,c("lwr", "upr"))
)

# Corrects the gene names found in the "cluster" column
fix_plasmid_names <- function(df) {
  genes <- unique(df$ref_name)
  new_names <- gsub("^(.*?)_.+", "\\1", genes)
  
  gene_names <- c()
  for (i in new_names) {
    p <- paste(tolower(substring(i, 1,3)),
               substring(i, 4),
               sep = "", 
               collapse = " ")
    gene_names <- c(gene_names,p)
  }
  df2 <- data.frame(genes,gene_names) %>%
    mutate(genes = as.character(genes)) %>%
    rename(ref_name = genes)
  
  df <- df %>%
    left_join(df2, by = "ref_name") %>%
    mutate(gene_names = as.character(gene_names),
           ref = gsub("(.*?)/report.tsv", "\\1", ref))
  
  return(df)
}

# Allowed flags from the ARIBA report
flag_selection <- c("19","27","147","155","403",
                    "411","915","923","787","795",
                    "531","539","659","667","787","795")  

# Function that returns info on flag selection
check_flags <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(flag_result = flag %in% flag_selection,
           flag_result = as.integer(flag_result)) %>%
    rename("gene" = gene_names)
  return(df)
}

# Function that handles plasmid data and returns a data frame with 
# presence/absence for plasmid inc types                        
create_plasmid_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, flag)) %>%
    gather(gene, result,-ref) %>%
    mutate(result_total = ifelse(result == "", 0, 1),
           result_total = as.character(result_total),
           type = "gene") %>%
    select(-result)
  return(df)
}

# generate plasmid report
create_plasmid_report <- function(df) {
  df <- df %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, type))
  return(df)
}

calc_stats <- function(df) {
  df <- df %>%
    group_by(gene, result_total) %>%
    count() %>%
    ungroup() %>%
    mutate(result_total = if_else(result_total == 1, 
                                  "Present",
                                  "Absent")) %>%
    spread(result_total, n, fill = 0) %>%
    rowwise() %>%
    mutate(
      Total = Present + Absent,
      Percent = round(Present / Total * 100, 1),
      lwr = round(get_binCI(Present, Total)[1], 1),
      upr = round(get_binCI(Present, Total)[2], 1)
    )
  return(df)
}

# -------------------------- Analysis ----------------------------

# Create output folder
dir.create(paste0(output_loc, "/plasmid/"), showWarnings = FALSE)
plasmid_output <- paste0(output_loc, "/plasmid/")

plasmid_data <- get_data(report_loc,
                         "^report.tsv") %>%
  fix_plasmid_names()

plasmid_flags <- check_flags(plasmid_data)
plasmid_table <- create_plasmid_table(plasmid_data)
plasmid_report <- create_plasmid_report(plasmid_table)

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