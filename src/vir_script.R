# ----------------------- Virulence analysis ----------------------

## This track of the ariba analysis script analyses virulence
## reports from ARIBA and gives summary statistic reports for 
## selected genes of interest.

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

# Confidence interval function
get_binCI <- function(x, n) as.numeric(
  setNames(binom.test(x,n)$conf.int*100,c("lwr", "upr"))
  )

# Identifies filenames in input folder
file_names_vir <- function(filepath) {
  files <- list.files(path = filepath,
                      pattern = "vir_report.tsv")
  return(files)
}

# Import virulence data from report.tsv
get_vir_data <- function(filepath) {
  files <- file_names_vir(filepath)
  
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

# Corrects the gene names found in the "cluster" column
fix_vir_names <- function(df) {
  genes <- unique(df$ref_name)
  new_names <- gsub("^(.+_[0-9]+)_.+", "\\1", genes)
  
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
           ref = gsub("(.*?)_vir_report.tsv", "\\1", ref))
  
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
           flag_result = as.integer(flag_result),
           ref = gsub("(.*?)\\_.+", "\\1", ref),
           ref = sub("^\\d*-?(\\d{4}-.*)", "\\1", ref),
           ref = sub("^(\\d{4}-\\d{2}-\\d*)-1", "\\1", ref)) %>%
    rename("gene" = gene_names)
  return(df)
}

# Function that handles vir data and returns a data frame with 
# presence/absence for virulence genes                           
create_vir_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_vir_report.tsv$",
                      "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, result,-ref) %>%
    mutate(result_total = ifelse(result == "", 0, 1),
           result_total = as.character(result_total)) %>%
    select(-result)
  return(df)
}

# calculates percentages for presence/absence
calc_stats <- function(df) {
  df <- df %>%
    group_by(gene, result_total) %>%
    count() %>%
    ungroup() %>%
    mutate(result_total = if_else(result_total == 1, 
                                  "Present",
                                  "Absent")) %>%
    spread(result_total, n, fill = 0) %>%
    mutate(Absent = if ("Absent" %in% names(.)) {
      return(Absent)
    } else {
      return(0) # create absent column if not there
    }) %>%
    rowwise() %>%
    mutate(
      Total = Present + Absent,
      Percent = round(Present / Total * 100, 1),
      lwr = round(get_binCI(Present, Total)[1], 1),
      upr = round(get_binCI(Present, Total)[2], 1)
    )
  return(df)
}

# generate virulence gene report
create_vir_report <- function(df) {
  df <- df %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-id)
  return(df)
}

create_summary_report <- function(df) {
  df <- df %>%
    mutate(gene = sub("(.*?)_.+", "\\1", gene)) %>%
    group_by(ref, gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-id) %>%
    mutate_at(vars(-ref),
              funs(ifelse(. == "0, 1", "1", "0")))
  return(df)
}

# --------------------------- Analysis ---------------------------

# Create output directory
dir.create(paste0(output_loc, "/vir/"), showWarnings = FALSE)
vir_output <- paste0(output_loc, "/vir/")

# Import data
vir_data <- get_vir_data(report_loc)

# Clean data
clean_vir_data <- fix_vir_names(vir_data)

# Check flags
vir_flags <- check_flags(clean_vir_data)

# Wrangle
vir_table <- create_vir_table(clean_vir_data)
vir_report <- create_vir_report(vir_table)
vir_summary <- create_summary_report(vir_table)
vir_quant <- calc_stats(vir_table)

# Write to output folder
write.table(vir_report,
            paste0(vir_output, "virulence_detailed_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(vir_summary,
            paste0(vir_output, "virulence_summary_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(vir_quant,
            paste0(vir_output, "virulence_stats.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(vir_flags,
            paste0(vir_output, "virulence_flags.tsv"),
            sep = "\t",
            row.names = FALSE)