# ------------- AMR gene analysis - Intrinsic genes ---------------

## This track of the ARIBA analysis script analyses AMR gene reports
## from ARIBA and generates result files based on user input of
## selected genes of interest.

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
in_report_loc <- args[1]
output_loc <- args[2]
in_genes <- args[3]

# adjust parameters for filtering
if (grepl("all", in_genes, ignore.case = TRUE) == TRUE) {
    in_genes <- "ALL"
  } else {
    in_genes <- unlist(strsplit(in_genes, ",", fixed = TRUE))
}

# ------------------------ Load libraries -------------------------

suppressPackageStartupMessages(if (!require("pacman")) 
  install.packages("pacman"))
suppressPackageStartupMessages(
  pacman::p_load(
    ggplot2,
    dplyr,
    tidyr,
    gridExtra,
    grid,
    forcats,
    purrr,
    stringr,
    kableExtra,
    knitr,
    IRdisplay,
    reprex,
    svglite
    )
)

# -------------------------- Functions ----------------------------

# Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

# Confidence interval function
get_binCI <- function(x, n) as.numeric(setNames(binom.test(x,n)$conf.int*100,
                                                c("lwr", "upr")))

# Identifies filenames in input folder
file_names <- function(filepath) {
  files <- list.files(path = filepath, pattern = "amr_report.tsv")
  return(files)
}

# Import ariba data from report.tsv from chosen database used in ariba
get_ariba_data <- function(filepath) {
  files <- file_names(filepath)
  
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
  data <- bind_rows(lapply(data_list, function(x) map(x, as.character)),
                    .id = "ref")
  return(data)
}

# Corrects the gene names found in the "cluster" column
fix_gene_names <- function(df) {
  genes <- unique(df$ref_name)
  new_names <- gsub("^(.*?)\\..*", "\\1", genes)
  new_names <- gsub("_", "", new_names, fixed = T)
  new_names <- gsub("-", "", new_names, fixed = T)
  
  gene_names <- c()
  for (i in new_names) {
    p <- paste(tolower(substring(i, 1,3)), substring(i, 4), sep = "", 
               collapse = " ")
    gene_names <- c(gene_names,p)
  }
  df2 <- data.frame(genes,gene_names) %>%
    mutate(genes = as.character(genes)) %>%
    rename(ref_name = genes)
  
  df <- df %>%
    left_join(df2, by = "ref_name") %>%
    mutate(gene_names = as.character(gene_names),
           ref = gsub("(.*?)_amr_report.tsv", "\\1", ref))
  
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

# Corrects or mutations in QRDR for gyrA, gyrB, parC and parE genes.
# The function returns columns of 1/0 values for whether or not the 
# mutations reported by ARIBA is within the QRDR in the gene:
# gyrA: AA 67 - 106
# gyrB: AA 333 - 481
# ParC: AA 51 - 170
# parE: AA 366 - 523
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142146/pdf/cjvr68pg229.pdf
fix_gyr_par_results <- function(df) {
  df <- df %>%
    mutate(gyrA_result = mut %>% # control for mutation within QRDR for gyrA
             str_extract_all("\\d+") %>% # from reprex package
             map(as.integer) %>% # converts all to integer
             map_lgl(~ any(.x >= 67L & .x <= 106L)), # returns TRUE/FALSE whether value is within range or not
           gyrA_result = if_else(gene != "gyrA", NA, gyrA_result), # filters out results for all other genes
           gyrA_result = as.integer(gyrA_result),
           gyrB_result = mut %>% # control for mutation within QRDR for gyrB
             str_extract_all("\\d+") %>%
             map(as.integer) %>%
             map_lgl(~ any(.x >= 333L & .x <= 481L)),
           gyrB_result = if_else(gene != "gyrB", NA, gyrB_result),
           gyrB_result = as.integer(gyrB_result),
           parC_result = mut %>% # control for mutation within QRDR for parC
             str_extract_all("\\d+") %>% 
             map(as.integer) %>% 
             map_lgl(~ any(.x >= 51L & .x <= 170L)),
           parC_result = if_else(gene != "parC", NA, parC_result),
           parC_result = as.integer(parC_result),
           parE_result = mut %>% # control for mutation within QRDR for parE
             str_extract_all("\\d+") %>% 
             map(as.integer) %>% 
             map_lgl(~ any(.x >= 366L & .x <= 523L)),
           parE_result = if_else(gene != "parE", NA, parE_result),
           parE_result = as.integer(parE_result)) %>%
    mutate(result_gyr_par = case_when(gene == "gyrA" ~ gyrA_result,
                                      gene == "gyrB" ~ gyrB_result,
                                      gene == "parC" ~ parC_result,
                                      gene == "parE" ~ parE_result)) %>%
    mutate(result_total = ifelse(gene %in% c("gyrA","gyrB","parC","parE"), result_gyr_par, result)) %>%
    select(-c(gyrA_result,
              gyrB_result,
              parC_result,
              parE_result,
              result_gyr_par,
              result))
  return(df)
}

# Function that returns a filtered dataframe based on the string "genes"
filter_in_table <- function(df) {
  report_genes <- unique(df$gene)
  grep_genes <- c()
  for (gene in in_genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  df <- df %>%
    filter(gene %in% grep_genes) %>%
    mutate(gene = as.character(gene))
  return(df)
}

# calculates percentage of present/absent mutations and genes
calc_stats <- function(df) {
  df <- df %>%
    group_by(gene, result_total) %>%
    count() %>%
    ungroup() %>%
    mutate(result_total = if_else(result_total == 1, "Present", "Absent")) %>%
    spread(result_total, n, fill = 0) %>%
    mutate(Absent = if ("Absent" %in% names(.)) {
      return(Absent)
    } else {
      return(0) # adds Absent column with all 0 if all isolates have genes
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

# Generate overview table of intrinsic genes
create_in_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    filter(flag %in% flag_selection) %>%
    mutate(id = 1:n()) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, flag)) %>%
    gather(gene, mut, -ref) %>%
    mutate(mut = ifelse(mut == "" | mut == "." | is.na(mut) == TRUE, 0, mut),
           result = ifelse(mut != 0, 1, 0),
           result = as.integer(result),
           type = "mut")
  return(df)
}        

# creates a data frame with one row per sample, and 1/0 results 
# for mutations in respective genes in columns
create_in_report <- function(df) {
  df <- df %>%
    select(-mut) %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, type))
  return(df)
}

# Creates a data frame with the actual mutations that were identified
# in the genes
create_mut_report <- function(df) {
  df <- df %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, mut) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, type, result_total))
  return(df)
}

# -------------------------- Analysis ----------------------------

# Create output directory
dir.create(paste0(output_loc, "/amr_in/"), showWarnings = FALSE)
amr_output <- paste0(output_loc, "/amr_in/")

## Intrinsic genes

in_data <- get_ariba_data(in_report_loc) %>%
  fix_gene_names()

in_flags <- check_flags(in_data)

in_table <- create_in_table(in_data) %>%
  fix_gyr_par_results()

if ("ALL" %in% in_genes) {
  in_table_filtered <- in_table
} else {
  in_table_filtered <- filter_in_table(in_table)
}

in_report <- create_in_report(in_table_filtered)
in_mut_report <- create_mut_report(in_table_filtered)
in_stats <- calc_stats(in_table_filtered)

## Write results to file

write.table(in_report,
            paste0(amr_output, "intrinsic_gene_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(in_mut_report,
            paste0(amr_output, "intrinsic_mut_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(in_flags,
            paste0(amr_output, "intrinsic_flag_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(in_stats,
            paste0(amr_output, "intrinsic_gene_stats.tsv"),
            sep = "\t",
            row.names = FALSE)