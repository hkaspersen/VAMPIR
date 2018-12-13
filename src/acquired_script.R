# --------------- AMR gene analysis - Acquired genes --------------

## This track of the ARIBA analysis script analyses AMR gene 
## reports from ARIBA and generates result files based on user
## input of selected genes of interest.

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
ac_report_loc <- args[1]
output_loc <- args[2]
ac_genes <- args[3]

# adjust parameters for filtering
if (grepl("all", ac_genes, ignore.case = TRUE) == TRUE) {
  ac_genes <- "ALL"
} else {
  ac_genes <- unlist(strsplit(ac_genes, ",", fixed = TRUE))
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
  data <- bind_rows(lapply(data_list, function(x) map(x, as.character)), .id = "ref")
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
    p <- paste(tolower(substring(i, 1,3)), substring(i, 4), sep = "", collapse = " ")
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
           flag_result = as.integer(flag_result)) %>%
    rename("gene" = gene_names)
  return(df)
}

# Function that returns a filtered dataframe based on the string "acquired_genes"
filter_ac_table <- function(df) {
  report_genes <- unique(df$gene)
  grep_genes <- c()
  for (gene in ac_genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  df <- df %>%
    filter(gene %in% grep_genes) %>%
    mutate(gene = gsub("_", "", gene),
           gene = as.character(gene))
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

# Generate overview table of acquired genes
create_ac_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_amr_report.tsv$", "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, result,-ref) %>%
    mutate(result_total = ifelse(result == "", 0, 1),
           result_total = as.character(result_total),
           type = "gene") %>%
    select(-result)
  return(df)
}

# creates a data frame with one row per sample, and 1/0 results
# for acquired genes in columns
create_ac_report <- function(df) {
  df <- df %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, type))
  return(df)
}
# -------------------------- Analysis ----------------------------

# Create output directory
dir.create(paste0(output_loc, "/amr_ac/"), showWarnings = FALSE)
amr_output <- paste0(output_loc, "/amr_ac/")

## Acquired genes

ac_data <- get_ariba_data(ac_report_loc) %>%
  fix_gene_names()

ac_flags <- check_flags(ac_data)

ac_table <- create_ac_table(ac_data)

if ("ALL" %in% ac_genes) {
  ac_table_filtered <- ac_table
} else {
  ac_table_filtered <- filter_ac_table(ac_table)
}

ac_report <- create_ac_report(ac_table_filtered)
ac_stats <- calc_stats(ac_table_filtered)

## Write results to file

write.table(ac_report,
            paste0(amr_output, "acquired_gene_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(ac_flags,
            paste0(amr_output, "acquired_flag_report.tsv"),
            sep = "\t",
            row.names = FALSE)

write.table(ac_stats,
            paste0(amr_output, "acquired_gene_stats.tsv"),
            sep = "\t",
            row.names = FALSE)