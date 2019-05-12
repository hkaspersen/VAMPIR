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

# adjust parameters for filtering
if (grepl("all", vir_genes, ignore.case = TRUE) == TRUE) {
  vir_genes <- "ALL"
} else {
  vir_genes <- unlist(strsplit(vir_genes, ",", fixed = TRUE))
}

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

# Corrects the gene names found in the "ref_name" column
# for data from the virfinder database
fix_virfinder_names <- function(df) {
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

# Corrects the gene names found in the "ref_name" column
# for data from the vfdb database
fix_vfdb_names <- function(df) {
  genes <- unique(df$ref_name)
  new_names <- sub("\\..+", "", genes)
  
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
           flag_result = as.integer(flag_result)) %>%
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
    summarise_all(list(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_vir_report.tsv$",
                      "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, result,-ref) %>%
    mutate(result_total = ifelse(result == "", 0, 1),
           result_total = as.character(result_total)) %>%
    select(-result)
  return(df)
}

# Handle data from virfinder
create_virfinder_table <- function(df) {
  df <- df %>%
    spread(`Virulence factor`, Identity) %>%
    group_by(ref) %>%
    summarise_all(list(func_paste)) %>%
    mutate(ref = sub("/results_tab.tsv", "", ref)) %>%
    select(-c(Database,
              `Query / Template length`,
              Contig,
              `Position in contig`,
              `Protein function`,
              `Accession number`,
              `<NA>`)) %>%
    mutate_at(vars(-ref),
              list(~ifelse(. == "", 0, 1))) %>%
    gather(gene, value, -ref)
  
  return(df)
}

# Filter genes in table based on user input
filter_vir_table <- function(df) {
  report_genes <- unique(df$gene)
  grep_genes <- c()
  for (gene in vir_genes) {
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
    mutate(Absent = if ("Absent" %in% names(.)){return(Absent)}else{return(0)}) %>%
    rowwise() %>%
    mutate(
      Total = Present + Absent,
      Percent = round(Present / Total * 100, 1),
      lwr = round(get_binCI(Present, Total)[1], 1),
      upr = round(get_binCI(Present, Total)[2], 1)
    )
  return(df)
}

calc_summary_stats <- function(df) {
  df <- df %>%
    gather(gene, result_total, -ref) %>%
    group_by(gene, result_total) %>%
    count() %>%
    ungroup() %>%
    mutate(result_total = if_else(result_total == 1, 
                                  "Present",
                                  "Absent")) %>%
    spread(result_total, n, fill = 0) %>%
    mutate(Absent = if ("Absent" %in% names(.)){return(Absent)}else{return(0)}) %>%
    rowwise() %>%
    mutate(
      Total = Present + Absent,
      Percent = round(Present / Total * 100, 1),
      lwr = round(get_binCI(Present, Total)[1], 1),
      upr = round(get_binCI(Present, Total)[2], 1)
    )
}

calc_virfinder_stats <- function(df) {
  df <- df %>%
    group_by(gene, value) %>%
    count() %>%
    ungroup() %>%
    mutate(value = if_else(value == 1, 
                           "Present",
                           "Absent")) %>%
    spread(value, n, fill = 0) %>%
    mutate(Absent = if ("Absent" %in% names(.)){return(Absent)}else{return(0)}) %>%
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
    summarise_all(list(func_paste)) %>%
    select(-id)
  return(df)
}

create_summary_report <- function(df) {
  df <- df %>%
    mutate(gene = gsub("[0-9]", "", gene)) %>%
    group_by(ref, gene, result_total) %>%
    summarise_all(list(func_paste)) %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(list(func_paste)) %>%
    select(-id) %>%
    mutate_at(vars(-ref),
              list(~ifelse(. == "0, 1", "1", .)))
  return(df)
}

# --------------------------- Analysis ---------------------------

# Create output directory
dir.create(paste0(output_loc, "/vir/"), showWarnings = FALSE)
vir_output <- paste0(output_loc, "/vir/")

# Clean data
if (vir_database == "virfinder") {
  # import data
  vir_data <- get_data(report_loc,
                       "vir_report.tsv",
                       convert = TRUE)
  
  clean_vir_data <- fix_virfinder_names(vir_data)
  vir_flags <- check_flags(clean_vir_data)
  vir_table <- create_vir_table(clean_vir_data)
  
  if ("ALL" %in% vir_genes) {
    vir_table_filtered <- vir_table
  } else {
    vir_table_filtered <- filter_vir_table(vir_table)
  }
  
  vir_report <- create_vir_report(vir_table_filtered)
  vir_summary <- create_summary_report(vir_table_filtered)
  vir_quant_detailed <- calc_stats(vir_table_filtered)
  vir_quant_summary <- calc_summary_stats(vir_summary)
  
  write.table(vir_report,
              paste0(vir_output, "virulence_report_detailed.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_summary,
              paste0(vir_output, "virulence_report_summary.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_quant_summary,
              paste0(vir_output, "virulence_stats_summary.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_quant_detailed,
              paste0(vir_output, "virulence_stats_detailed.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_flags,
              paste0(vir_output, "virulence_flags.tsv"),
              sep = "\t",
              row.names = FALSE)
}

if (vir_database %in% c("vfdb", "vfdb_core")) {
  # Import data
  vir_data <- get_data(report_loc,
                       "vir_report.tsv",
                       convert = TRUE)
  
  clean_vir_data <- fix_vfdb_names(vir_data)
  vir_flags <- check_flags(clean_vir_data)
  vir_table <- create_vir_table(clean_vir_data)
  
  if ("ALL" %in% vir_genes) {
    vir_table_filtered <- vir_table
  } else {
    vir_table_filtered <- filter_vir_table(vir_table)
  }
  
  vir_report <- create_vir_report(vir_table_filtered)
  vir_quant_detailed <- calc_stats(vir_table_filtered)
  
  write.table(vir_report,
              paste0(vir_output, "virulence_report_detailed.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_quant_detailed,
              paste0(vir_output, "virulence_stats_detailed.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_flags,
              paste0(vir_output, "virulence_flags.tsv"),
              sep = "\t",
              row.names = FALSE)
}

if (vir_database == "virfinderDtu") {
  
  vir_data <- get_data(report_loc,
                       "results_tab.tsv")
  
  vir_table <- create_virfinder_table(vir_data)
  
  if ("ALL" %in% vir_genes) {
    vir_table_filtered <- vir_table
  } else {
    vir_table_filtered <- filter_vir_table(vir_table)
  }
  
  vir_quant <- calc_virfinder_stats(vir_table_filtered)
  presence_absence <- spread(vir_table_filtered, gene, value)
  
  
  write.table(presence_absence,
              paste0(vir_output, "virfinder_results.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_quant,
              paste0(vir_output, "virfinder_stats.tsv"),
              sep = "\t",
              row.names = FALSE)
  
  write.table(vir_data,
              paste0(vir_output, "virfinder_results_full.tsv"),
              sep = "\t",
              row.names = FALSE)
  
}