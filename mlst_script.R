# ------------------------- MLST analysis -------------------------

## This track of the ariba analysis script analyses the MLST
## reports from ARIBA and gives an UPGMA tree based on the 
## distances betweeen alleles, as well as an overall report of 
## sequence types.

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
report_loc <- args[1]
output_loc <- args[2]

# ------------------------ Load libraries -------------------------

suppressPackageStartupMessages(
  if (!require("pacman")) install.packages("pacman"))
suppressPackageStartupMessages(
  pacman::p_load(
    ggplot2,
    dplyr,
    tidyr,
    stringr,
    svglite,
    phangorn,
    ggtree,
    tibble,
    purrr
  )
)

# -------------------------- Functions ----------------------------

file_names_mlst <- function(filepath) {
  files <- list.files(path = filepath,
                      pattern = "mlst_summarized_results.tsv")
  return(files)
}

get_mlst_data <- function(filepath) {
  files <- file_names_mlst(filepath)
  
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
    )
  )
  return(data)
}


# -------------------------- Analysis ----------------------------

# Create output folder
dir.create(paste0(output_loc, "/mlst/"), showWarnings = FALSE)
mlst_output <- paste0(output_loc, "/mlst/")

# Import and analyse data
mlst_data <- get_mlst_data(report_loc)

mlst_genes <- unlist(strsplit(names(mlst_data)[-1],
                              split = ".",
                              fixed = TRUE))

split_column <- names(mlst_data)[2]

allele_matrix <- suppressWarnings(mlst_data %>%
  mutate_at(vars(-1),
            funs(trimws)) %>%
  tidyr::separate(split_column, into = mlst_genes, extra = "merge") %>%
  rename("ref" = header) %>%
  { . ->> mlst_report } %>%
  select(-ST) %>%
  mutate_at(vars(-ref),
            funs(sub("*", "", .))) %>%
  mutate_at(vars(-ref),
            funs(as.integer(.))) %>% 
  mutate(test = complete.cases(.)) %>%
  filter(test == TRUE) %>%
  select(-test) %>%
  column_to_rownames("ref"))

# Calculate distance matrix from sequence typing alleles
d.mlst.distances <- matrix(0, ncol=nrow(allele_matrix), 
                           nrow=nrow(allele_matrix))
for (i in 1:(nrow(allele_matrix)-1)) {
  for (j in (i+1):nrow(allele_matrix)){
    d.mlst.distances[i,j] <- sum(
      allele_matrix[i, ] != allele_matrix[j, ])
    d.mlst.distances[j,i] <- sum(
      allele_matrix[i, ] != allele_matrix[j, ])
  }
}

# Create tree from distance matrix
tree <- upgma(d.mlst.distances)

tree$tip.label <- rownames(allele_matrix)

p <- suppressWarnings(ggtree(tree) +
  geom_treescale() +
  geom_tiplab(size = 1,
              align = TRUE))

# Save output
ggsave(paste0(mlst_output, "NJTree_mlst.svg"),
       p,
       device = "svg",
       dpi = 300,
       height = 30,
       width = 25,
       units = "cm")

write.table(mlst_report,
            paste0(mlst_output, "mlst_report.tsv"),
            sep = "\t",
            row.names = FALSE)
