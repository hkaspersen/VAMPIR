# ------------------------- MLST analysis -------------------------

## This track of the ariba analysis script analyses the MLST
## reports from ARIBA and gives an UPGMA tree based on the 
## distances betweeen alleles, as well as an overall report of 
## sequence types.

# ------------------------- Parameters ----------------------------

args <- commandArgs(trailingOnly = TRUE)
report_loc <- args[1]
output_loc <- args[2]
ending <- args[3]

# ------------------------ Load libraries -------------------------

packages <-
  c(
    "dplyr",
    "tidyr",
    "ggtree",
    "tibble",
    "ape",
    "impoRt",
    "distanceR",
    "vampfunc"
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
dir.create(paste0(output_loc, "/mlst/"), showWarnings = FALSE)
mlst_output <- paste0(output_loc, "/mlst/")

# Import and analyse data
mlst_data <- get_data(report_loc,
                      ending,
                      convert = TRUE) %>%
  select(-ref)

perc_ST <- percent_presence(mlst_data, "ST")

allele_matrix <- mlst_data %>%
  select(-ST) %>%
  mutate_at(vars(-header),
            ~as.factor(as.character(.))) %>%
  column_to_rownames("header")

# Calculate distance matrix from sequence typing alleles and create tree
tree <- calc_dist(allele_matrix)

p <- plot_tree(tree)

# Save output
ggsave(paste0(mlst_output, "mlst_tree.png"),
       p,
       device = "png",
       dpi = 300,
       height = 20,
       width = 22,
       units = "cm")

write.table(mlst_report,
            paste0(mlst_output, "mlst_report.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(perc_ST,
            paste0(mlst_output, "mlst_percent.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.tree(phy = tree, file = paste0(mlst_output,"mlst_tree.newick"))