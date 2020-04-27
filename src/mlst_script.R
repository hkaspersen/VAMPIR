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
    "ggplot2",
    "dplyr",
    "tidyr",
    "stringr",
    "svglite",
    "phangorn",
    "ggtree",
    "tibble",
    "purrr",
    "cluster",
    "ape",
    "impoRt"
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
                      "mlst_summarized_results.tsv",
                      convert = TRUE) %>%
  select(-ref)

mlst_genes <- unlist(strsplit(names(mlst_data)[-1],
                              split = " ",
                              fixed = TRUE))

mlst_genes <- mlst_genes[-1]

split_column <- names(mlst_data)[2]

allele_matrix <- suppressWarnings(mlst_data %>%
  mutate_at(vars(-1),
            list(trimws)) %>%
  tidyr::separate(split_column, into = mlst_genes, extra = "merge") %>%
  rename("ref" = header) %>%
  { . ->> mlst_report } %>%
  select(-ST) %>%
  mutate_at(vars(-ref),
            list(sub("*", "", .))) %>%
  mutate_at(vars(-ref),
            list(as.integer(.))) %>% 
  mutate(test = complete.cases(.)) %>%
  filter(test == TRUE) %>%
  select(-test) %>%
  column_to_rownames("ref"))

# Calculate distance matrix from sequence typing alleles and create tree
tree <- as.phylo(hclust(daisy(allele_matrix,
                              metric = "gower"),
                        method = "average"))

tree$tip.label <- rownames(allele_matrix)

p <- suppressWarnings(ggtree(tree) +
  geom_treescale() +
  geom_tiplab(size = 1,
              align = TRUE))

# Save output
ggsave(paste0(mlst_output, "mlst_tree.svg"),
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

write.tree(phy = tree, file = paste0(mlst_output,"mlst_tree.newick"))