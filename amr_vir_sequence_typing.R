#!/usr/bin/env Rscript

###########################################################
# AMR, Virulence, Plasmid- and Multilocus sequence typing #
###########################################################

# This script takes multiple input from the software
# ARIBA (github.com/sanger-pathogens/ariba) and
# create summary reports and visualizations. One
# can select specific analyses based on the data you
# have: AMR, virulence, MLST or plasmids.

# Author: Håkon Kaspersen

##########################################################
##########################################################
##########################################################

# ------------------- Libraries --------------------
suppressPackageStartupMessages(pacman::p_load(optparse))
parser <- OptionParser(usage = "Usage: %prog [options]")

# Create command line options
parser <- add_option(parser,
                     c("-u", "--mut"),
                     action = "store",
                     help = "Location of intrinsic gene reports.")
parser <- add_option(parser,
                     c("-a", "--acq"),
                     action = "store",
                     help = "Location of acquired gene reports.")
parser <- add_option(parser,
                     c("-i", "--intrinsic"),
                     action = "store",
                     help = "List of intrinsic genes of interest. 
                     Type 'all' for including all reported genes.")
parser <- add_option(parser,
                     c("-c", "--acquired"),
                     action = "store",
                     help = "List of acquired genes of interest.
                     Type 'all' for including all reported genes.")
parser <- add_option(parser,
                     c("-v", "--vir"),
                     action = "store",
                     help = "Location of ARIBA virulence reports.")
parser <- add_option(parser,
                     c("-m", "--mlst"),
                     action = "store",
                     help = "Location of ARIBA MLST reports.")
parser <- add_option(parser,
                     c("-p", "--plasmid"),
                     action = "store",
                     help = "Location of ARIBA plasmid reports.")
parser <- add_option(parser,
                     c("-o", "--output"),
                     action = "store",
                     help = "Output directory location. 
                     One folder for each analysis will be created
                     at given location.")
opt <- parse_args(parser)

# Check if output folder is specified
if (is.null(opt$output)) {
  print("Please specify an output directory.")
  stop()
}

## ------------------- Tracks ----------------------
## Intrinsic AMR genes track
if (!is.null(opt$mut)) {
  if (is.null(opt$intrinsic)) {
    print("Please specify genes of interest.")
    stop()
  } else {
    print(paste0(
      "Running intrinsic AMR gene summary analysis. Reports location: ",
    opt$mut,
    ". Output location: ",
    opt$out))
  system(paste("Rscript /work/projects/nn9305k/vi_src/amr_vir_sequence_typing/intrinsic_script.R",
               opt$mut, 
               opt$out, 
               opt$intrinsic))
  }
}

## Acquired AMR genes track
if (!is.null(opt$acq)) {
  if (is.null(opt$acquired)) {
    print("Please specify genes of interest.")
    stop()
  } else {
    print(paste0(
      "Running acquired AMR gene summary analysis. Reports location: ",
    opt$acq,
    ". Output location: ",
    opt$out))
  system(paste("Rscript /work/projects/nn9305k/vi_src/amr_vir_sequence_typing/acquired_script.R",
               opt$acq, 
               opt$out, 
               opt$acquired))
  }
}

## Virulence gene track
if (!is.null(opt$vir)) {
  print(paste0(
    "Running virulence gene summary analysis. Reports location: ",
    opt$vir,
    ". Output location: ",
    opt$out))
  system(paste("Rscript /work/projects/nn9305k/vi_src/amr_vir_sequence_typing/vir_script.R", opt$vir, opt$out))
}

## MLST track
if (!is.null(opt$mlst)) {
  print(paste0(
    "Running MLST summary analysis. Reports location: ", 
    opt$mlst,
    ". Output location: ",
    opt$out))
  system(paste("Rscript /work/projects/nn9305k/vi_src/amr_vir_sequence_typing/mlst_script.R", opt$mlst, opt$out))
}

## Plasmid typing track
if (!is.null(opt$plasmid)) {
  print(paste0(
    "Running plasmid summary analysis. Reports location: ",
    opt$plasmid,
    ". Output location: ",
    opt$out))
  system(paste("Rscript /work/projects/nn9305k/vi_src/amr_vir_sequence_typing/plasmid_script.R", opt$plasmid, opt$out))
}
