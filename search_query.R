#!/usr/bin/Rscript
library(optparse)

source("./network_generator_lib.R")

args <- commandArgs(trailingOnly = TRUE)
# args <- parser_arguments(args = c('--search', '1_173143867'))
# args <- parser_arguments(args = c('--search', '6:52155590-52158317'))
# args <- parser_arguments(args = c('--search', 'asdfasdfa'))
# args <- parser_arguments(args = c('--search', 'Hoxa1'))
args <- parser_arguments(args)

# Load Rdata
load("garnet.Rdata")

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(igraph))

if (str_detect(args$search, "((1?[0-9])|([XY])):\\d+(-\\d+)?$")) {
  # Only load GenomicRanges if the string searched is a range
  suppressPackageStartupMessages(library(GenomicRanges))
}

# Search the required subnetwork
required_subnet <- search_subnetwork(args$search, args$expand, args$nearest, net, curated_chrs_vertex)

if (is.null(required_subnet)) {
  cat("{}")
} else {
  library(rjson)
  # Convert the required subnetwork to Cytoscape Json format
  required_subnet_json <- generate_cytoscape_json(required_subnet)
  cat(required_subnet_json)
}
