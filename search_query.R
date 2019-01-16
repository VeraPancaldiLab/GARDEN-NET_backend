#!/usr/bin/env Rscript
library(rjson)
library(optparse)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(tidyverse))

source("./network_generator_lib.R")

args <- commandArgs(trailingOnly = TRUE)
args <- parser_arguments(args)

# Load Rdata
load('garnet.Rdata')

# Search the required subnetwork
required_subnet <- search_subnetwork(args$search, args$expand, args$nearest, net, curated_chrs_vertex)

# Convert the required subnetwork to Cytoscape Json format
required_subnet_json <- generate_cytoscape_json(required_subnet)
cat(required_subnet_json)
