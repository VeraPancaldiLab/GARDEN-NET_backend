#!/usr/bin/Rscript
library(optparse)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(igraph))

source("./network_generator_lib.R")

args <- commandArgs(trailingOnly = TRUE)
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt","--search", "1_173143867"))
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt","--search", "6:52155590-52158317"))
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt","--search", "asdfasdfa"))
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "Hoxa1"))
args <- parser_arguments(args)

PCHiC <- load_PCHiC(args$PCHiC)

PCHiC <- filter_by_threshold(PCHiC, args$wt_threshold)

if (!is.null(args$chromosome)) {
  PCHiC <- filter_by_chromosome(PCHiC, args$chromosome)
}

curated_PCHiC_vertex <- generate_vertex(PCHiC)

# Finally add all features to their corresponding fragments
if (!is.null(args$features)) {
  curated_PCHiC_vertex <- generate_features(args$features)
}

curated_PCHiC_edges <- generate_edges(PCHiC)

## ------------------------------------------------------------------------
# Generate the network
net <-
  graph_from_data_frame(curated_PCHiC_edges, directed = F, curated_PCHiC_vertex)

# Add additional network metadata
V(net)$total_degree <- degree(net)

# Search the required subnetwork
if (!is.null(args$search)) {
  if (str_detect(args$search, "((1?[0-9])|([XYxy])):\\d+(-\\d+)?$")) {
    # Only load GenomicRanges if the string searched is a range
    suppressPackageStartupMessages(library(GenomicRanges))
  }

  required_subnet <- search_subnetwork(args$search, args$expand, args$nearest, net, curated_PCHiC_vertex)
} else {
  required_subnet <- net
}

# Convert the required subnetwork to Cytoscape Json format
if (is.null(required_subnet)) {
  cat("{}")
} else {
  library(rjson)

  # Convert the required subnetwork to Cytoscape Json format
  required_subnet_json <- generate_cytoscape_json(required_subnet)
  if (is.null(args$pipeline)) {
    cat(required_subnet_json)
  }
}

# Plotting example
# plot(
#  required_subnet,
#  vertex.label = V(required_subnet)$curated_gene_name,
#  vertex.size = 25 + 2 * degree(required_subnet),
#  vertex.color = c("gray", "lightgreen")[1 + V(required_subnet)$EZH2],
#  vertex.shape = ifelse(V(required_subnet)$type == "bait", "square", "circle"),
#  edge.color = "black"
# )
