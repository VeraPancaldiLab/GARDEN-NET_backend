#!/usr/bin/env Rscript
library(rjson)
library(optparse)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(tidyverse))

source("./network_generator_lib.R")

# Argparse
parser <- OptionParser(description = "Separated values file to cytoscape json mapper")
parser <- add_option(parser, "--PCHiC", help = "Separated values file PCHiC as input file")
parser <- add_option(parser, "--wt_threshold",
                    type = "double",
                    default = 5.0,
                    help = "The minimun value for considering the edge [default: %default]")
parser <- add_option(parser, "--features",
                    help = "Separated values file of features as input file")
parser <- add_option(parser, "--search",
                    help = "Search node by name or fragment position in the graph to generate a neighborhood subgraph")
parser <- add_option(parser, "--chromosome",
                    help = "Filter by chromosome")
parser <- add_option(parser, "--no-features-binarization",
                    action = "store_true",
                    default = F,
                    help = "Features will be binarized by default")
parser <- add_option(parser, "--nearest",
                    action = "store_true",
                    default = F,
                    help = "Search the nearest range")
parser <- add_option(parser, "--expand",
                    type = "integer",
                    default = 0,
                    help = "Number of bases to expand the search by range")

# Differents examples of parameters
# args <- c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt")
# args <- c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "6:52155590-52158317", "--expand", "20000")
# args <- c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "6:52155590-52158317")
# args <- c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "6:52155590-52158317", "--nearest")
#args <- c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "Hoxa1")
# args <- c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "C10orf54")
# args <- c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "Hoxa13")
args <- parse_args(parser, convert_hyphens_to_underscores = T)

# Load PCHiC
chrs <-
  suppressMessages(read_tsv(
    file = args$PCHiC,
    col_types = cols(baitChr = col_character())
  ))

# Filter by threshold
chrs_wt <- chrs[chrs$mESC_wt > args$wt_threshold,]

# Filter by chromosome
if (!is.null(args$chromosome)) {
  chrs_wt <- chrs_wt[which(chrs_wt$baitChr == args$chromosome
                           & chrs_wt$oeChr == args$chromosome),]
}

## ------------------------------------------------------------------------
# Join chr number with the start position
# Also join bait and oe in the same column
fragment <- c(
  paste(chrs_wt$baitChr, chrs_wt$baitStart, sep = "_"),
  paste(chrs_wt$oeChr, chrs_wt$oeStart, sep = "_")
)
# Extract bait and oe names and join them to the same column
gene_names <- c(chrs_wt$baitName, chrs_wt$oeName)
# Remove repeted nodes
chr <- c(chrs_wt$baitChr, chrs_wt$oeChr)
start <- c(chrs_wt$baitStart, chrs_wt$oeStart)
end <- c(chrs_wt$baitEnd, chrs_wt$oeEnd)
curated_chrs_vertex <-
  distinct(tibble(fragment, gene_names, chr, start, end))
# Only use the first name
curated_gene_name <-
  str_split_fixed(curated_chrs_vertex$gene_names, ",", n = 2)[, 1]
# Remove from the last dash to the end of the name
curated_gene_name <- str_replace(curated_gene_name, "-[^-]+$", "")
curated_gene_name <-
  ifelse(curated_gene_name != ".", curated_gene_name, "")
curated_chrs_vertex <-
  mutate(curated_chrs_vertex, curated_gene_name)
# Add gene names in array form splitted by ,
curated_chrs_vertex$gene_list <-
  str_split(curated_chrs_vertex$gene_names, ",")
# Add the type for bait and oes
# Be careful because in the oe column there are many baits
# So oe is only oe if it not exist in the bait column
baits <- paste(chrs_wt$baitChr, chrs_wt$baitStart, sep = "_")
curated_chrs_vertex$type <-
  ifelse(curated_chrs_vertex$fragment %in% baits, "bait", "oe")
# Finally all all features to their corresponding fragments
# Load the features
if (!is.null(args$features)) {
  features <-
    suppressMessages(read_tsv(file = args$features))
  # Remove chr prefix from the fragment column
  features$fragment <- str_sub(features$fragment, start = 4)
  # Binarize all the features
  if (!args$no_features_binarization) {
    features[, -1] <- ifelse(features[, -1] == 0.0, 0, 1)
  }
  curated_chrs_vertex <-
    left_join(curated_chrs_vertex, features, by = "fragment")
}
## ------------------------------------------------------------------------
# Generate a dataframe with the extremes of the edges
oes <- paste(chrs_wt$oeChr, chrs_wt$oeStart, sep = "_")
curated_chrs_edges <- tibble(source = baits, target = oes)

## ------------------------------------------------------------------------
# Generate the network
net <-
  graph_from_data_frame(curated_chrs_edges, directed = F, curated_chrs_vertex)


# Search the required subnetwork
required_subnet <- search_subnetwork(args$search, args$expand, args$nearest, net, curated_chrs_vertex)

# Convert the required subnetwork to Cytoscape Json format
required_subnet_json <- generate_cytoscape_json(required_subnet)
cat(required_subnet_json)

# Plotting example
#plot(
#  required_subnet,
#  vertex.label = V(required_subnet)$curated_gene_name,
#  vertex.size = 25 + 2 * degree(required_subnet),
#  vertex.color = c("gray", "lightgreen")[1 + V(required_subnet)$EZH2],
#  vertex.shape = ifelse(V(required_subnet)$type == "bait", "square", "circle"),
#  edge.color = "black"
#)
