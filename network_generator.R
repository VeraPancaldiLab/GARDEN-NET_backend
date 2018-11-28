#!/usr/bin/env Rscript
library(rjson)
library(argparse)
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(tidyverse))

# Argparse
parser <-
  ArgumentParser(description = "Separated values file to cytoscape json mapper")
parser$add_argument("PCHiC file", nargs = 1, help = "Separated values file PCHiC as input file")
parser$add_argument("--wt_threshold",
  type = "double",
  default = "5.0",
  help = "The minimun value for considering the edge"
)
parser$add_argument("--features",
  help = "Separated values file of features as input file"
)
parser$add_argument("--search",
  help = "Search node by name or fragment position in the graph to generate a neighborhood subgraph"
)
parser$add_argument("--chromosome",
  help = "Filter by chromosome"
)
parser$add_argument("--no-features-binarization",
  action = "store_true",
  help = "Features will be binarized by default"
)

# args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
args <- parser$parse_args()

# Load PCHiC
chrs <-
  suppressMessages(read_tsv(
    file = args$`PCHiC file`,
    col_types = cols(baitChr = col_character())
  ))
# Filter by threshold
chrs_wt <- chrs[chrs$mESC_wt > args$wt_threshold, ]

# Filter by chromosome
if (!is.null(args$chromosome)) {
  chrs_wt <- chrs_wt[which(chrs_wt$baitChr == args$chromosome
  & chrs_wt$oeChr == args$chromosome), ]
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
curated_chrs_vertex <- distinct(data_frame(fragment, gene_names))
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
curated_chrs_edges <- data.frame(source = baits, target = oes)

## ------------------------------------------------------------------------
net <-
  graph_from_data_frame(curated_chrs_edges, directed = F, curated_chrs_vertex)

## ------------------------------------------------------------------------
search_vertex <- function(vertex, graph) {
  # Detect if we are searching by positiona (we are working with mouse chromosomes by now) o by name
  # Always return NULL if it doesn't exist the vertex in the graph
  if (str_detect(vertex, "^((1?[0-9])|([XY]))_\\d+$")) {
    if (!any(V(net)$name == vertex)) {
      return(NULL)
    }
    return(V(net)[vertex])
  } else {
    searched_vertex_index <- sapply(
      V(net)$gene_list,
      function(gene_list)
        if (vertex %in% gene_list) {
          T
        } else {
          F
        }
    )
    if (!any(searched_vertex_index)) {
      return(NULL)
    }
    return(V(net)[searched_vertex_index])
  }
}

## Generate required subnetwork
if (!is.null(args$search)) {
  required_vertex <- search_vertex(args$search, net)
  if (!is.null(required_vertex)) {
    # make_ego_graph always returns a list
    required_subnet <- make_ego_graph(net, nodes = required_vertex)[[1]]
  } else {
    required_subnet <- NULL
  }
} else {
  required_subnet <- net
}


## Generate Cytoscape JSON
generate_cytoscape_json <- function(required_subnet) {
  # Recover vertices neighbourhood subnetwork dataframe from the graph
  vertices_df <-
    igraph::as_data_frame(required_subnet, what = "vertices")
  # Remove row names
  row.names(vertices_df) <- NULL
  # _ in column names is not valid in Cytoscape JSON
  vertices_df <-
    rename(vertices_df, !!c(id = "name", names = "gene_names"))
  # lists are not a valid supported type in Cytoscape JSON
  vertices_df$gene_list <- NULL
  # Nest all vertice rows inside data key and add the group type, both required by Cytoscape JSON
  vertices_df <-
    apply(vertices_df, 1, function(vertice_row) {
      list(data = vertice_row, group = "nodes")
    })
  # Recover edges neighbourhood subnetwork dataframe from the graph
  edges_df <- igraph::as_data_frame(required_subnet, what = "edges")
  # Remove row names
  row.names(edges_df) <- NULL
  # Rename edge extremes with the Cytoscape JSON squema
  edges_df <- rename(edges_df, !!c(source = "from", target = "to"))
  # Add id to the edges
  edges_df$id <- paste(edges_df$source, edges_df$target, sep = "~")
  # Nest all edge rows inside data key and add the group type, both required by Cytoscape JSON
  edges_df <-
    apply(edges_df, 1, function(edge_row) {
      list(data = edge_row, group = "edges")
    })
  # Join vertices and edges inside the same dataframe
  JSON_df <- c(vertices_df, edges_df)
  # Write JSON dataframe to a file
  # write_lines(toJSON(JSON_df, indent = 2), 'neighboord.json')
  cat(toJSON(JSON_df, indent = 2))
}
if (!is.null(required_subnet)) {
  generate_cytoscape_json(required_subnet)
} else {
  cat("[]")
}

# plot(required_subnet, vertex.label = ifelse(V(required_subnet)$name == required_vertex$name, required_vertex$name, ''))