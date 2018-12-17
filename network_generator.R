#!/usr/bin/env Rscript
library(rjson)
library(argparse)
suppressPackageStartupMessages(library(GenomicRanges))
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
parser$add_argument("--nearest",
  action = "store_true",
  help = "Search the nearest range"
)
parser$add_argument("--expand",
  type = "integer",
  default = "0",
  help = "Number of bases to expand the search by range"
)


args <- parser$parse_args()
# Differents examples of parameters
#args <- parser$parse_args(c("~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "6:52155590-52158317", "--expand", "20000"))
# args <- parser$parse_args(c("~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "6:52155590-52158317"))
#args <- parser$parse_args(c("~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "6:52155590-52158317", "--nearest"))
#args <- parser$parse_args(c("~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt", "--search", "Hoxa1"))

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
chr <- c(chrs_wt$baitChr, chrs_wt$oeChr)
start <- c(chrs_wt$baitStart, chrs_wt$oeStart)
end <- c(chrs_wt$baitEnd, chrs_wt$oeEnd)
curated_chrs_vertex <- distinct(data_frame(fragment, gene_names, chr, start, end))
# Only use the first name
curated_gene_name <- str_split_fixed(curated_chrs_vertex$gene_names, ",", n = 2)[, 1]
# Remove from the last dash to the end of the name
curated_gene_name <- str_replace(curated_gene_name, "-[^-]+$", '')
curated_gene_name <- ifelse(curated_gene_name != '.', curated_gene_name, '')
curated_chrs_vertex <- mutate(curated_chrs_vertex, curated_gene_name)
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
  # Detect if we are searching by position (we are working with mouse chromosomes by now) or by name
  # Always return NULL if it doesn't exist the vertex in the graph
  if (str_detect(vertex, "^((1?[0-9])|([XY]))_\\d+$")) {
    if (!any(V(net)$name == vertex)) {
      return(NULL)
    }
    return(V(net)[vertex])
  } else {
    searched_vertex_index <- curated_chrs_vertex$curated_gene_name == vertex
    if (!any(searched_vertex_index)) {
      return(NULL)
    }
    return(V(net)[searched_vertex_index])
  }
}

## Generate required subnetwork
if (!is.null(args$search)) {
  if (str_detect(args$search, "((1?[0-9])|([XY])):\\d+(-\\d+)?$")) {
    # We are working with a range
    curated_chrs_vertex_ranges <- makeGRangesFromDataFrame(curated_chrs_vertex, keep.extra.columns = T, ignore.strand = FALSE)
    required_range <- GRanges(args$search)
    # Expand the selected range if it is required
    if (args$expand != 0) {
      start(required_range) <- start(required_range) - args$expand
      end(required_range) <- end(required_range) + args$expand
    }
    # Work with the nearest if it is required
    if (args$nearest) {
      nearest_range_index <- nearest(required_range, curated_chrs_vertex_ranges)
      required_vertex <- curated_chrs_vertex_ranges[nearest_range_index]$fragment
      # make_ego_graph always returns a list
      required_subnet <- make_ego_graph(net, nodes = required_vertex)[[1]]
    } else {
      # Work with overlaps instead
      overlaps_index <- subjectHits(findOverlaps(required_range, curated_chrs_vertex_ranges))
      required_vertex <- curated_chrs_vertex_ranges[overlaps_index]$fragment
      required_vertex_with_neighbours <- names(unlist(lapply(required_vertex, function(rv) {
        neighbors(net, rv)
      })))
      # Add the overlapping vertex
      required_vertex_with_neighbours <- unique(c(required_vertex_with_neighbours, required_vertex))
      required_subnet <- induced_subgraph(net, vids = required_vertex_with_neighbours)
    }
  } else {
    required_vertex <- search_vertex(args$search, net)
    # make_ego_graph always returns a list
    required_subnet <- make_ego_graph(net, nodes = required_vertex)[[1]]
  }
  if (is.null(required_vertex)) {
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
# Convert the required subnetwork to Cytoscape Json format
if (!is.null(required_subnet)) {
  generate_cytoscape_json(required_subnet)
} else {
  cat("[]")
}

# Plotting example
#plot(required_subnet, vertex.label = V(required_subnet)$curated_gene_name, vertex.size = 5 + 2 * degree(required_subnet), vertex.color = c("gray", "lightgreen")[1 + V(required_subnet)$EZH2], vertex.shape = ifelse(V(required_subnet)$type == "bait", "circle", "square"), edge.color = "black")

