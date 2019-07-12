#!/usr/bin/Rscript
library(optparse)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(chaser))

# Never use scientific format for numbers
options(scipen = 999)

source("./network_generator_lib.R")

args <- commandArgs(trailingOnly = TRUE)
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt","--search", "1_173143867"))
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt","--search", "6:52155590-52158317"))
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC.txt","--search", "asdfasdfa"))
# args <- parser_arguments(args = c("--PCHiC", "~/R_DATA/ChAs/PCHiC_interaction_map.txt", "--features", "~/R_DATA/ChAs/Features_mESC_ALL.tsv", "--search", "Hoxa1"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Mus_musculus-Embryonic_stem_cells.tsv", "--features", "./input_datasets/Mus_musculus-Embryonic_stem_cells.features", "--chromosome", "2"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Mus_musculus-Embryonic_stem_cells.tsv", "--features", "./input_datasets/Mus_musculus-Embryonic_stem_cells.features", "--only_pp_interactions"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Homo_sapiens-aCD4.tsv", "--chromosome", "1"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Homo_sapiens-aCD4.tsv"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Homo_sapiens-Mon.tsv", "--alias", "./alias_databases/Homo_sapiens.tsv", "--intronic_regions", "intronic_regions.tsv"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Homo_sapiens-Neu.tsv", "--alias", "./alias_databases/Homo_sapiens.tsv", "--intronic_regions", "intronic_regions.tsv", "--bait_names", "./HindIII_annotation_ens37.txt"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Homo_sapiens-Mon.tsv", "--alias", "./alias_databases/Homo_sapiens.tsv", "--intronic_regions", "intronic_regions.tsv", "--bait_names", "./HindIII_annotation_ens37.txt"))
# args <- parser_arguments(args = c("--PCHiC", "../GARDEN-NET_utils_related/PCHiC_peak_matrix_cutoff5.tsv", "--alias", "./alias_databases/Homo_sapiens.tsv", "--intronic_regions", "intronic_regions.tsv", "--bait_names", "./HindIII_annotation_ens37.txt"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Mus_musculus-Embryonic_stem_cells.tsv", "--features", "./input_datasets/Mus_musculus-Embryonic_stem_cells.features", "--chromosome", "1"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Mus_musculus-Embryonic_stem_cells.tsv", "--features", "./input_datasets/Mus_musculus-Embryonic_stem_cells.features", "--alias", "./alias_databases/Mus_musculus.tsv"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Mus_musculus-Embryonic_stem_cells.tsv", "--features", "./input_datasets/Mus_musculus-Embryonic_stem_cells.features", "--alias", "./alias_databases/Mus_musculus.tsv", "--search", "hoxa6"))
# args <- parser_arguments(args = c("--PCHiC", "./input_datasets/Mus_musculus-Embryonic_stem_cells.tsv", "--features", "./input_datasets/Mus_musculus-Embryonic_stem_cells.features", "--alias", "./alias_databases/Mus_musculus.tsv", "--chromosome", "PP"))

args <- parser_arguments(args)

PCHiC <- load_PCHiC(args$PCHiC)

PCHiC <- filter_by_threshold(PCHiC, args$wt_threshold)

PCHiC <- add_PCHiC_types(PCHiC)

PCHiC_ALL <- PCHiC

if (!is.null(args$chromosome)) {
  if (args$chromosome != "PP") {
    PCHiC <- filter_by_chromosome(PCHiC, args$chromosome)
  } else {
    PCHiC <- PCHiC[PCHiC$type == "P-P", ]
  }
}
# No promoters promoters interaction network
if (nrow(PCHiC) == 0) {
  cat("{}")

  quit(status = 0)
}

curated_PCHiC_vertex <- generate_vertex(PCHiC)

if (!is.null(args$bait_names)) {
  suppressPackageStartupMessages(library(GenomicRanges))
  curated_PCHiC_vertex <- generate_real_bait_names(curated_PCHiC_vertex, args$bait_names)
}

if (!is.null(args$alias)) {
  suppressPackageStartupMessages(library(GenomicRanges))
  organism <- str_split(basename(args$PCHiC), fixed("-"))[[1]][1]
  alias <- read_tsv(args$alias, col_types = cols(chr = col_character()))
  if (organism == "Mus_musculus") {
    curated_PCHiC_vertex <- generate_alias_mus(curated_PCHiC_vertex, alias)
  } else if (organism == "Homo_sapiens") {
    curated_PCHiC_vertex <- generate_alias_homo(curated_PCHiC_vertex, alias)
  }
}

if (!is.null(args$intronic_regions)) {
  suppressPackageStartupMessages(library(GenomicRanges))
  curated_PCHiC_vertex <- generate_intronics_regions(curated_PCHiC_vertex, args$intronic_regions)
}

# Finally add all features to their corresponding fragments
initial_features <- NULL
initial_features_position <- NULL
if (!is.null(args$features)) {
  features <- suppressMessages(read_tsv(file = args$features))
  # Remove chr prefix from the fragment column
  features$fragment <- str_sub(features$fragment, start = 4)
  # Binarize all the features
  # initial_features_not_binarized <- features
  # if ("V2" %in% colnames(features)) {
  #   features["V2"] <- ifelse(features["V2"] <= 0.5, 0, 1)
  # }
  # features[, -1] <- ifelse(features[, -1] == 0.0, 0, 1)
  initial_features <- features
  curated_PCHiC_vertex <- merge_features(curated_PCHiC_vertex, initial_features)
  initial_features_position <- which(colnames(curated_PCHiC_vertex) == colnames(initial_features)[2])
}

curated_PCHiC_edges <- generate_edges(PCHiC)

## ------------------------------------------------------------------------
# Generate the network
# Use fragment id as igraph id
curated_PCHiC_vertex$name <- curated_PCHiC_vertex$fragment
net <- graph_from_data_frame(curated_PCHiC_edges, directed = F, curated_PCHiC_vertex)
# Remove igraph id after the network generation
curated_PCHiC_vertex$name <- NULL
# Remove repeated edges and self loops
# When the repeated edges are being removing take the first type (always is c("P-P, "P-P"))
net <- simplify(net, edge.attr.comb = "first")

# Add additional network metadata
V(net)$total_degree <- degree(net)

# Search the required subnetwork
if (!is.null(args$search)) {
  if (str_detect(args$search, "(([12]?[0-9])|([XYxy])):\\d+(-\\d+)?$")) {
    # Only load GenomicRanges if the string searched is a range
    suppressPackageStartupMessages(library(GenomicRanges))
  }

  required_subnet <- search_subnetwork(args$search, args$expand, args$nearest, net, curated_PCHiC_vertex, ensembl2name)
} else {
  required_subnet <- net
}

# Convert the required subnetwork to Cytoscape Json format
if (is.null(required_subnet)) {
  cat("{}")
} else {
  library(rjson)

  if (!is.null(args$pipeline)) {

    # Prepare folder structure
    output_folder <- args$pipeline
    # First extract the filename from the path and then remove the suffix
    filename <- str_remove(basename(args$PCHiC), "\\..+$")
    organism <- str_split(filename, "-")[[1]][1]
    cell_type <- str_split(filename, "-")[[1]][2]

    # Generate graph metadata
    graph_metadata <- generate_graph_metadata(net)

    # Save graph metadata
    if (!is.null(args$chromosome)) {
      write(toJSON(graph_metadata), file = file.path(output_folder, organism, cell_type, "metadata", str_c("chr", args$chromosome, ".json")))
    } else {
      write(toJSON(graph_metadata), file = file.path(output_folder, organism, cell_type, "metadata.json"))
    }

    features <- NULL
    # Generate again all the network but without removing by chromosome
    if (!is.null(args$chromosome) && args$chromosome == "1") {
      # We need to take all the network for statistics insteand of chromosome network
      if (!is.null(args$chromosome)) {
        curated_PCHiC_vertex <- generate_vertex(PCHiC_ALL)

        if (!is.null(args$bait_names)) {
          suppressPackageStartupMessages(library(GenomicRanges))
          curated_PCHiC_vertex <- generate_real_bait_names(curated_PCHiC_vertex, args$bait_names)
        }

        ensembl2name <- NULL
        if (!is.null(args$alias)) {
          suppressPackageStartupMessages(library(GenomicRanges))
          alias <- read_tsv(args$alias, col_types = cols(chr = col_character()))
          organism <- str_split(basename(args$PCHiC), fixed("-"))[[1]][1]
          if (organism == "Mus_musculus") {
            curated_PCHiC_vertex <- generate_alias_mus(curated_PCHiC_vertex, alias)
          } else if (organism == "Homo_sapiens") {
            curated_PCHiC_vertex <- generate_alias_homo(curated_PCHiC_vertex, alias)
          }
          ensembl2name <- alias$`Gene name`
          names(ensembl2name) <- alias$`Ensembl gene ID`
        }

        if (!is.null(args$intronic_regions)) {
          suppressPackageStartupMessages(library(GenomicRanges))
          curated_PCHiC_vertex <- generate_intronics_regions(curated_PCHiC_vertex, args$intronic_regions)
        }

        if (!is.null(args$features)) {
          curated_PCHiC_vertex <- merge_features(curated_PCHiC_vertex, initial_features)
          # Generate features
          features <- sort(colnames(curated_PCHiC_vertex[initial_features_position:length(curated_PCHiC_vertex)]))
        } else {
          features <- list()
        }
        curated_PCHiC_edges <- generate_edges(PCHiC_ALL)
        net <-
          graph_from_data_frame(curated_PCHiC_edges, directed = F, curated_PCHiC_vertex)
        V(net)$total_degree <- degree(net)
      }
      # Only generate network metadata for the first chromosome because is the same for all
      # Generate chromosomes
      chromosomes <- unique(curated_PCHiC_vertex$chr)
      # Remove MT mouse chromosome
      # Add promoter-promoter only networks
      chromosomes <- c(chromosomes, "PP")
      chromosomes <- str_sort(chromosomes[chromosomes != "MT"], numeric = T)

      # Generate suggestions
      suggestions <- generate_suggestions(net)

      # Generate gchas
      chaser_input_PCHiC <- generate_input_chaser_PCHiC(PCHiC_ALL)
      chaser_net <- make_chromnet(chaser_input_PCHiC)

      if (!is.null(args$features)) {
        chaser_input_features <- generate_input_chaser_features(curated_PCHiC_vertex, initial_features_position)
        chaser_net <- chaser::load_features(chaser_net, chaser_input_features, type = "features_on_nodes", missingv = 0)

        # All network
        baits <- unique(chaser::export(chaser_net, "originaldf")$node_from)
        net_features_metadata <- generate_features_metadata(chaser_net, randomize = 100, preserve.nodes = baits)
        # PP network only
        chaser_net_bb <- chaser::subset_chromnet(chaser_net, method = "nodes", nodes1 = baits)
        pp_net_features_metadata <- generate_features_metadata(chaser_net_bb)
        # PO network only
        all_oes <- unique(chaser::export(chaser_net, "originaldf")$node_to)
        oes <- all_oes[!(all_oes %in% baits)]
        chaser_net_bo <- chaser::subset_chromnet(chaser_net, method = "nodes", nodes1 = baits, nodes2 = oes)
        po_net_features_metadata <- generate_features_metadata(chaser_net_bo)
        features_metadata <- list(net = net_features_metadata, pp = pp_net_features_metadata, po = po_net_features_metadata)
        write(toJSON(features_metadata), file = file.path(output_folder, organism, cell_type, "features_metadata.json"))
        curated_PCHiC_vertex[, initial_features_position:length(curated_PCHiC_vertex)] <- round(curated_PCHiC_vertex[, initial_features_position:length(curated_PCHiC_vertex)], 2)
      }

      # Save all metadata to their folders
      # Save suggestions
      write(toJSON(suggestions), file = file.path(output_folder, organism, cell_type, "suggestions.json"))
      # Save chromosomes
      write(toJSON(chromosomes), file = file.path(output_folder, organism, "chromosomes.json"))
      # Save features
      write(toJSON(features), file = file.path(output_folder, organism, cell_type, "features.json"))
      # Save search cache
      save(net, curated_PCHiC_vertex, ensembl2name, file = file.path(output_folder, organism, cell_type, "search_cache.Rdata"), compress = F)
      # Save features generation cache
      save(chaser_net, file = file.path(output_folder, organism, cell_type, "merge_features_cache.Rdata"), compress = F)
    }
  }
  # Convert the required subnetwork to Cytoscape Json format
  required_subnet_json <- generate_cytoscape_json(required_subnet)
  cat(required_subnet_json)
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
