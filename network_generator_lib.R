parser_arguments <- function(args) {
  parser <- OptionParser(description = "Separated values file to cytoscape json mapper")
  parser <- add_option(parser, "--PCHiC", help = "Separated values file PCHiC as input file")
  parser <- add_option(parser, "--wt_threshold",
    type = "double",
    default = 5.0,
    help = "The minimun value for considering the edge [default: %default]"
  )
  parser <- add_option(parser, "--features",
    help = "Separated values file of features as input file"
  )
  parser <- add_option(parser, "--search",
    help = "Search node by name or fragment position in the graph to generate a neighborhood subgraph"
  )
  parser <- add_option(parser, "--chromosome",
    help = "Filter by chromosome"
  )
  parser <- add_option(parser, "--no-features-binarization",
    action = "store_true",
    default = F,
    help = "Features will be binarized by default"
  )
  parser <- add_option(parser, "--nearest",
    action = "store_true",
    default = F,
    help = "Search the nearest range"
  )
  parser <- add_option(parser, "--expand",
    type = "integer",
    default = 0,
    help = "Number of bases to expand the search by range"
  )
  parser <- add_option(parser, "--pipeline",
    metavar = "folder",
    help = "Run the pipeline mode:
                       \t\tIt generates the folder structure datasets/Organism/Cell_type\n
                       \t\tIt takes the PCHiC file Organism_Cell_type.PCHiC and the features file Organism_Cell_type.features if exists\n
                       \t\tThe metadata files: datasets/Organism/Cell_type/{search.Rdata, suggestions.json, features.json}\n
                       \t\tThe chromosomes: datasets/Organism/Cell_type/chromosomes/chrNN.json (according to the organisms Homo_sapiens and Mus_musculus)"
  )
  parser <- add_option(parser, "--organism",
    help = "Select an organism, only for searcher_query.R"
  )
  parser <- add_option(parser, "--cell_type",
    help = "Select a cell_type, only for searcher_query.R"
  )

  parser <- add_option(parser, "--only_pp_interactions",
    action = "store_true",
    default = F,
    help = "Use only promoter-promoter interactions from the network"
  )
  return(parse_args(parser, args, convert_hyphens_to_underscores = T))
}
## ------------------------------------------------------------------------
search_vertex_by_name <-
  function(vertex, net) {
    # Detect if we are searching by position (we are working with mouse chromosomes by now) or by name
    # Always return NULL if it doesn't exist the vertex in the graph
    if (str_detect(vertex, "^(([12]?[0-9])|([XYxy]))_\\d+$")) {
      # Always use upper case here
      vertex <- str_to_upper(vertex)
      if (!vertex %in% V(net)$name) {
        return(NULL)
      }
      required_vertex <- V(net)[vertex]
      required_subnet <-
        make_ego_graph(net, nodes = required_vertex)[[1]]
      required_subnet <- set_vertex_attr(required_subnet, "searched", value = "false")
      required_subnet <- set_vertex_attr(required_subnet, "searched", index = V(net)[required_vertex]$name, value = "true")
      return(required_subnet)
    } else {
      # Always search in lowercase
      vertex <- str_to_lower(vertex)

      searched_vertex_index <- which(str_detect(V(net)$gene_names, fixed(vertex)))

      if (length(searched_vertex_index) == 0) {
        return(NULL)
      }

      required_vertex <- V(net)[searched_vertex_index]

      # Multiple fragments here
      required_subnet <- make_ego_graph(net, nodes = required_vertex)

      required_union_subnet <- union_graphs_with_attributes(required_subnet)

      required_subnet <- set_vertex_attr(required_union_subnet, "searched", value = "false")
      required_subnet <- set_vertex_attr(required_subnet, "searched", index = V(net)[searched_vertex_index]$name, value = "true")
      return(required_subnet)
    }
  }

nearest_subnetwork <- function(required_range, net, curated_chrs_vertex_ranges) {
  nearest_range_index <-
    nearest(required_range, curated_chrs_vertex_ranges)
  required_vertex <-
    curated_chrs_vertex_ranges[nearest_range_index]$fragment
  if (is.null(required_vertex)) {
    required_subnet <- NULL
  } else {
      # Multiple fragments here
      required_subnet <- make_ego_graph(net, nodes = required_vertex)

      required_union_subnet <- union_graphs_with_attributes(required_subnet)

      required_subnet <- set_vertex_attr(required_union_subnet, "searched", value = "false")
      required_subnet <- set_vertex_attr(required_subnet, "searched", index = V(net)[searched_vertex_index]$name, value = "true")
  }
  return(required_subnet)
}

search_vertex_by_range <- function(search, expand, nearest, net, curated_chrs_vertex) {
  curated_chrs_vertex_ranges <-
    makeGRangesFromDataFrame(
      curated_chrs_vertex,
      keep.extra.columns = T,
      ignore.strand = FALSE
    )
  required_range <- GRanges(search)
  # Expand the selected range if it is required
  if (expand != 0) {
    start(required_range) <- start(required_range) - expand
    end(required_range) <- end(required_range) + expand
  }
  # Work with the nearest if it is required
  if (nearest) {
    required_subnet <- nearest_subnetwork(required_range, net, curated_chrs_vertex_ranges)
  } else {
    # Work with overlaps instead
    overlaps_index <-
      subjectHits(findOverlaps(required_range, curated_chrs_vertex_ranges))
    required_vertex <-
      curated_chrs_vertex_ranges[overlaps_index]$fragment
    required_vertex_with_neighbours <-
      names(unlist(lapply(required_vertex, function(rv) {
        neighbors(net, rv)
      })))
    # Add the overlapping vertex
    required_vertex_with_neighbours <-
      unique(c(required_vertex_with_neighbours, required_vertex))
    if (length(required_vertex_with_neighbours) == 0) {
      required_subnet <- nearest_subnetwork(required_range, net, curated_chrs_vertex_ranges)
    } else {
      required_subnet <-
        induced_subgraph(net, vids = required_vertex_with_neighbours)
      required_subnet <- set_vertex_attr(required_subnet, "searched", value = "false")
      required_subnet <- set_vertex_attr(required_subnet, "searched", index = required_vertex, value = "true")
    }
  }
  return(required_subnet)
}
## Generate required subnetwork
search_subnetwork <- function(search, expand, nearest, net, curated_chrs_vertex) {
  if (!is.null(search)) {
    if (str_detect(search, "(([12]?[0-9])|([XYxy])):\\d+(-\\d+)?$")) {
      # We are working with a range
      required_subnet <-
        search_vertex_by_range(search, expand, nearest, net, curated_chrs_vertex)
    } else {
      required_subnet <- search_vertex_by_name(search, net)
    }
    if (!is.null(required_subnet)) {
      # Always recalculate degrees for each neighborhood
      V(required_subnet)$degree <- degree(required_subnet)
    }
  } else {
    required_subnet <- net
  }
  return(required_subnet)
}

## Generate Cytoscape JSON
generate_cytoscape_json <- function(required_subnet) {
  if (is.null(required_subnet)) {
    return("{}")
  }
  # Recover vertices neighbourhood subnetwork dataframe from the graph
  vertices_df <- as_tibble(igraph::as_data_frame(required_subnet, what = "vertices"))
  # _ in column names is not valid in Cytoscape JSON
  vars <- c(id = "name", names = "gene_names")
  vertices_df <- dplyr::rename(vertices_df, !!vars)
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
  colnames(edges_df)[1] <- "source"
  colnames(edges_df)[2] <- "target"
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
  return(toJSON(JSON_df, indent = 2))
}

# Load PCHiC
load_PCHiC <- function(PCHiC_file) {
  suppressMessages(read_tsv(
    file = PCHiC_file,
    col_types = cols(baitChr = col_character(), oeChr = col_character())
  ))
}

# Filter by threshold (column 12)
filter_by_threshold <- function(PCHiC, threshold) {
  PCHiC[PCHiC[12] > threshold, ]
}

# Filter by chromosome
filter_by_chromosome <- function(PCHiC, chromosome) {
  PCHiC[which(PCHiC$baitChr == chromosome
  | PCHiC$oeChr == chromosome), ]
}

# Generate vertex from a PCHiC file
generate_vertex <- function(PCHiC) {
  ## ------------------------------------------------------------------------
  # Join chr number with the start position
  # Also join bait and oe in the same column
  fragment <- c(
    paste(PCHiC$baitChr, PCHiC$baitStart, sep = "_"),
    paste(PCHiC$oeChr, PCHiC$oeStart, sep = "_")
  )
  # Extract bait and oe names and join them to the same column
  gene_names <- c(PCHiC$baitName, PCHiC$oeName)
  # Remove repeted nodes
  chr <- c(PCHiC$baitChr, PCHiC$oeChr)
  start <- c(PCHiC$baitStart, PCHiC$oeStart)
  end <- c(PCHiC$baitEnd, PCHiC$oeEnd)
  # Type
  bait_types <- sapply(PCHiC$type, function(type) {
    str_split(type, fixed("-"))[[1]][1]
  })
  oe_types <- sapply(PCHiC$type, function(type) {
    str_split(type, fixed("-"))[[1]][2]
  })
  type <- c(bait_types, oe_types)
  # Only use lowercase in gene names
  gene_names <- str_to_lower(gene_names)
  # Uniform "." and NA name to empty strings
  gene_names <-
    ifelse(gene_names == "." | is.na(gene_names), "", gene_names)
  # Remove duplicated vertex
  curated_PCHiC_vertex <-
    distinct(tibble(fragment, gene_names, chr, start, end, type))
  # Replace separators by one space
  curated_PCHiC_vertex$gene_names <-
    str_replace_all(curated_PCHiC_vertex$gene_names, "[,;]", " ")
  curated_PCHiC_vertex
}

# Load the features
generate_features <- function(curated_PCHiC_vertex, features_file, binarization = T) {
  features <-
    suppressMessages(read_tsv(file = features_file))
  # Remove chr prefix from the fragment column
  features$fragment <- str_sub(features$fragment, start = 4)
  # Binarize all the features
  if (!args$no_features_binarization && binarization) {
    if ("V2" %in% colnames(features)) {
      features["V2"] <- ifelse(features["V2"] <= 0.5, 0, 1)
    }
    features[, -1] <- ifelse(features[, -1] == 0.0, 0, 1)
  }

  left_join(curated_PCHiC_vertex, features, by = "fragment")
}


# Generate a dataframe with the extremes of the edges
generate_edges <- function(PCHiC) {
  baits <- paste(PCHiC$baitChr, PCHiC$baitStart, sep = "_")
  oes <- paste(PCHiC$oeChr, PCHiC$oeStart, sep = "_")
  curated_PCHiC_edges <- tibble(source = baits, target = oes, type = PCHiC$type)
  curated_PCHiC_edges
}

# Generate gene name list
generate_suggestions <- function(net) {
  suggestions <- sort(unique(unlist(sapply(V(net)$gene_names, function(gene_names) {str_split(gene_names, fixed(" "))}))))
  if (suggestions[1] == "") {
    suggestions <- suggestions[-1]
  }
  suggestions
}

generate_graph_metadata <- function(net) {
  nodes <- length(V(net))
  degree_average <- round(mean(degree(net)), 2)
  edges <- length(E(net))
  connected_components <- components(net)$no
  largest_connected_component <- sort(components(net)$csize, decreasing = T)[1]
  nodes_in_largest_connected_component <- paste0(round(largest_connected_component / nodes * 100, 2), "%")
  network_diameter <- diameter(net)
  promoters <- sum(V(net)$type == "P")
  other_ends <- sum(V(net)$type == "O")
  PP_edges <- sum(E(net)$type == "P-P")
  PO_edges <- sum(E(net)$type == "P-O")
  edge_ends <- ends(net, E(net))
  edge_ends_source <- edge_ends[, 1]
  edge_ends_target <- edge_ends[, 2]
  edge_ends_source_chromosome <- unname(sapply(edge_ends_source, function(edge_end_source) {
    str_split(edge_end_source, "_")[[1]][1]
  }))
  edge_ends_target_chromosome <- unname(sapply(edge_ends_target, function(edge_end_target) {
    str_split(edge_end_target, "_")[[1]][1]
  }))
  interchromosomal_interactions <- sum(edge_ends_source_chromosome != edge_ends_target_chromosome)
  clustering_coefficient <- round(transitivity(net), 2)

  network_properties <- list(
    Nodes = nodes, Edges = edges, Promoters = promoters,
    "Other ends" = other_ends, "PP edges" = PP_edges, "PO edges" = PO_edges,
    "Interchromosomal interactions" = interchromosomal_interactions
  )

  network_statistics <- list(
    "Degree average" = degree_average, "Connected components (CC)" = connected_components,
    "Nodes in largest CC" = nodes_in_largest_connected_component,
    "Network diameter" = network_diameter, "Clustering coefficient" = clustering_coefficient
  )

  graph_metadata <- list(network_properties = network_properties, network_statistics = network_statistics)

  graph_metadata
}

add_PCHiC_types <- function(PCHiC) {
  # Add the type for bait and oes
  # Be careful because in the oe column there are many baits
  # So oe is only oe if it not exist in the bait column
  baits <- paste(PCHiC$baitChr, PCHiC$baitStart, PCHiC$baitEnd, sep = "_")
  oes <- paste(PCHiC$oeChr, PCHiC$oeStart, PCHiC$oeEnd, sep = "_")
  PCHiC$type <- ifelse(oes %in% baits, "P-P", "P-O")
  PCHiC
}

generate_input_data_gchas <- function(PCHiC, curated_PCHiC_vertex, randomize = F) {
  chaser_PCHiC <- PCHiC[, c(1:3, 6:8)]
  if (any(grepl("MT", chaser_PCHiC$baitChr))) {
    chaser_PCHiC <- chaser_PCHiC[-grep("MT", chaser_PCHiC$baitChr), ]
  }
  chaser_PCHiC$baitChr <- paste0("chr", chaser_PCHiC$baitChr)
  chaser_PCHiC$oeChr <- paste0("chr", chaser_PCHiC$oeChr)
  chaser_PCHiC_df <- as.data.frame(chaser_PCHiC)
  chaser_features <- curated_PCHiC_vertex
  chaser_features$fragment <- paste(curated_PCHiC_vertex$chr, paste(curated_PCHiC_vertex$start, curated_PCHiC_vertex$end, sep = "-"), sep = ":")
  chaser_features$fragment <- paste0("chr", chaser_features$fragment)
  chaser_features <- select(chaser_features, c(1, 7:length(curated_PCHiC_vertex)))
  chaser_features_df <- as.data.frame(chaser_features)
  rownames(chaser_features_df) <- chaser_features_df[, 1]
  chaser_features_df[, 1] <- NULL

  list(chaser_PCHiC = chaser_PCHiC_df, chaser_features = chaser_features_df)
}

generate_gchas <- function(chaser_net, features, randomize = F) {
  if (randomize) {
    chaser_net <- chaser::randomize(chaser_net)
  }

  chas <- sapply(features, function(feature) {
    round(gchas(chaser_net, feature), 2)
  })

  chas
}

generate_features_metadata <- function(PCHiC) {
  curated_PCHiC_vertex <- generate_vertex(PCHiC)
  # Always without binarization to calcule the gchas number
  curated_PCHiC_vertex_not_binarized <- generate_features(curated_PCHiC_vertex, args$features, binarization = F)
  curated_PCHiC_vertex_binarized <- generate_features(curated_PCHiC_vertex, args$features, binarization = T)

  gchas_input <- generate_input_data_gchas(PCHiC, curated_PCHiC_vertex_not_binarized)
  chaser_net <- chaser::chromnet_of_data_frames(gchas_input$chaser_PCHiC, gchas_input$chaser_features)

  features <- sort(colnames(curated_PCHiC_vertex_binarized[7:length(curated_PCHiC_vertex_binarized)]))

  chas <- generate_gchas(chaser_net, features)

  curated_PCHiC_edges <- generate_edges(PCHiC)
  net_not_binarized <- graph_from_data_frame(curated_PCHiC_edges, directed = F, curated_PCHiC_vertex_not_binarized)
  net_binarized <- graph_from_data_frame(curated_PCHiC_edges, directed = F, curated_PCHiC_vertex_binarized)

  # Calculate random ChAs
  random_chas_list <- list()
  for (i in 1:100) {
    random_chas_list[[i]] <- generate_gchas(chaser_net, features, randomize = T)
  }

  random_chas_min <- c()
  random_chas_max <- c()
  for (feature in features) {
    random_chas_feature <- sapply(1:100, function(i) {
      random_chas_list[[i]][feature]
    })
    random_chas_min <- c(random_chas_min, min(random_chas_feature))
    random_chas_max <- c(random_chas_max, max(random_chas_feature))
  }

  random_chas <- paste(random_chas_min, random_chas_max, sep = ",")
  names(random_chas) <- features

  # mean degree of nodes with one specific feature
  mean_degree <- sapply(features, function(feature) {
    round(mean(degree(net_binarized)[vertex_attr(net_binarized)[[feature]] != 0], na.rm = T), 2)
  })
  # abundance of each feature
  abundance <- sapply(features, function(feature) {
    round(mean(vertex_attr(net_not_binarized)[[feature]], na.rm = T), 2)
  })
  list("Abundance" = abundance, "ChAs" = chas, "Random ChAs interval" = random_chas, "Mean degree" = mean_degree)
}

# Adapted to N graphs from only 2, see https://stackoverflow.com/a/46338136
union_graphs_with_attributes <- function(graph_list) {

  # Internal function that cleans the names of a given attribute and merge them
  merge_attributes <- function(union_graph, component) {
    # get component names
    gNames <- parse(text = (paste0(component, "_attr_names(union_graph)"))) %>% eval()
    # find names that have a "_1", "_2" ... "_N" at the end
    AttrNeedsCleaning <- grepl("(_\\d)$", gNames)
    # Suffix number list to find the max to the loop counter
    suffix_number_list <- str_extract(gNames, "\\d+$")
    if (length(suffix_number_list) == 0) {
      return(union_graph)
    }
    max_suffix_number <- max(as.numeric(suffix_number_list), na.rm = T)
    # remove the _N ending
    StemName <- gsub("(_\\d)$", "", gNames)

    NewnNames <- unique(StemName[AttrNeedsCleaning])
    # replace attribute name for all attributes
    for (i in NewnNames) {
      attr_list <- list()
      # Save in a list all attribute values
      for (j in 1:max_suffix_number) {
        attr_list <- append(attr_list, list(parse(text = (paste0(component, "_attr(union_graph,'", paste0(i, "_", j), "')"))) %>% eval()))
        union_graph <- parse(text = (paste0("delete_", component, "_attr(union_graph,'", paste0(i, "_", j), "')"))) %>% eval()
      }

      # Collapse values replacing the NA with corresponding existing value from the different graphs
      values <- do.call(pmin, c(attr_list, na.rm = T))

      union_graph <- parse(text = (paste0("set_", component, "_attr(union_graph, i, value = values)"))) %>% eval()

    }

    return(union_graph)
  }


  union_graph <- do.call(igraph::union, graph_list)
  # loop through each attribute type in the graph and clean
  for (component in c("graph", "edge", "vertex")) {
    union_graph <- merge_attributes(union_graph, component)
  }

  return(union_graph)
}
