#' @title parser_arguments
#' @description Convert input parameters to a named list
#' @param args input parameters
#' @return named list with the input parameters
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
  parser <- add_option(parser, "--bait_names",
    help = "Separated values file of bait names as input file"
  )
  parser <- add_option(parser, "--alias",
    help = "Separated values file of alias as input file"
  )
  parser <- add_option(parser, "--intronic_regions",
    help = "Separated values file of intronics as input file"
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

#' @title search_vertex_by_names
#' @description search vertex by names
#' @param vertex one or many vertex
#' @param net igraph network
#' @param ensembl2name named list with ensembl id to name translation
#' @return specific neighborhood for that set of vertex
search_vertex_by_names <- function(vertex, net, ensembl2name) {
  # Detect if we are searching by position (we are working with mouse chromosomes by now) or by name
  # Always return NULL if it doesn't exist the vertex in the graph
  if (str_detect(vertex, "^(([12]?[0-9])|([XYxy]))_\\d+_\\d+$")) {
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
    curated_vertex_list <- c()
    for (v in vertex) {
      if (str_detect(v, "^ens(mus)?g\\d+")) {
        curated_vertex_list <- c(curated_vertex_list, str_to_lower(ensembl2name[str_to_upper(v)]))
        if (is.na(v)) {
          return(NULL)
        }
      } else {
        curated_vertex_list <- c(curated_vertex_list, v)
      }
    }

    all_gene_names_together <- NULL
    if ("alias" %in% vertex_attr_names(net)) {
      all_gene_names_together <- str_c(V(net)$gene_names, V(net)$alias, sep = " ")
    } else {
      all_gene_names_together <- V(net)$gene_names
    }
    searched_vertex_index <- unique(unlist(sapply(curated_vertex_list, function(vertex) {
      str_which(str_to_lower(all_gene_names_together), regex(str_c("\\b", vertex, "\\b")))
    })))

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

#' @title nearest_subnetwork
#' @description Find the nearest subnetwork closer to the required range
#' @param required_range Required range reference
#' @param net igraph network
#' @param curated_PCHiC_vertex_ranges curated_PCHiC_vertex with ranges
#' @return the neighborhood closer to the required range
nearest_subnetwork <- function(required_range, net, curated_PCHiC_vertex_ranges) {
  nearest_range_index <- nearest(required_range, curated_PCHiC_vertex_ranges)
  searched_vertex_index <- curated_PCHiC_vertex_ranges[nearest_range_index]$fragment
  if (is.null(searched_vertex_index)) {
    required_subnet <- NULL
  } else {
    # Multiple fragments here
    required_subnet <- make_ego_graph(net, nodes = searched_vertex_index)

    required_union_subnet <- union_graphs_with_attributes(required_subnet)

    required_subnet <- set_vertex_attr(required_union_subnet, "searched", value = "false")
    required_subnet <- set_vertex_attr(required_subnet, "searched", index = V(net)[searched_vertex_index]$name, value = "true")
  }
  return(required_subnet)
}

#' @title search_vertex_by_range
#' @description Search neighborhood by range
#' @param search the range to be used
#' @param expand expand the range
#' @param nearest search the nearest neighborhood always used if there is not neighborhood
#' @param net igraph network
#' @param curated_PCHiC_vertex  curated vertex from PCHiC
#' @return the found neighborhood
search_vertex_by_range <- function(search, expand, nearest, net, curated_PCHiC_vertex) {
  curated_PCHiC_vertex_ranges <-
    makeGRangesFromDataFrame(
      curated_PCHiC_vertex,
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
    required_subnet <- nearest_subnetwork(required_range, net, curated_PCHiC_vertex_ranges)
  } else {
    # Work with overlaps instead
    overlaps_index <- subjectHits(findOverlaps(required_range, curated_PCHiC_vertex_ranges))
    required_vertex <- curated_PCHiC_vertex_ranges[overlaps_index]$fragment
    required_vertex_with_neighbours <- required_vertex
    # required_vertex_with_neighbours <-
    #   names(unlist(lapply(required_vertex, function(rv) {
    #     neighbors(net, rv)
    #   })))
    # # Add the overlapping vertex
    # required_vertex_with_neighbours <- unique(c(required_vertex_with_neighbours, required_vertex))
    if (length(required_vertex_with_neighbours) == 0) {
      required_subnet <- nearest_subnetwork(required_range, net, curated_PCHiC_vertex_ranges)
      required_subnet <- set_vertex_attr(required_subnet, "searched", value = "false")
      required_subnet <- set_vertex_attr(required_subnet, "searched", index = required_vertex, value = "true")
    } else {
      required_subnet <- induced_subgraph(net, vids = required_vertex_with_neighbours)
      required_subnet <- set_vertex_attr(required_subnet, "searched", value = "false")
      required_subnet <- set_vertex_attr(required_subnet, "searched", index = required_vertex, value = "true")
    }
  }
  return(required_subnet)
}

#' @title search_subnetwork
#' @description Search the subnetwork for the required string
#' @param search search string
#' @param expand expand the range search
#' @param nearest find the nearest range
#' @param net igraph network
#' @param curated_PCHiC_vertex curated vertex from PCHiC
#' @param ensembl2name conversion between ensembl id to gene name
#' @return the found neighborhood
search_subnetwork <- function(search, expand, nearest, net, curated_PCHiC_vertex, ensembl2name) {
  if (!is.null(search)) {
    if (str_detect(search, "(([12]?[0-9])|([XYxy])):\\d+(-\\d+)?$")) {
      # We are working with a range
      required_subnet <-
        search_vertex_by_range(search, expand, nearest, net, curated_PCHiC_vertex)
    } else if (str_detect(search, "(\\w+,\\w+)+")) {
      gene_list_to_search <- str_split(search, "[, \\t]")[[1]]
      required_subnet <- search_vertex_by_names(gene_list_to_search, net, ensembl2name)
    } else {
      required_subnet <- search_vertex_by_names(search, net, ensembl2name)
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

#' @title generate_cytoscape_json
#' @description Convert the igraph network to Cytoscape json format
#' @param required_subnet the network to be converted
#' @return the json string
#' @seealso
#'  \code{\link[igraph]{as_data_frame}}
#'  \code{\link[dplyr]{select}}
#' @importFrom igraph as_data_frame
#' @importFrom dplyr rename
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
  edges_df$id <- str_c(edges_df$source, edges_df$target, sep = "~")
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

#' @title load_PCHiC
#' @description Load PCHiC data
#' @param PCHiC_file PCHiC file
#' @return PCHiC as tibble
load_PCHiC <- function(PCHiC_file) {
  suppressMessages(read_tsv(
    file = PCHiC_file,
    col_types = cols(baitChr = col_character(), oeChr = col_character())
  ))
}

#' @title filter_by_threshold
#' @description Filter PCHiC data using the threshold
#' @param PCHiC PCHiC data
#' @param threshold Requested threshold (always column 12)
#' @return Filtered PCHiC
filter_by_threshold <- function(PCHiC, threshold) {
  PCHiC[PCHiC[12] > threshold, ]
}

#' @title filter_by_chromosome
#' @description filter PCHiC data by chromosome
#' @param PCHiC PCHiC data
#' @param chromosome Chromosome used as filter
#' @return PCHiC filtered which contains interchromosomal interactions
filter_by_chromosome <- function(PCHiC, chromosome) {
  PCHiC[which(PCHiC$baitChr == chromosome | PCHiC$oeChr == chromosome), ]
}

#' @title generate_vertex
#' @description Generate vertex data from PCHiC data
#' @param PCHiC data
#' @return generated vertex
generate_vertex <- function(PCHiC) {
  ## ------------------------------------------------------------------------
  # Join chr number with the start position
  # Also join bait and oe in the same column
  fragment <- c(
    str_c(PCHiC$baitChr, PCHiC$baitStart, PCHiC$baitEnd, sep = "_"),
    str_c(PCHiC$oeChr, PCHiC$oeStart, PCHiC$oeEnd, sep = "_")
  )
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
  # Extract bait and oe names and join them to the same column
  gene_names <- c(PCHiC$baitName, PCHiC$oeName)
  # Remove duplicated vertex
  curated_PCHiC_vertex <- distinct(tibble(fragment, gene_names, chr, start, end, type))
  # Uniform "." and NA name to empty strings
  curated_PCHiC_vertex$gene_names <- ifelse(curated_PCHiC_vertex$gene_names == "." | is.na(curated_PCHiC_vertex$gene_names), "", curated_PCHiC_vertex$gene_names)
  # Delete transcript names
  curated_PCHiC_vertex$gene_names <- sapply(curated_PCHiC_vertex$gene_names, function(gene_name) {
    str_trim(str_remove_all(gene_name, "-\\d+\\b"))
  })
  # Collapse unique names
  curated_PCHiC_vertex$gene_names <- sapply(curated_PCHiC_vertex$gene_names, function(gene_name) {
    str_c(unique(str_split(gene_name, "[ ;,]")[[1]]), collapse = " ", sep = " ")
  })
  distinct(curated_PCHiC_vertex)
}

#' @title merge_features
#' @description Merge features to curated vertex data
#' @param curated_PCHiC_vertex curated vertex data
#' @param features data
#' @return curated vertex with features
merge_features <- function(curated_PCHiC_vertex, features) {
  if (str_detect(features$fragment[1], "(([12]?[0-9])|([XYxy]))_\\d+_\\d+$")) {
    curated_PCHiC_vertex <- left_join(curated_PCHiC_vertex, features, by = "fragment")
  } else if (str_detect(features$fragment[1], "(([12]?[0-9])|([XYxy]))_\\d+$")) {
    curated_PCHiC_vertex$fragment_tmp <- str_c(curated_PCHiC_vertex$chr, "_", curated_PCHiC_vertex$start)
    features$fragment_tmp <- features$fragment
    features$fragment <- NULL
    curated_PCHiC_vertex <- left_join(curated_PCHiC_vertex, features, by = "fragment_tmp")
    curated_PCHiC_vertex$fragment_tmp <- NULL
  }
  curated_PCHiC_vertex
}

#' @title generate_edges
#' @description Generate edges from PCHiC data
#' @param PCHiC PCHiC data
#' @return Edges extracted from PCHiC data
generate_edges <- function(PCHiC) {
  baits <- str_c(PCHiC$baitChr, PCHiC$baitStart, PCHiC$baitEnd, sep = "_")
  oes <- str_c(PCHiC$oeChr, PCHiC$oeStart, PCHiC$oeEnd, sep = "_")
  curated_PCHiC_edges <- tibble(source = baits, target = oes, type = PCHiC$type)
  curated_PCHiC_edges
}

#' @title generate_suggestions
#' @description Generate genes names list from all network to be used as suggestions for garden-net frontend
#' @param net igraph network
#' @return Genes names suggestions list
generate_suggestions <- function(net) {
  suggestions <- sort(unique(unlist(sapply(V(net)$gene_names, function(gene_names) {
    str_split(gene_names, fixed(" "))
  }))))
  if (suggestions[1] == "") {
    suggestions <- suggestions[-1]
  }
  suggestions
}

#' @title generate_graph_metadata
#' @description Generate graph metadata from an igraph network
#' @param net igraph network
#' @return named list with many network properties and statistics
generate_graph_metadata <- function(net) {
  nodes <- length(V(net))
  degree_average <- round(mean(degree(net)), 2)
  edges <- length(E(net))
  connected_components <- components(net)$no
  largest_connected_component <- sort(components(net)$csize, decreasing = T)[1]
  nodes_in_largest_connected_component <- str_c(round(largest_connected_component / nodes * 100, 2), "%")
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

#' @title add_PCHiC_types
#' @description Add promoter-promoter or promoter other-end interaction type to PCHiC data
#' @param PCHiC PCHiC data
#' @return PCHiC data with types
add_PCHiC_types <- function(PCHiC) {
  # Add the type for bait and oes
  # Be careful because in the oe column there are many baits
  # So oe is only oe if it not exist in the bait column
  baits <- str_c(PCHiC$baitChr, PCHiC$baitStart, PCHiC$baitEnd, sep = "_")
  oes <- str_c(PCHiC$oeChr, PCHiC$oeStart, PCHiC$oeEnd, sep = "_")
  PCHiC$type <- ifelse(oes %in% baits, "P-P", "P-O")
  PCHiC
}

#' @title generate_input_chaser_PCHiC
#' @description Convert PCHiC to chaser input data
#' @param PCHiC PCHiC data
#' @return PCHiC
generate_input_chaser_PCHiC <- function(PCHiC) {
  chaser_PCHiC <- PCHiC[, c(1:3, 6:8)]
  if (any(grepl("MT", chaser_PCHiC$baitChr))) {
    chaser_PCHiC <- chaser_PCHiC[-grep("MT", chaser_PCHiC$baitChr), ]
  }
  chaser_PCHiC$baitChr <- str_c("chr", chaser_PCHiC$baitChr)
  chaser_PCHiC$oeChr <- str_c("chr", chaser_PCHiC$oeChr)
  chaser_PCHiC_df <- as.data.frame(chaser_PCHiC)
  chaser_PCHiC_df
}

#' @title generate_input_chaser_features
#' @description Extract features from curated vertex to adapt to chaser format
#' @param curated_PCHiC_vertex curated vertex
#' @param initial_features_position first feature position in curated vertex tibble
#' @return features in chaser format
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @importFrom dplyr select
generate_input_chaser_features <- function(curated_PCHiC_vertex, initial_features_position) {
  chaser_features <- curated_PCHiC_vertex
  chaser_features$fragment <- str_c(curated_PCHiC_vertex$chr, str_c(curated_PCHiC_vertex$start, curated_PCHiC_vertex$end, sep = "-"), sep = ":")
  chaser_features$fragment <- str_c("chr", chaser_features$fragment)
  chaser_features <- dplyr::select(chaser_features, c(1, initial_features_position:length(curated_PCHiC_vertex)))
  chaser_features_df <- as.data.frame(chaser_features)
  rownames(chaser_features_df) <- chaser_features_df[, 1]
  chaser_features_df[, 1] <- NULL
  chaser_features_df
}

#' @title generate_features_metadata
#' @description Generate ChAs, random ChAs, abundance and mean degree for each feature
#' @param chaser_net chaser network
#' @param randomize number of random ChAs networks to be used to calculate the random range, Default: 0
#' @param preserve.nodes which nodes to preserve during the randomization process, Default: NULL
#' @return metadata features information
#' @seealso
#'  \code{\link[chaser]{randomize}}
#' @importFrom chaser randomize
generate_features_metadata <- function(chaser_net, randomize = 0, preserve.nodes = NULL) {
  features <- colnames(chaser_net$features)
  chas <- chas(chaser_net)

  if (randomize != 0) {
    # Calculate random ChAs
    random_chaser_net_list <- chaser::randomize(chaser_net, nrandom = randomize, preserve.nodes)
    random_chaser_net_list_chas <- lapply(random_chaser_net_list, function(random_chaser_net) {
      chas(random_chaser_net)
    })

    random_chas_min <- c()
    random_chas_max <- c()
    for (feature_index in 1:length(features)) {
      random_chas_feature <- sapply(1:randomize, function(i) {
        random_chaser_net_list_chas[[i]][feature_index]
      })
      random_chas_min <- c(random_chas_min, min(random_chas_feature))
      random_chas_max <- c(random_chas_max, max(random_chas_feature))
    }

    random_chas <- str_c(round(random_chas_min, 3), round(random_chas_max, 3), sep = ",")
    names(random_chas) <- features
  }

  # mean degree of nodes with one specific feature
  mean_degree <- sapply(features, function(feature) {
    round(mean(chaser_net$degree[chaser_net$features[, feature] != 0], na.rm = T), 2)
  })

  # abundance of each feature
  abundance <- sapply(features, function(feature) {
    round(mean(chaser_net$features[, feature], na.rm = T), 2)
  })

  features_metadata <- list()
  features_metadata[["Abundance"]] <- abundance
  features_metadata[["ChAs"]] <- chas
  if (randomize != 0) {
    features_metadata[["Random ChAs interval"]] <- random_chas
  }
  features_metadata[["Mean degree"]] <- mean_degree
  features_metadata
}

# Adapted to N graphs from only 2, see https://stackoverflow.com/a/46338136
#' @title union_graphs_with_attributes
#' @description Join all attributes from two igraph networks
#' @param graph_list list of igraph to be joined
#' @return igraph network with all attributed merged
#' @seealso
#'  \code{\link[igraph]{union}}
#' @importFrom igraph union
union_graphs_with_attributes <- function(graph_list) {

  # Internal function that cleans the names of a given attribute and merge them
  merge_attributes <- function(union_graph, component) {
    # get component names
    gNames <- parse(text = (str_c(component, "_attr_names(union_graph)"))) %>% eval()
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
        attr_list <- append(attr_list, list(parse(text = (str_c(component, "_attr(union_graph,'", str_c(i, "_", j), "')"))) %>% eval()))
        union_graph <- parse(text = (str_c("delete_", component, "_attr(union_graph,'", str_c(i, "_", j), "')"))) %>% eval()
      }

      # Collapse values replacing the NA with corresponding existing value from the different graphs
      values <- do.call(pmin, c(attr_list, na.rm = T))

      union_graph <- parse(text = (str_c("set_", component, "_attr(union_graph, i, value = values)"))) %>% eval()
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

#' @title generate_alias_homo
#' @description Add homo sapiens alias to the curated vertex
#' @param curated_PCHiC_vertex curated vertex
#' @param alias alias data
#' @return curated vertex with homo sapiens alias
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @importFrom dplyr select
generate_alias_homo <- function(curated_PCHiC_vertex, alias) {
  # First overlap others ends and annotate them
  curated_PCHiC_vertex_with_alias <- NULL
  there_are_other_ends <- any(curated_PCHiC_vertex$type == "O")
  if (there_are_other_ends) {
    alias_grange <- makeGRangesFromDataFrame(alias, keep.extra.columns = T)
    vertex_grange <- makeGRangesFromDataFrame(curated_PCHiC_vertex[curated_PCHiC_vertex$type == "O", ], keep.extra.columns = T)
    merged_overlaps <- mergeByOverlaps(vertex_grange, alias_grange)
    # concatenate fragments with multiple overlaps
    overlaps_tibble <- tibble(range = str_c(seqnames(merged_overlaps$vertex_grange), as.character(ranges(merged_overlaps$vertex_grange)), sep = ":"))
    overlaps_tibble$gene_type <- merged_overlaps$`Gene type`
    overlaps_tibble$ensembl <- merged_overlaps$`Ensembl gene ID`
    overlaps_tibble$name <- merged_overlaps$`Gene name`
    overlaps_tibble$alias <- merged_overlaps$Alias
    overlaps_tibble$hgnc <- merged_overlaps$`HGNC ID`

    # Collapse multiple gene names and annotations in one string
    collapsed_overlaps <- overlaps_tibble %>%
      group_by(range) %>%
      summarise(
        ensembl = str_c(str_replace_na(ensembl), collapse = " "),
        name = str_c(str_replace_na(name), collapse = " "),
        alias = str_c(str_replace_na(alias), collapse = " "),
        hgnc = str_c(str_replace_na(hgnc), collapse = " "),
        gene_type = str_c(str_replace_na(gene_type), collapse = " ")
      )
    # Join the new annotations the original data
    curated_PCHiC_vertex_ranges <- curated_PCHiC_vertex %>%
      mutate(range = str_c(chr, str_c(start, end, sep = "-"), sep = ":"))
    curated_PCHiC_vertex_with_alias <- left_join(curated_PCHiC_vertex_ranges, collapsed_overlaps, by = "range")
    curated_PCHiC_vertex_with_alias$gene_names <- if_else(curated_PCHiC_vertex_with_alias$type == "O", curated_PCHiC_vertex_with_alias$name, curated_PCHiC_vertex_with_alias$gene_names)
    curated_PCHiC_vertex_with_alias <- curated_PCHiC_vertex_with_alias %>% select(-c(name, range))
  } else {
    curated_PCHiC_vertex_with_alias <- curated_PCHiC_vertex
  }
  # Then use promoters transcript names to add the missing annotations
  all_bait_names <- curated_PCHiC_vertex_with_alias %>%
    filter(type == "P") %>%
    select(gene_names) %>%
    mutate(all_bait_names = str_to_upper(str_c(gene_names, sep = " "))) %>%
    pull(all_bait_names)
  curated_bait_names <- sapply(all_bait_names, function(bait_name) {
    str_trim(str_remove_all(bait_name, "-\\d+\\b"))
  })
  curated_bait_names_unique <- sapply(curated_bait_names, function(curated_bait_name) {
    unique(str_split(curated_bait_name, fixed(" "))[[1]])
  })
  names(curated_bait_names_unique) <- NULL
  # curated_bait_names_unique <- sapply(curated_bait_names, function(curated_bait_name){ length(unique(str_split(curated_bait_name, fixed(" "))[[1]]))})
  # names(curated_bait_names_unique) <- NULL
  # Original data with a list of promoters names
  # curated_PCHiC_vertex_with_alias$gene_names <- ifelse(curated_PCHiC_vertex_with_alias$type=="P",curated_bait_names_unique, curated_PCHiC_vertex_with_alias$gene_names)
  # Extract all promoters data to do unnest of the names and then the annotation with the symbols
  curated_bait_names_unique_df <- curated_PCHiC_vertex_with_alias %>% filter(type == "P")
  curated_bait_names_unique_df$all_bait_names <- curated_bait_names_unique
  curated_bait_names_unique_df_unnested <- NULL
  if (there_are_other_ends) {
    curated_bait_names_unique_df_unnested <- curated_bait_names_unique_df %>%
      dplyr::select(-c(hgnc, alias, ensembl, gene_type)) %>%
      unnest()
  } else {
    curated_bait_names_unique_df_unnested <- curated_bait_names_unique_df %>% unnest()
  }
  alias_promoters <- alias %>%
    rename(`Gene name` = "all_bait_names") %>%
    mutate(all_bait_names = str_to_upper(all_bait_names)) %>%
    dplyr::select(-c(chr, start, end))
  promoters_merged_alias <- curated_bait_names_unique_df_unnested %>%
    left_join(alias_promoters, by = "all_bait_names")
  promoters_merged_alias_collapsed <- promoters_merged_alias %>%
    group_by(fragment) %>%
    summarise(
      ensembl = str_c(str_replace_na(`Ensembl gene ID`), collapse = " "),
      gene_type = str_c(str_replace_na(`Gene type`), collapse = " "),
      hgnc = str_c(str_replace_na(`HGNC ID`), collapse = " "),
      alias = str_c(str_replace_na(Alias), collapse = " "),
    )

  # Order promoters according promoters_merged_alias_collapsed
  curated_PCHiC_vertex_with_alias <- curated_PCHiC_vertex_with_alias %>% arrange(fragment)

  curated_PCHiC_vertex_with_alias[curated_PCHiC_vertex_with_alias$type == "P", "alias"] <- promoters_merged_alias_collapsed$alias
  curated_PCHiC_vertex_with_alias[curated_PCHiC_vertex_with_alias$type == "P", "hgnc"] <- promoters_merged_alias_collapsed$hgnc
  curated_PCHiC_vertex_with_alias[curated_PCHiC_vertex_with_alias$type == "P", "ensembl"] <- promoters_merged_alias_collapsed$ensembl
  curated_PCHiC_vertex_with_alias[curated_PCHiC_vertex_with_alias$type == "P", "gene_type"] <- promoters_merged_alias_collapsed$gene_type

  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$alias) | curated_PCHiC_vertex_with_alias$alias == "NA", "alias"] <- ""
  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$hgnc) | curated_PCHiC_vertex_with_alias$hgnc == "NA", "hgnc"] <- ""

  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$ensembl) | curated_PCHiC_vertex_with_alias$ensembl == "NA", "ensembl"] <- ""
  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$gene_type) | curated_PCHiC_vertex_with_alias$gene_type == "NA", "gene_type"] <- ""
  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$gene_names), "gene_names"] <- ""
  curated_PCHiC_vertex_with_alias
}

#' @title generate_alias_mus
#' @description Add mus musculus alias to the curated_vertex
#' @param curated_PCHiC_vertex curated vertex
#' @param alias alias data
#' @return curated vertex alias with mus musculus alias
#' @seealso
#'  \code{\link[dplyr]{select}}
#'  \code{\link[tidyr]{unnest}}
#' @importFrom dplyr select
#' @importFrom tidyr unnest
generate_alias_mus <- function(curated_PCHiC_vertex, alias) {
  # First overlap others ends and annotate them
  curated_PCHiC_vertex_with_alias <- NULL
  there_are_other_ends <- any(curated_PCHiC_vertex$type == "O")
  if (there_are_other_ends) {
    alias_grange <- makeGRangesFromDataFrame(alias, keep.extra.columns = T)
    vertex_grange <- makeGRangesFromDataFrame(curated_PCHiC_vertex[curated_PCHiC_vertex$type == "O", ], keep.extra.columns = T)
    merged_overlaps <- mergeByOverlaps(vertex_grange, alias_grange)
    # concatenate fragments with multiple overlaps
    overlaps_tibble <- tibble(range = str_c(seqnames(merged_overlaps$vertex_grange), as.character(ranges(merged_overlaps$vertex_grange)), sep = ":"))
    overlaps_tibble$ensembl <- merged_overlaps$`Ensembl gene ID`
    overlaps_tibble$name <- merged_overlaps$`Gene name`
    overlaps_tibble$gene_type <- merged_overlaps$`Gene type`
    overlaps_tibble$mgi <- merged_overlaps$`MGI ID`

    # Collapse multile gene names and annotations in one string
    collapsed_overlaps <- overlaps_tibble %>%
      group_by(range) %>%
      summarise(
        ensembl = str_c(str_replace_na(ensembl), collapse = " "),
        name = str_c(str_replace_na(name), collapse = " "),
        gene_type = str_c(str_replace_na(gene_type), collapse = " "),
        mgi = str_c(str_replace_na(mgi), collapse = " ")
      )
    # Join the new annotations the original data
    curated_PCHiC_vertex_ranges <- curated_PCHiC_vertex %>%
      mutate(range = str_c(chr, str_c(start, end, sep = "-"), sep = ":"))
    curated_PCHiC_vertex_with_alias <- left_join(curated_PCHiC_vertex_ranges, collapsed_overlaps, by = "range")
    curated_PCHiC_vertex_with_alias$gene_names <- if_else(curated_PCHiC_vertex_with_alias$type == "O", str_to_upper(curated_PCHiC_vertex_with_alias$name), curated_PCHiC_vertex_with_alias$gene_names)

    curated_PCHiC_vertex_with_alias$mgi <- str_remove_all(curated_PCHiC_vertex_with_alias$mgi, fixed("MGI:"))
    curated_PCHiC_vertex_with_alias <- curated_PCHiC_vertex_with_alias %>% select(-c(name, range))
  } else {
    curated_PCHiC_vertex_with_alias <- curated_PCHiC_vertex
  }
  # Then use promoters transcript names to add the missing annotations
  all_bait_names <- curated_PCHiC_vertex_with_alias %>%
    filter(type == "P") %>%
    select(gene_names) %>%
    mutate(all_bait_names = str_to_upper(gene_names)) %>%
    pull(all_bait_names)
  curated_bait_names <- sapply(all_bait_names, function(bait_name) {
    str_trim(str_remove_all(bait_name, "-\\d+\\b"))
  })

  curated_bait_names_unique <- sapply(curated_bait_names, function(curated_bait_name) {
    unique(str_split(curated_bait_name, fixed(" "))[[1]])
  })
  # curated_bait_names_unique <- sapply(curated_bait_names, function(curated_bait_name){ length(unique(str_split(curated_bait_name, fixed(" "))[[1]]))})
  names(curated_bait_names_unique) <- NULL
  # Original data with a list of promoters names
  # curated_PCHiC_vertex_with_alias$gene_names <- ifelse(curated_PCHiC_vertex_with_alias$type=="P",curated_bait_names_unique, curated_PCHiC_vertex_with_alias$gene_names)
  # Extract all promoters data to do unnest of the names and then the annotation with the symbols
  curated_bait_names_unique_df <- curated_PCHiC_vertex_with_alias %>% filter(type == "P")
  curated_bait_names_unique_df$all_bait_names <- curated_bait_names_unique
  curated_bait_names_unique_df_unnested <- NULL
  if (there_are_other_ends) {
    curated_bait_names_unique_df_unnested <- curated_bait_names_unique_df %>%
      dplyr::select(-c(mgi, ensembl, gene_type)) %>%
      unnest()
  } else {
    curated_bait_names_unique_df_unnested <- curated_bait_names_unique_df %>% tidyr::unnest()
  }
  alias_promoters <- alias %>%
    rename(`Gene name` = "all_bait_names") %>%
    mutate(all_bait_names = str_to_upper(all_bait_names)) %>%
    dplyr::select(-c(chr, start, end))
  promoters_merged_alias <- curated_bait_names_unique_df_unnested %>%
    mutate(gene_names = str_to_upper(all_bait_names)) %>%
    left_join(alias_promoters, by = "all_bait_names") %>%
    mutate(gene_names = str_to_sentence(all_bait_names))
  promoters_merged_alias_collapsed <- promoters_merged_alias %>%
    group_by(fragment) %>%
    summarise(
      gene_names = str_c(str_replace_na(gene_names), collapse = " "),
      ensembl = str_c(str_replace_na(`Ensembl gene ID`), collapse = " "),
      gene_type = str_c(str_replace_na(`Gene type`), collapse = " "),
      mgi = str_c(str_replace_na(`MGI ID`), collapse = " ")
    )

  promoters_merged_alias_collapsed$mgi <- str_remove_all(promoters_merged_alias_collapsed$mgi, fixed("MGI:"))

  # Order promoters according promoters_merged_alias_collapsed
  curated_PCHiC_vertex_with_alias <- curated_PCHiC_vertex_with_alias %>% arrange(fragment)

  curated_PCHiC_vertex_with_alias[curated_PCHiC_vertex_with_alias$type == "P", "mgi"] <- promoters_merged_alias_collapsed$mgi
  curated_PCHiC_vertex_with_alias[curated_PCHiC_vertex_with_alias$type == "P", "ensembl"] <- promoters_merged_alias_collapsed$ensembl
  curated_PCHiC_vertex_with_alias[curated_PCHiC_vertex_with_alias$type == "P", "gene_type"] <- promoters_merged_alias_collapsed$gene_type

  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$mgi) | curated_PCHiC_vertex_with_alias$mgi == "NA", "mgi"] <- ""

  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$ensembl) | curated_PCHiC_vertex_with_alias$ensembl == "NA", "ensembl"] <- ""
  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$gene_type) | curated_PCHiC_vertex_with_alias$gene_type == "NA", "gene_type"] <- ""
  curated_PCHiC_vertex_with_alias[is.na(curated_PCHiC_vertex_with_alias$gene_names), "gene_names"] <- ""

  curated_PCHiC_vertex_with_alias
}

#' @title generate_intronics_regions
#' @description Add intronic_region flag to curated vertex
#' @param curated_PCHiC_vertex curated vertex
#' @param intronic_regions_file intronic regions file
#' @return curated vertex with intronic regions flag
generate_intronics_regions <- function(curated_PCHiC_vertex, intronic_regions_file) {
  intronic_regions <- read_tsv(intronic_regions_file, col_types = cols(chr = col_character()))
  intronic_regions_grange <- makeGRangesFromDataFrame(intronic_regions)
  vertex_grange <- makeGRangesFromDataFrame(curated_PCHiC_vertex)
  overlaps <- findOverlaps(vertex_grange, intronic_regions_grange)
  query_hits <- unique(queryHits(overlaps))
  curated_PCHiC_vertex$intronic_regions <- F
  curated_PCHiC_vertex$intronic_regions[query_hits] <- ifelse(curated_PCHiC_vertex$type[query_hits] == "O", T, F)
  curated_PCHiC_vertex
}

#' @title generate_real_bait_names
#' @description Add real bait names to curated vertex
#' @param curated_PCHiC_vertex curated vertex
#' @param real_bait_names_file real bait names file
#' @return curated vertex with real bait names
generate_real_bait_names <- function(curated_PCHiC_vertex, real_bait_names_file) {
  real_bait_names <- read_tsv(real_bait_names_file, col_types = cols(Chr = col_character()))
  real_bait_names <- real_bait_names %>% mutate(gene_id = sapply(gene_id, function(g) {
    str_remove_all(g, "-\\d+\\b")[[1]][1]
  }))
  real_bait_names <- real_bait_names %>% mutate(gene_id = sapply(gene_id, function(g) {
    str_c(unique(str_split(g, fixed(","))[[1]]), collapse = " ")
  }))
  real_bait_names <- real_bait_names %>% mutate(fragment = str_c(Chr, Start, End, sep = "_"))

  merged_curated_PCHiC_vertex_with_baits <- left_join(curated_PCHiC_vertex, real_bait_names, by = "fragment")
  merged_curated_PCHiC_vertex_with_baits[merged_curated_PCHiC_vertex_with_baits$type == "P", "gene_names"] <- merged_curated_PCHiC_vertex_with_baits[merged_curated_PCHiC_vertex_with_baits$type == "P", "gene_id"]
  merged_curated_PCHiC_vertex_with_baits <- merged_curated_PCHiC_vertex_with_baits %>% select(-c(gene_id, ensembl_id, region, Chr, Start, End))

  merged_curated_PCHiC_vertex_with_baits
}
