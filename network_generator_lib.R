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
  return(parse_args(parser, args, convert_hyphens_to_underscores = T))
}
## ------------------------------------------------------------------------
search_vertex_by_name <-
  function(vertex, net, curated_chrs_vertex) {
    # Detect if we are searching by position (we are working with mouse chromosomes by now) or by name
    # Always return NULL if it doesn't exist the vertex in the graph
    if (str_detect(vertex, "^((1?[0-9])|([XY]))_\\d+$")) {
      if (!vertex %in% curated_chrs_vertex$fragment) {
        return(NULL)
      }
      required_vertex <- V(net)[vertex]
      required_subnet <-
        make_ego_graph(net, nodes = required_vertex)[[1]]
      return(required_subnet)
    } else {
      # Always search in lowercase
      vertex <- str_to_lower(vertex)
      if (!vertex %in% curated_chrs_vertex$curated_gene_name) {
        return(NULL)
      }
      searched_vertex_index <- curated_chrs_vertex$curated_gene_name == vertex
      required_vertex <- V(net)[searched_vertex_index]
      required_subnet <-
        make_ego_graph(net, nodes = required_vertex)[[1]]
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
    # make_ego_graph always returns a list
    required_subnet <-
      make_ego_graph(net, nodes = required_vertex)[[1]]
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
    }
  }
  return(required_subnet)
}
## Generate required subnetwork
search_subnetwork <- function(search, expand, nearest, net, curated_chrs_vertex) {
  if (!is.null(args$search)) {
    if (str_detect(args$search, "((1?[0-9])|([XY])):\\d+(-\\d+)?$")) {
      # We are working with a range
      required_subnet <-
        search_vertex_by_range(search, expand, nearest, net, curated_chrs_vertex)
    } else {
      required_subnet <-
        search_vertex_by_name(search, net, curated_chrs_vertex)
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
  vertices_df <-
    igraph::as_data_frame(required_subnet, what = "vertices")
  # Remove row names
  row.names(vertices_df) <- NULL
  # _ in column names is not valid in Cytoscape JSON
  colnames(vertices_df)[1] <- "id"
  colnames(vertices_df)[2] <- "names"
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
