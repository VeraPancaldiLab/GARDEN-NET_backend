## ------------------------------------------------------------------------
search_vertex_by_name <-
  function(vertex, net, curated_chrs_vertex) {
    # Detect if we are searching by position (we are working with mouse chromosomes by now) or by name
    # Always return NULL if it doesn't exist the vertex in the graph
    if (str_detect(vertex, "^((1?[0-9])|([XY]))_\\d+$")) {
      if (!any(V(net)$name == vertex)) {
        return(NULL)
      }
      required_vertex <- V(net)[vertex]
      required_subnet <-
        make_ego_graph(net, nodes = required_vertex)[[1]]
      return(required_subnet)
    } else {
      searched_vertex_index <-
        curated_chrs_vertex$curated_gene_name == vertex
      if (!any(searched_vertex_index)) {
        return(NULL)
      }
      required_vertex <- V(net)[searched_vertex_index]
      required_subnet <-
        make_ego_graph(net, nodes = required_vertex)[[1]]
      return(required_subnet)
    }
  }

search_vertex_by_range <- function(search, expand, nearest, net, curated_chrs_vertex) {
  curated_chrs_vertex_ranges <-
    makeGRangesFromDataFrame(curated_chrs_vertex,
                             keep.extra.columns = T,
                             ignore.strand = FALSE)
  required_range <- GRanges(search)
  # Expand the selected range if it is required
  if (expand != 0) {
    start(required_range) <- start(required_range) - expand
    end(required_range) <- end(required_range) + expand
  }
  # Work with the nearest if it is required
  if (nearest) {
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
    required_subnet <-
      induced_subgraph(net, vids = required_vertex_with_neighbours)
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
  return(toJSON(JSON_df, indent = 2))
}
