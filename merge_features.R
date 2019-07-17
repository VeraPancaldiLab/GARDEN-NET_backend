#!/usr/bin/env Rscript
library(optparse)
library(stringr)
library(tibble)
library(dplyr)
library(rjson)
suppressPackageStartupMessages(library(chaser))

# Never use scientific format for numbers
options(scipen = 999)

source("./network_generator_lib.R")

parser <- OptionParser(description = "Merge features script")
parser <- add_option(parser, "--fifo_file", help = "Fifo shared with celery queue job system")
parser <- add_option(parser, "--organism", help = "Organism")
parser <- add_option(parser, "--cell_type", help = "Cell type")
parser <- add_option(parser, "--features_file", help = "Features file")
parser <- add_option(parser, "--features_file_type", help = "Features file type")
parser <- add_option(parser, "--feature_format_function", default = "", help = "Features format function")
args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(parser, args, convert_hyphens_to_underscores = T)
# args <- parse_args(parser, args = c("--organism", "Mus_musculus", "--cell_type", "Embryonic_stem_cells", "--features_file", "/mnt/SERVER-CRCT-STORAGE/CRCT21/Private/Common/forGARDEN-NET/mean_efficiency_wt_mm9_good_1.bed", "--features_file_type", "bed6"))
# args <- parse_args(parser, args = c("--organism", "Mus_musculus", "--cell_type", "Embryonic_stem_cells", "--features_file", "/tmp/ngoc/mon/S01RHSH1.ERX1305388.H3K27me3.bwa.GRCh38.broad.20160630.bed", "--features_file_type", "macs2"))
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Mon", "--features_file", "/tmp/ngoc/mon/S01RHSH1.ERX1305388.H3K27me3.bwa.GRCh38.broad.20160630.bed", "--features_file_type", "macs2"))

# Load Rdata from merge_features cache
load(file.path("data", args$organism, args$cell_type, "merge_features_cache.Rdata"))

if (!is.null(args$fifo_file)) {
  tmp_dir_path <- dirname(args$fifo_file)
}

counter <- 1

total <- 5

# All network
if (!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Adding features metadata for all network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}

# Remove old features first
chaser_net$features <- NULL
# Never use dash - in feature names (forbidden in cytoscape attributes)
feature_name <- str_replace_all(str_split(basename(args$features_file), fixed("."))[[1]][1], fixed("-"), fixed("_"))

if (args$feature_format_function == "" || args$feature_format_function == "None") {
  args$feature_format_function <- NULL
}
return_status <- 0
tryCatch({
  chaser_net <- chaser::load_features(chaser_net, args$features_file, featname = feature_name, type = args$features_file_type, missingv = 0, auxfun = args$feature_format_function)

  counter <- counter + 1

  # All network
  if (!is.null(args$fifo_file)) {
    con <- pipe(paste("echo 'Generating features metadata for whole network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
    close(con)
  }

  baits <- unique(chaser::export(chaser_net, "originaldf")$node_from)
  net_features_metadata <- generate_features_metadata(chaser_net, randomize = 10, preserve.nodes = baits)
  # PP network only
  counter <- counter + 1

  if (!is.null(args$fifo_file)) {
    con <- pipe(paste("echo 'Generating features metadata for PP only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
    close(con)
  }
  chaser_net_bb <- chaser::subset_chromnet(chaser_net, method = "nodes", nodes1 = baits)
  pp_net_features_metadata <- generate_features_metadata(chaser_net_bb)

  counter <- counter + 1

  # PO network only
  if (!is.null(args$fifo_file)) {
    con <- pipe(paste("echo 'Generating features metadata for PO only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
    close(con)
  }
  all_oes <- unique(chaser::export(chaser_net, "originaldf")$node_to)
  oes <- all_oes[!(all_oes %in% baits)]
  chaser_net_bo <- chaser::subset_chromnet(chaser_net, method = "nodes", nodes1 = baits, nodes2 = oes)
  po_net_features_metadata <- generate_features_metadata(chaser_net_bo)

  features_metadata <- list(net = net_features_metadata, pp = pp_net_features_metadata, po = po_net_features_metadata)

  features <- as_tibble(chaser_net$features, rownames = "fragment")
  if (args$features_file_type == "macs2" || args$features_file_type == "bed6") {
    colnames(features)[2] <- feature_name
  }

  features$fragment <- sapply(features$fragment, function(fragment) {
    str_remove(str_replace_all(fragment, "[-:]", fixed("_")), fixed("chr"))
  })

  features <- features %>% select(fragment, everything())

  to_json <- list()
  for (index in 2:length(colnames(features))) {
    feature <- colnames(features)[index]
    features_for_json <- pull(features, feature)
    names(features_for_json) <- features$fragment
    to_json[[feature]] <- features_for_json
  }


  write(toJSON(to_json), file.path(tmp_dir_path, "features.json"))
  write(toJSON(features_metadata), file.path(tmp_dir_path, "features_metadata.json"))
},
error = function(error) {
  print(error)
  return_status <- 1
},
finally = {
  if (!is.null(args$fifo_file)) {
    pipe(paste("echo QUIT >", args$fifo_file, sep = " "), "w")
  }
  quit(status = return_status)
}
)
