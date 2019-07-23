#!/usr/bin/env Rscript

# Load all libraries at the beginning
library(optparse)
library(stringr)
library(tibble)
library(dplyr)
library(rjson)
suppressPackageStartupMessages(library(chaser))

# Never use scientific format for numbers
options(scipen = 999)

source("./network_generator_lib.R")

# Define input command line parameters
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
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Neu", "--features_file", "/tmp/ngoc/mon/S01RHSH1.ERX1305388.H3K27me3.bwa.GRCh38.broad.20160630.bed", "--features_file_type", "macs2"))
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Neu", "--features_file", "/tmp/feat_ZZZ3.bed", "--features_file_type", "bed6"))
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Mac2", "--features_file", "/mnt/SERVER-CRCT-STORAGE/CRCT21/Private/Ngoc/Ngoc/Carrillo2017_chromatinstatesBP/SEGMENTATION/C002TWH2_11_segments.bed.gz", "--features_file_type", "chromhmm"))
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Mac2", "--features_file", "/mnt/SERVER-CRCT-STORAGE/CRCT21/Private/Common/forGARDEN-NET/BPfeatures/mono27acproc.bed.gz", "--features_file_type", "bed3"))
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Neu", "--features_file", "/tmp/mono-meth.bed6.gz", "--features_file_type", "bed6"))


# Load Rdata from merge_features cache
load(file.path("data", args$organism, args$cell_type, "merge_features_cache.Rdata"))

# Only use and write to the fifo_file when it exists
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

# If the feature format function is empty or is define from the frontend as None that means the default NULL value has to be used
if (args$feature_format_function == "" || args$feature_format_function == "None") {
  args$feature_format_function <- NULL
}

# Define good default return status at the beginning and define a bad status in the catch block
return_status <- 0
tryCatch({
  # Load user features to the chaser network with the parameters customized by the user
  chaser_net <- chaser::load_features(chaser_net, args$features_file, featname = feature_name, type = args$features_file_type, missingv = 0, auxfun = args$feature_format_function)

  counter <- counter + 1

  # All network
  if (!is.null(args$fifo_file)) {
    con <- pipe(paste("echo 'Generating features metadata for whole network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
    close(con)
  }

  # Define bait nodes for the random preservation and for the promoter-promoter only statistics
  baits <- chaser::export(chaser_net, "baits")
  net_features_metadata <- generate_features_metadata(chaser_net, randomize = 10, preserve.nodes = baits)
  counter <- counter + 1

  if (!is.null(args$fifo_file)) {
    con <- pipe(paste("echo 'Generating features metadata for PP only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
    close(con)
  }
  # PP network only
  chaser_net_bb <- chaser::subset_chromnet(chaser_net, method = "nodes", nodes1 = baits)
  pp_net_features_metadata <- generate_features_metadata(chaser_net_bb)

  counter <- counter + 1

  if (!is.null(args$fifo_file)) {
    con <- pipe(paste("echo 'Generating features metadata for PO only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
    close(con)
  }
  # PO network only
  all_oes <- chaser::export(chaser_net, "nodes")$name
  oes <- all_oes[!(all_oes %in% baits)]
  chaser_net_bo <- chaser::subset_chromnet(chaser_net, method = "nodes", nodes1 = baits, nodes2 = oes)
  po_net_features_metadata <- generate_features_metadata(chaser_net_bo)

  features_metadata <- list(net = net_features_metadata, pp = pp_net_features_metadata, po = po_net_features_metadata)

  features <- as_tibble(chaser_net$features, rownames = "fragment")
  # Use feature name when the chaser format provides it
  if (args$features_file_type == "macs2" || args$features_file_type == "bed3" || args$features_file_type == "bed6") {
    colnames(features)[2] <- feature_name
  }

  # Convert from genomic positions to node positions
  features$fragment <- sapply(features$fragment, function(fragment) {
    str_remove(str_replace_all(fragment, "[-:]", fixed("_")), fixed("chr"))
  })

  # Always the fragment column has to be the first
  features <- features %>% select(fragment, everything())

  # Convert features tibble to json
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
