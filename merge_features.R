#!/usr/bin/env Rscript
library(optparse)
library(stringr)
library(tibble)
library(dplyr)
library(rjson)
suppressPackageStartupMessages(library(chaser))

source("./network_generator_lib.R")

parser <- OptionParser(description = "Merge features script")
parser <- add_option(parser, "--fifo_file", help = "Fifo shared with celery queue job system")
parser <- add_option(parser, "--organism", help = "Organism")
parser <- add_option(parser, "--cell_type", help = "Cell type")
parser <- add_option(parser, "--features_file", help = "Features file")
parser <- add_option(parser, "--features_file_type", help = "Features file type")
args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(parser, args, convert_hyphens_to_underscores = T)
# args <- parse_args(parser, args = c("--organism", "Mus_musculus", "--cell_type", "Embryonic_stem_cells"))
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Mon", "--features_file", "/tmp/ngoc/mon/S01RHSH1.ERX1305388.H3K27me3.bwa.GRCh38.broad.20160630.bed", "--features_file_type", "broad_peaks"))

# Load Rdata from merge_features cache
load(file.path("data", args$organism, args$cell_type, "merge_features_cache.Rdata"))
# Load Rdata from search_query cache
load(file.path("data", args$organism, args$cell_type, "search_cache.Rdata"))

if (!is.null(args$fifo_file)) {
  tmp_dir_path <- dirname(args$fifo_file)
}

counter <- 1

total <- 4

# All network
if (!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Adding features metadata for all network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}

# Remove old features first
chaser_net$features <- NULL
feature_name <- str_split(basename(args$features_file), fixed("."))[[1]][1]
chaser_net <- chaser::load_features(chaser_net, args$features_file, featname = feature_name, type = args$features_file_type, missingv = 0)

counter <- counter + 1

# All network
if (!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Generating features metadata for all network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}

net_features_metadata <- generate_features_metadata(chaser_net, randomize = 10)
# PP network only
counter <- counter + 1

if (!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Generating features metadata for PP only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}
pp_net_features_metadata <- generate_features_metadata(chaser::subset_chromnet(chaser_net, method = "bb"))

counter <- counter + 1

# PO network only
if (!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Generating features metadata for PO only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}
po_net_features_metadata <- generate_features_metadata(chaser::subset_chromnet(chaser_net, method = "bo"))

features_metadata <- list(net = net_features_metadata, pp = pp_net_features_metadata, po = po_net_features_metadata)

features <- as_tibble(chaser_net$features, rownames = "fragment")
colnames(features)[2] <- feature_name

features$fragment <- sapply(features$fragment, function(fragment) {
  str_remove(str_replace(str_split(fragment, fixed("-"))[[1]][1], fixed(":"), fixed("_")), fixed("chr"))
})

features <- features %>% select(fragment, everything())

features_for_json <- pull(features, feature_name)
names(features_for_json) <- features$fragment
json <- list(features_for_json)
names(json) <- feature_name

write(toJSON(json), file.path(tmp_dir_path, "features.json"))
write(toJSON(features_metadata), file.path(tmp_dir_path, "features_metadata.json"))

if (!is.null(args$fifo_file)) {
  pipe(paste("echo QUIT >", args$fifo_file, sep = " "), "w")
}
