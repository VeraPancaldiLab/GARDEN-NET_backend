#!/usr/bin/env Rscript
library(optparse)
suppressPackageStartupMessages(library(chaser))
# suppressPackageStartupMessages(library(doParallel))

source("./network_generator_lib.R")

parser <- OptionParser(description = "Merge features script")
parser <- add_option(parser, "--fifo_file", help = "Fifo shared with celery queue job system")
parser <- add_option(parser, "--organism", help = "Organism")
parser <- add_option(parser, "--cell_type", help = "Cell type")
parser <- add_option(parser, "--features_file", help = "Features file")
args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(parser, args, convert_hyphens_to_underscores = T)
# args <- parse_args(parser, args = c("--organism", "Mus_musculus", "--cell_type", "Embryonic_stem_cells"))
# args <- parse_args(parser, args = c("--organism", "Homo_sapiens", "--cell_type", "Mon", "--features_file", "/tmp/ngoc/mono/S01RHSH1.ERX1305388.H3K27me3.bwa.GRCh38.broad.20160630.bed"))

# Load Rdata from merge_features cache
load(file.path("/tmp/merge_features_cache", args$organism, args$cell_type, "merge_features_cache.Rdata"))
# Load Rdata from search_query cache
load(file.path("/tmp/merge_features_cache", args$organism, args$cell_type, "search_cache.Rdata"))

if(!is.null(args$fifo_file)) {
  tmp_dir_path <- dirname(args$fifo_file)
}

# Generate chromosomes
chromosomes <- c("X", "Y", "PP")

if (args$organism == "Mus_musculus") {
  chromosomes <- c(1:19, chromosomes)
} else if (args$organism == "Homo_sapiens") {
  chromosomes <- c(1:22, chromosomes)
}

counter <- 1

total <- length(chromosomes) + 3

# All network
if(!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Generating features metadata for all network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}

# TODO: select right features format file
# chaser_net <- chaser::load_features(chaser_net, args$features_file, featname="H3K27me3", type="broad_peaks", missingv=0)
# chaser_net <- chaser::load_features(chaser_net, args$features_file , type="data.frame", missingv=0)

net_features_metadata <- generate_features_metadata(chaser_net, randomize = 10)
# PP network only
counter <- counter + 1

if(!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Generating features metadata for PP only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}
pp_net_features_metadata <- generate_features_metadata(chaser::subset_chromnet(chaser_net, method="bb"))
# PO network only
counter <- counter + 1

if(!is.null(args$fifo_file)) {
  con <- pipe(paste("echo 'Generating features metadata for PO only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
}
po_net_features_metadata <- generate_features_metadata(chaser::subset_chromnet(chaser_net, method="bo"))

features_metadata <- list(net = net_features_metadata, pp = pp_net_features_metadata, po = po_net_features_metadata)
# write(toJSON(features_metadata), file = file.path(tmp_dir_path, args$organism, args$cell_type, "features_metadata.json"))

counter <- counter + 1
# cl <- makeCluster(detectCores()-1)
# registerDoParallel(cl)
# foreach (chromosome=chromosomes, .packages=c("tibble", "igraph", "rjson","stringr")) %dopar% {
for (chromosome in chromosomes) {
  if(!is.null(args$fifo_file)) {
    con <- pipe(paste("echo 'Processing chromosome ", chromosome, ":", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
    close(con)
  }
  if (chromosome != "PP") {
    chaser_net_chr <- chaser::subset_chromnet(chaser_net, chrom=paste0("chr", chromosome), interchrom=T)
  } else {
    chaser_net_chr <- chaser::subset_chromnet(chaser_net, method="bb")
  }

  # TODO: generate node + features file
  # print(chaser_net_chr$features)
  # chaser_net_json <- generate_cytoscape_json(chaser::export(chaser_net, type="igraph"))
  # chaser_net_json <- str_replace_all(chaser_net_json, fixed("\n"), "")
  # chaser_net_json <- str_replace_all(chaser_net_json, fixed(", "), "")
  # chaser_net_json <- str_replace_all(chaser_net_json, "\\{?\\s+\\}?", "")
  # write(chaser_net_json, file = file.path(tmp_dir_path, args$organism, args$cell_type, "chromosomes", paste0("chr", chromosome, ".json")))
  counter <- counter + 1
}

if(!is.null(args$fifo_file)) {
  pipe(paste("echo QUIT >", args$fifo_file, sep = " "), "w")
}
