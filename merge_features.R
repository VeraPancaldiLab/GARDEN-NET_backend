#!/usr/bin/env RScript
library(optparse)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(chaser))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rjson))
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

# Load Rdata from merge_features cache
load(file.path("/tmp/merge_features_cache", args$organism, args$cell_type, "merge_features_cache.Rdata"))
# Load Rdata from search_query cache
load(file.path("/tmp/merge_features_cache", args$organism, args$cell_type, "search_cache.Rdata"))

tmp_dir_path <- dirname(args$fifo_file)

# Calculate metadata
## Generate chas metadata
## Generate features metadata
# If there are  merged features
initial_features <- NULL
if (ncol(curated_PCHiC_vertex) > 6) {
  initial_features <- curated_PCHiC_vertex[, c(1, 7:length(curated_PCHiC_vertex))]
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
con <- pipe(paste("echo 'Generating features metadata for all network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
close(con)

net_features_metadata <- generate_features_metadata(PCHiC_ALL, curated_PCHiC_vertex, randomize = 10)
# PP network only
counter <- counter + 1
con <- pipe(paste("echo 'Generating features metadata for PP only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
close(con)
PCHiC_PP <- PCHiC_ALL[PCHiC_ALL$type == "P-P",]
  PCHiC_PP_fragment <- c(
    paste(PCHiC_PP$baitChr, PCHiC_PP$baitStart, sep = "_"),
    paste(PCHiC_PP$oeChr, PCHiC_PP$oeStart, sep = "_")
)
curated_PCHiC_vertex_PP <- curated_PCHiC_vertex[curated_PCHiC_vertex$fragment %in% PCHiC_PP_fragment,]
pp_net_features_metadata <- generate_features_metadata(PCHiC_PP, curated_PCHiC_vertex_PP)
# PO network only
counter <- counter + 1
con <- pipe(paste("echo 'Generating features metadata for PO only network:", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
close(con)
PCHiC_PO <- PCHiC_ALL[PCHiC_ALL$type == "P-O",]
  PCHiC_PO_fragment <- c(
    paste(PCHiC_PO$baitChr, PCHiC_PO$baitStart, sep = "_"),
    paste(PCHiC_PO$oeChr, PCHiC_PO$oeStart, sep = "_")
)
curated_PCHiC_vertex_PO <- curated_PCHiC_vertex[curated_PCHiC_vertex$fragment %in% PCHiC_PO_fragment,]
po_net_features_metadata <- generate_features_metadata(PCHiC_PO, curated_PCHiC_vertex_PO)

features_metadata <- list(net = net_features_metadata, pp = pp_net_features_metadata, po = po_net_features_metadata)
write(toJSON(features_metadata), file = file.path(tmp_dir_path, args$organism, args$cell_type, "features_metadata.json"))

counter <- counter + 1
# cl <- makeCluster(detectCores()-1)
# registerDoParallel(cl)
# foreach (chromosome=chromosomes, .packages=c("tibble", "igraph", "rjson","stringr")) %dopar% {
for (chromosome in chromosomes) {
  con <- pipe(paste("echo 'Processing chromosome ", chromosome, ":", counter / total * 100, "'>", args$fifo_file, sep = " "), "w")
  close(con)
  if (chromosome != "PP") {
    PCHiC_chr <- filter_by_chromosome(PCHiC_ALL, chromosome)
  } else {
    PCHiC_chr <- PCHiC_ALL[PCHiC_ALL$type == "P-P", ]
  }

  PCHiC_chr_fragment <- c(
    paste(PCHiC_chr$baitChr, PCHiC_chr$baitStart, sep = "_"),
    paste(PCHiC_chr$oeChr, PCHiC_chr$oeStart, sep = "_")
  )

  curated_PCHiC_vertex_chr <- curated_PCHiC_vertex[curated_PCHiC_vertex$fragment %in% PCHiC_chr_fragment,]
  
  json_curated_PCHiC_vertex <- generate_curated_PCHiC_vertex_json(curated_PCHiC_vertex_chr)
  json_curated_PCHiC_vertex <- str_replace_all(json_curated_PCHiC_vertex, fixed("\n"), "")
  json_curated_PCHiC_vertex <- str_replace_all(json_curated_PCHiC_vertex, fixed(", "), "")
  json_curated_PCHiC_vertex <- str_replace_all(json_curated_PCHiC_vertex, "\\{?\\s+\\}?", "")
  write(json_curated_PCHiC_vertex, file = file.path(tmp_dir_path, args$organism, args$cell_type, "chromosomes", paste0("chr", chromosome, ".json")))
  counter <- counter + 1
}

pipe(paste("echo QUIT >", args$fifo_file, sep = " "), "w")
