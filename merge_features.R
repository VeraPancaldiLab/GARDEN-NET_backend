#!/usr/bin/env RScript
library(optparse)

parser <- OptionParser(description = "Merge features script")
parser <- add_option(parser, "--fifo_path", help = "Fifo shared with celery queue job system")
args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(parser, args, convert_hyphens_to_underscores = T)

print(args$fifo_path)
for (i in 1:10) {
  Sys.sleep(1)
  pipe(paste("echo 'TEST", i, "'>", args$fifo_path, sep = " "), "w")
}
pipe(paste("echo QUIT >", args$fifo_path, sep = " "), "w")
