#!/usr/bin/env RScript
library(optparse)

parser <- OptionParser(description = "Merge features script")
parser <- add_option(parser, "--temp_file", help = "Temporal file shared with celery queue job system")
args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(parser, args, convert_hyphens_to_underscores = T)


print(args$temp_file)
for (i in 1:10) {
 Sys.sleep(1)
 write(paste("TEST", i, "...", sep =" "), file = args$temp_file)
}
