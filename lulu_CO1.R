#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

args[1]
args[2]
args[3]
args[4]

library(lulu)

otutable <- read.csv(args[1], sep = '\t', header = TRUE, as.is = TRUE, row.names = 1)

colnames(otutable) <- gsub(".CO1", "", gsub("X", "site", colnames(otutable)))

matchlist <- read.table(args[2], header = FALSE, as.is = TRUE, stringsAsFactors = FALSE)

lulu_result <- lulu(otutable, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 90, minimum_relative_cooccurence = 0.95)

write.csv(lulu_result$curated_otus, file = args[3])

write.csv(lulu_result$curated_table, file = args[4])
