## Chunyu Zhao 2020-12-09
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)


library(tidyverse)
library(magrittr)


true_sites_from_aligned <- function(aligned_sites_fp) {
  ## Compute true_sites from aligned_sites

  aligned_sites <- read_delim(aligned_sites_fp, delim = "\t", col_types = cols())

  ## ref sites that aligned to more than one qry sites. haven't looked into the reason for this.
  dupsites <- aligned_sites %>% group_by(refid, pos1) %>% filter(n() > 1) %>% ungroup() %>% select(refid, pos1) %>% unique()

  ## generate the true positive set
  true_sites <- aligned_sites %>% left_join(dupsites %>% mutate(isdup = T), by=c("refid", "pos1")) %>% filter(is.na(isdup)) %>% select(-isdup)

  ## ignore the dupsites for now...
  dupsites <- left_join(dupsites, aligned_sites, by=c("refid", "pos1")) %>% arrange(refid, pos1)

  return(true_sites)
}


aligned_sites_fp <- args[1]
true_sites_fp <- args[2]


true_sites <- true_sites_from_aligned(aligned_sites_fp) %>% filter(sitetype != "mism")
true_sites %>% write.table(true_sites_fp, sep="\t", quote = F, row.names = F)
