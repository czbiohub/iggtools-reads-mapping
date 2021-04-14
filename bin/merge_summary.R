## Chunyu Zhao 2020-12-10
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(magrittr)


specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))


midas_dir <- args[1]
pstats_dir <- args[2]
o_sites_summary_fp <- args[3]
o_reads_summary_fp <- args[4]
species_under_investigation <- args[5]


bt2_levels <- data.frame(bt2_db = list.files(midas_dir, "^ani.*")) %>%
  separate(bt2_db, sep="\\.", into = c("_", "bt2_ani"), remove = F) %>% mutate(bt2_ani = as.numeric(bt2_ani)) %>%
  mutate(bt2_db = as.factor(bt2_db)) %>% mutate(bt2_db = fct_reorder(bt2_db, bt2_ani, .desc=TRUE)) %>% .$bt2_db %>% levels()


list_of_sites_summary <- list()
list_of_reads_summary <- list()

for (sim_cov in 1:50) {
  curr_reads_list <- list()
  curr_sites_list <- list()
  for (bt2_db in bt2_levels) {
    sample_name <- paste("art_", bt2_db, "_", sim_cov, "X", sep="")
    
    midas_summary_fp <- file.path(midas_dir, bt2_db, "out", sample_name, "snps", "snps_summary.tsv")
    summary <- read_delim(midas_summary_fp, delim = "\t", col_types= cols()) %>% mutate(sim_cov = sim_cov)
    curr_reads_list[[bt2_db]] <- summary
    
    pileup_stats_fp <- file.path(pstats_dir, paste(sample_name, ".tsv", sep=""))
    stats <- read_delim(pileup_stats_fp, delim = "\t", col_types= cols()) %>% mutate(sim_cov = sim_cov)
    curr_sites_list[[bt2_db]] <- stats
  }
  list_of_sites_summary[[as.character(sim_cov)]] <- bind_rows(curr_sites_list, .id = "bt2_db")
  list_of_reads_summary[[as.character(sim_cov)]] <- bind_rows(curr_reads_list, .id = "bt2_db")
}


## Add total simulated read counts
reads_summary <- bind_rows(list_of_reads_summary, .id = "sim_cov") %>% mutate(sim_cov = as.numeric(sim_cov)) %>%
  mutate(bt2_db = factor(bt2_db, levels = bt2_levels)) %>% mutate(species_id = as.character(species_id))
genome_length = unique(reads_summary %>% filter(species_id == species_under_investigation) %>% .$genome_length)
read_length = 125
readcounts_df <- tibble(sim_cov = 1:50, species_id = as.character(species_under_investigation)) %>% mutate(sim_readcounts = ceiling(genome_length * sim_cov / read_length))


reads_summary %<>%
  left_join(readcounts_df %>% select(sim_cov, sim_readcounts), by=c("sim_cov")) %>%
  mutate(percentage_aligned_reads = specify_decimal(aligned_reads / sim_readcounts, 3)) %>%
  mutate(percentage_piled_reads = specify_decimal(mapped_reads / sim_readcounts, 3)) %>%
  mutate(relative_vertical_coverage = specify_decimal(mean_coverage / sim_cov, 3)) %>%
  dplyr::rename(fraction_covered_sites  = fraction_covered)
reads_summary %>% write.table(o_reads_summary_fp, sep = "\t", quote=F, row.names = F)


bind_rows(list_of_sites_summary, .id = "sim_cov") %>%
  write.table(o_sites_summary_fp, sep="\t", quote = F, row.names = F)
