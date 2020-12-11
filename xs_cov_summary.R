## Chunyu Zhao 2020-12-10
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(magrittr)


merge_summary_xs_coverage <- function(r_outdir, species_under_investigation, bt2_levels, midas_dir) {
  
  species_outdir <- file.path(r_outdir, species_under_investigation)
  
  o_sites_summary_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_summary.tsv", sep=""))
  o_reads_summary_fp <- file.path(r_outdir, paste(species_under_investigation, "_reads_summary.tsv", sep=""))
  
  list_of_sites_summary <- list()
  list_of_reads_summary <- list()
  
  for (sim_cov in 1:50) {
    
    site_sum_fp <- file.path(species_outdir, paste("sites_summary_", sim_cov, "X.tsv", sep=""))
    list_of_sites_summary[[as.character(sim_cov)]] <- read_delim(site_sum_fp, delim = "\t", col_types = cols())
  
    #### Read in the summary file from iggtools midas_run_snps
    curr_sum_list <- list()
    for (bt2_db in bt2_levels) {
      sample_name <- paste("art_", bt2_db, "_", sim_cov, "X", sep="")
      summary_fp <- file.path(midas_dir, bt2_db, "out", sample_name, "snps", "snps_summary.tsv")
      summary <- read_delim(summary_fp, delim = "\t") %>% mutate(sim_cov = sim_cov)
      curr_sum_list[[bt2_db]] <- summary
    }
    list_of_reads_summary[[as.character(sim_cov)]] <- bind_rows(curr_sum_list, .id = "bt2_db")
  }
  
  
  reads_summary <- bind_rows(list_of_reads_summary, .id = "sim_cov") %>% mutate(sim_cov = as.numeric(sim_cov)) %>% 
    mutate(bt2_db = factor(bt2_db, levels = bt2_levels)) %>% mutate(species_id = as.character(species_id))
  ## Add total simulated read counts
  genome_length = unique(reads_summary %>% filter(species_id == species_under_investigation) %>% .$genome_length)
  read_length = 125
  readcounts_df <- tibble(sim_cov = 1:50, species_id = as.character(species_under_investigation)) %>% mutate(sim_readcounts = ceiling(genome_length * sim_cov / read_length))
  
  reads_summary %<>% 
    left_join(readcounts_df %>% select(sim_cov, sim_readcounts), by=c("sim_cov")) %>% 
    mutate(percentage_aligned_reads = aligned_reads / sim_readcounts) %>% 
    mutate(percentage_piled_reads = mapped_reads / sim_readcounts) %>% 
    mutate(relative_vertical_coverage = mean_coverage / sim_cov) %>% 
    dplyr::rename(fraction_covered_sites  = fraction_covered)
  
  reads_summary %>% write.table(o_reads_summary_fp, sep = "\t", quote=F, row.names = F)
  
  bind_rows(list_of_sites_summary, .id="sim_cov") %>% write.table(o_sites_summary_fp, sep="\t", quote = F, row.names = F)
  
}


####################### main functions #######################
if (length(args) > 0 ) {
  datadir <- args[1]
  species_under_investigation <- args[2]
  
  midas_dir <- file.path(datadir, "4_midas")
  r_outdir <- file.path(datadir, "6_rdata")
  
  
  bt2_levels <- data.frame(bt2_db = list.files(file.path(datadir, "4_midas"))) %>% 
    separate(bt2_db, sep="\\.", into = c("_", "bt2_ani"), remove = F) %>% mutate(bt2_ani = as.numeric(bt2_ani)) %>% 
    mutate(bt2_db = as.factor(bt2_db)) %>% mutate(bt2_db = fct_reorder(bt2_db, bt2_ani, .desc=TRUE)) %>% .$bt2_db %>% levels()
  
  merge_summary_xs_coverage(r_outdir, species_under_investigation, bt2_levels, midas_dir)
} else {
  datadir <- "~/Desktop/103703"
  species_under_investigation <- "103703"
}

