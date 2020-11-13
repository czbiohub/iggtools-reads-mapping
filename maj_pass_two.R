## Chunyu Zhao 2020-10-21
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(magrittr)



specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))



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


maj_pass_two <- function(r_outdir, species_under_investigation, bt2_levels, midas_dir, nucmer_dir) {
  #' MERGE maj_wide => calculate stats
  
  species_outdir <- file.path(r_outdir, species_under_investigation)
  
  o_det_power_fp <- file.path(r_outdir, paste(species_under_investigation, "_detection_power.tsv", sep=""))
  o_site_presence_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_presence.tsv", sep=""))
  o_site_abuncorr_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_abuncorr.tsv", sep=""))
  o_site_abun_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_abundance.tsv", sep=""))
  o_site_af_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_allelefreq.tsv", sep=""))
  o_reads_summary_fp <- file.path(r_outdir, paste(species_under_investigation, "_reads_summary.tsv", sep=""))
  o_true_sites_fp <- file.path(r_outdir, paste(species_under_investigation, "_true_sites.tsv", sep=""))
  
  
  
  if (! file.exists(o_true_sites_fp)) {
    aligned_sites_tp <- file.path(nucmer_dir, paste(species_under_investigation, "_aligned_sites.tsv", sep=""))
    true_sites <- true_sites_from_aligned(aligned_sites_tp) %>% filter(sitetype != "mism")
    true_sites %>% write.table(o_true_sites_fp, sep="\t", quote = F, row.names = F)
  }
  
  
  
  list_of_detection_power <- list()
  site_presence_list <- list()
  site_abuncorr_list <- list()
  site_abun_list <- list()
  site_af_list <- list()
  all_sum_list <- list() # Read in the summary file from iggtools midas_run_snps

  
  for (sim_cov in 1:50) {
    
    stats_fp <- file.path(species_outdir, paste("maj.stats_", sim_cov, "X.tsv", sep=""))
    site_presence_fp <- file.path(species_outdir, paste("sites_presence_", sim_cov, "X.tsv", sep=""))
    site_abuncorr_fp <- file.path(species_outdir, paste("sites_abuncorr_", sim_cov, "X.tsv", sep=""))
    site_abun_fp <- file.path(species_outdir, paste("sites_abundance_", sim_cov, "X.tsv", sep=""))
    site_af_fp <- file.path(species_outdir, paste("sites_allelefreq_", sim_cov, "X.tsv", sep=""))
    
    stopifnot(file.exists(stats_fp))
    stopifnot(file.exists(site_af_fp))
    
    
    
    #### Read in the summary file from iggtools midas_run_snps
    curr_sum_list <- list()
    for (bt2_db in bt2_levels) {
      sample_name <- paste("art_", bt2_db, "_", sim_cov, "X", sep="")
      summary_fp <- file.path(midas_dir, bt2_db, "out", sample_name, "snps", "snps_summary.tsv")
      summary <- read_delim(summary_fp, delim = "\t") %>% mutate(sim_cov = sim_cov)
      curr_sum_list[[bt2_db]] <- summary
    }
    all_sum_list[[as.character(sim_cov)]] <- bind_rows(curr_sum_list, .id = "bt2_db")
    
    
    
    curr_stats <- read_delim(stats_fp, delim = "\t", col_types = cols())
    list_of_detection_power[[as.character(sim_cov)]] <- curr_stats
    
    sites_abun <- read_delim(site_abun_fp, delim = "\t", col_types = cols())
    site_abun_list[[as.character(sim_cov)]] <- sites_abun
    
    
    sites_corr <- read_delim(site_abuncorr_fp, delim = "\t", col_types = cols())
    site_abuncorr_list[[as.character(sim_cov)]] <- sites_corr
    
    
    sites_presence <- read_delim(site_presence_fp, delim = "\t", col_types = cols())
    site_presence_list[[as.character(sim_cov)]] <- sites_presence
    
    
    sites_af <- read_delim(site_af_fp, delim = "\t", col_types = cols())
    site_af_list[[as.character(sim_cov)]] <- sites_af
  }
  
  
  reads_summary <- bind_rows(all_sum_list, .id = "sim_cov") %>%
    mutate(sim_cov = as.numeric(sim_cov)) %>% mutate(bt2_db = factor(bt2_db, levels = bt2_levels)) %>% mutate(species_id = as.character(species_id))
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
  
  
  bind_rows(list_of_detection_power, .id="sim_cov") %>% write.table(o_det_power_fp, sep="\t", quote = F, row.names = F)
  
  
  bind_rows(site_abuncorr_list, .id = "sim_cov") %>% write.table(o_site_abuncorr_fp, sep="\t", quote = F, row.names = F)
  
  
  bind_rows(site_abun_list, .id = "sim_cov") %>% write.table(o_site_abun_fp, sep="\t", quote = F, row.names = F)
  
  
  bind_rows(site_presence_list, .id = "sim_cov") %>% write.table(o_site_presence_fp, sep="\t", quote = F, row.names = F)
  
  
  bind_rows(site_af_list, .id = "sim_cov") %>% write.table(o_site_af_fp, sep="\t", quote = F, row.names = F)
  
}


####################### main functions #######################
if (length(args) > 0 ) {
  datadir <- args[1]
  species_under_investigation <- args[2]
  
  
  midas_dir <- file.path(datadir, "4_midas")
  nucmer_dir <- file.path(datadir, "5_nucmer")
  r_outdir <- file.path(datadir, "7_rdata")
  
  
  bt2_levels <- data.frame(bt2_db = list.files(file.path(datadir, "4_midas"))) %>% 
    separate(bt2_db, sep="\\.", into = c("_", "bt2_ani"), remove = F) %>% mutate(bt2_ani = as.numeric(bt2_ani)) %>% 
    mutate(bt2_db = as.factor(bt2_db)) %>% mutate(bt2_db = fct_reorder(bt2_db, bt2_ani, .desc=TRUE)) %>% .$bt2_db %>% levels()
  
  
  maj_pass_two(r_outdir, species_under_investigation, bt2_levels, midas_dir, nucmer_dir)
  
} else {
  datadir <- "~/Desktop/done/100107"
  species_under_investigation <- "100107"
}

