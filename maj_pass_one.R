## Chunyu Zhao 2020-10-21
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)


library(tidyverse)
library(docstring)
library(readr)
library(stringr)
library(reshape2)
library(readxl)
library(magrittr)
library(forcats)
library(data.table)

library(ape)
library(vegan)
library(seqinr)
library("Biostrings")
library(DescTools)
library(data.table)
library(assertthat)
library(ggsci)


library(foreach)
library(doParallel)



pileup_to_major <- function(pileup_fp, min_allele_depth=2) {
  #' Read in MIDAS's pileup result and generate the major alleles for each covered sites
  #' min_allele_depth: define callable alleles as allele being covered by at least 2 reads
  
  if (! file.exists(pileup_fp)) {
    system(paste("lz4 -d -c", paste(pileup_fp, ".lz4", sep=""), " > ", pileup_fp, ep="")) 
  }
  
  pile.long <- read_delim(pileup_fp, delim = "\t", col_types = cols()) %>% gather(allele, counts, count_a:count_t) %>% filter(counts >= min_allele_depth)
  
  sites <- pile.long %>% group_by(ref_id, ref_pos) %>% dplyr::count() %>% ungroup() %>% dplyr::rename(allele_counts = n)
  sites %>% dplyr::count(allele_counts)
  
  
  major <- pile.long %>% group_by(ref_id, ref_pos) %>% filter(counts == max(counts))
  
  tags <- major %>% dplyr::count() %>% ungroup() %>% dplyr::rename(max_counts = n) %>% 
    mutate(site_type = ifelse(max_counts == 2, "undetm", "detm")) # <-- site_type: tied max allele counts
  tags %>% dplyr::count(max_counts, site_type)
  
  
  
  ## randomly assign allele for undetermined sites
  major %<>% filter(row_number() == 1) %>% ungroup() %>% arrange(ref_id, ref_pos) %>% mutate(allele = toupper(sub("count_", "", allele)))
  
  
  major %<>% left_join(sites, by=c("ref_id", "ref_pos"))
  stopifnot(nrow(filter(major, is.na(allele_counts))) == 0)
  
  
  major %<>% left_join(tags %>% select(-max_counts), by=c("ref_id", "ref_pos"))
  stopifnot(nrow(filter(major, is.na(site_type))) == 0)
  
  return(major)    
}



true_sites_from_aligned <- function(aligned_sites_tp) {
  ## Compute true_sites from aligned_sites
  
  aligned_sites <- read_delim(aligned_sites_tp, delim = "\t", col_types = cols())
  
  ## ref sites that aligned to more than one qry sites. haven't looked into the reason for this.
  dupsites <- aligned_sites %>% group_by(refid, pos1) %>% filter(n() > 1) %>% ungroup() %>% select(refid, pos1) %>% unique()
  
  ## generate the true positive set
  true_sites <- aligned_sites %>% left_join(dupsites %>% mutate(isdup = T), by=c("refid", "pos1")) %>% filter(is.na(isdup)) %>% select(-isdup)
  
  ## ignore the dupsites for now...
  dupsites <- left_join(dupsites, aligned_sites, by=c("refid", "pos1")) %>% arrange(refid, pos1)
  
  return(true_sites)
}



gen_maj_per_simcov <- function(sim_cov, bt2_levels, cal_maj_fp, stats_fp, min_allele_depth=2) {
  
  ######################################## pileup => maj.long
  maj_list <- list()
  for (bt2_db in bt2_levels) {
    sample_name <- paste("art_", bt2_db, "_", sim_cov, "X", sep="")
    pileup_fp <- file.path(datadir, "4_midas", bt2_db, "out", sample_name, "snps", paste(species_under_investigation, ".snps.tsv", sep=""))
    maj_list[[as.character(bt2_db)]] <- pileup_to_major(pileup_fp, min_allele_depth)
  }
  
  maj_long <- bind_rows(maj_list, .id = "bt2_db") %>%
    select(bt2_db, everything()) %>% filter(site_type == "detm") # <-- no longer to compute the site_type inside pileup_to_major() <= TODO
  
  
  site_stats <- maj_long %>% dplyr::count(bt2_db, site_type) %>% spread(bt2_db, n, fill = 0)
  
  
  #### WHEN major_allele != altallele, could be two possibilities: (1) true_sites were wrong (2) pileup errors
  dis_sites <- left_join(maj_long, true_sites, by=c("ref_id" = "refid", "ref_pos" = "pos1", "ref_allele" = "refallele")) %>% filter(allele != altallele)
  dis_sites_types <- dis_sites %>% dplyr::count(sitetype) %>% mutate(sitetype = paste("disagree_", sitetype, sep="")) %>% spread(sitetype, n)
  if (nrow(dis_sites_types) == 0 ) {
    dis_sites_types <- data.frame(disagree_ref = 0, disagree_snp = 0)
  }
  
  ######################################## maj.long => maj.wide: potentially this can be ported into MIDAS
  maj_wide <- maj_long %>% mutate(pid = paste(ref_id, ref_pos, ref_allele, sep="|")) %>% select(bt2_db, pid, ref_id, ref_pos, ref_allele, counts) %>% spread(bt2_db, counts, fill = 0)
  
  ## IGNORE: reported sites that are not **in** the true sites set, could from non-uniquely aligned regions (repeats).
  sites_ignored <- nrow(maj_wide %>% left_join(true_sites, by=c("ref_id"="refid", "ref_pos" = "pos1", "ref_allele" = "refallele")) %>% filter(is.na(sitetype)))
  
  ## AMONG all the detected true sites by at least one bt2db
  maj_wide <- left_join(true_sites, maj_wide, by=c("refid"="ref_id", "pos1" = "ref_pos", "refallele" = "ref_allele"))
  
  ## Detection Power <- snp detection power to start with is 50% ... 
  site_stats <- cbind(site_stats, maj_wide %>% filter(!is.na(pid)) %>% dplyr::count(sitetype) %>% spread(sitetype, n) %>% mutate(sites_outside = as.numeric(sites_ignored))) %>% 
    mutate(sites_disagree = as.numeric(dis_sites %>% select(ref_id, ref_pos) %>% unique() %>% nrow()))
  site_stats <- cbind(site_stats, dis_sites_types)
  
  maj_wide %<>% filter(!is.na(pid)) %>% select(-pid)
  
  #### WHAT TO DO with those disagreed_sites?
  discard_sites <- left_join(maj_wide, dis_sites %>% group_by(ref_id, ref_pos) %>% summarise(num_bt2dbs = n()) %>% ungroup(), by=c("refid" = "ref_id", "pos1" = "ref_pos")) %>% filter(! is.na(num_bt2dbs))
  
  site_stats$sites_discard <- discard_sites %>% select(refid, pos1) %>% unique() %>% nrow()
  
  maj_wide <- left_join(maj_wide, dis_sites %>% group_by(ref_id, ref_pos) %>% summarise(num_bt2dbs = n()) %>% ungroup(), by=c("refid" = "ref_id", "pos1" = "ref_pos")) %>% filter( is.na(num_bt2dbs)) %>% select(-num_bt2dbs)
  
  
  site_stats %>% write.table(stats_fp, sep="\t", row.names = F, quote = F)
  maj_wide %>% write.table(cal_maj_fp, sep="\t", row.names = F, quote = F)
}






####################### main functions #######################
## Input
datadir <- args[1]
species_under_investigation <- args[2]
aligned_sites_fp <- args[3]
sim_cov <- args[4]
min_allele_depth = 2


## Output
r_tempdir <- file.path(datadir, "6_rtemp")
r_outdir <- file.path(datadir, "7_rdata")
species_outdir <- file.path(r_outdir,species_under_investigation)

cal_maj_fp <- file.path(r_tempdir, paste("maj.wide_", sim_cov, "X.tsv", sep=""))
stats_fp <- file.path(species_outdir, paste("maj.stats_", sim_cov, "X.tsv", sep=""))


## Prepare
bt2_levels <- data.frame(bt2_db = list.files(file.path(datadir, "4_midas"))) %>% 
  separate(bt2_db, sep="\\.", into = c("_", "bt2_ani"), remove = F) %>% mutate(bt2_ani = as.numeric(bt2_ani)) %>% 
  mutate(bt2_db = as.factor(bt2_db)) %>% mutate(bt2_db = fct_reorder(bt2_db, bt2_ani, .desc=TRUE)) %>% .$bt2_db %>% levels()

true_sites <- true_sites_from_aligned(aligned_sites_fp)

## Compute
gen_maj_per_simcov(sim_cov, bt2_levels, cal_maj_fp, stats_fp, min_allele_depth=2) 



######################## PROTOTYPE WITH NOTES #############################

maj_pass_one <- function(r_tempdir, species_outdir, species_under_investigation, bt2_levels, aligned_sites_fp, numCores=4, min_allele_depth=2) {
  #' Accumulate useful information
  #' Analysis flow: pileup => major_long => major_wide
  
  true_sites <- true_sites_from_aligned(aligned_sites_fp)
  
  for (sim_cov in 50:1) {
    
    print(paste(species_under_investigation, ":", sim_cov))
    
    cal_maj_fp <- file.path(r_tempdir, paste("maj.wide_", sim_cov, "X.tsv", sep=""))
    stats_fp <- file.path(species_outdir, paste("maj.stats_", sim_cov, "X.tsv", sep=""))
    
    
    maj_list <- list()
    for (bt2_db in bt2_levels) {
      sample_name <- paste("art_", bt2_db, "_", sim_cov, "X", sep="")
      pileup_fp <- file.path(datadir, "4_midas", bt2_db, "out", sample_name, "snps", paste(species_under_investigation, ".snps.tsv", sep=""))
      maj_list[[as.character(bt2_db)]] <- pileup_to_major(pileup_fp, min_allele_depth)
    }
    
    
    
    maj_long <- bind_rows(maj_list, .id = "bt2_db") %>%
      select(bt2_db, everything()) %>% filter(site_type == "detm") # <-- no longer to compute the site_type inside pileup_to_major() <= TODO
    
    
    site_stats <- maj_long %>% dplyr::count(bt2_db, site_type) %>% spread(bt2_db, n, fill = 0)
    
    
    #### WHEN major_allele != altallele, could be two possibilities: (1) true_sites were wrong (2) pileup errors
    dis_sites <- left_join(maj_long, true_sites, by=c("ref_id" = "refid", "ref_pos" = "pos1", "ref_allele" = "refallele")) %>% filter(allele != altallele)
    dis_sites_types <- dis_sites %>% dplyr::count(sitetype) %>% mutate(sitetype = paste("disagree_", sitetype, sep="")) %>% spread(sitetype, n)
    
    
    ######################################## maj.long => maj.wide: potentially this can be ported into MIDAS
    maj_wide <- maj_long %>% mutate(pid = paste(ref_id, ref_pos, ref_allele, sep="|")) %>% select(bt2_db, pid, ref_id, ref_pos, ref_allele, counts) %>% spread(bt2_db, counts, fill = 0)
    
    ## IGNORE: reported sites that are not **in** the true sites set, could from non-uniquely aligned regions (repeats).
    sites_ignored <- nrow(maj_wide %>% left_join(true_sites, by=c("ref_id"="refid", "ref_pos" = "pos1", "ref_allele" = "refallele")) %>% filter(is.na(sitetype)))
    
    ## AMONG all the detected true sites by at least one bt2db
    maj_wide <- left_join(true_sites, maj_wide, by=c("refid"="ref_id", "pos1" = "ref_pos", "refallele" = "ref_allele"))
    
    ## Detection Power <- snp detection power to start with is 50% ... 
    site_stats <- cbind(site_stats, maj_wide %>% filter(!is.na(pid)) %>% dplyr::count(sitetype) %>% spread(sitetype, n) %>% mutate(sites_outside = as.numeric(sites_ignored))) %>% 
      mutate(sites_disagree = as.numeric(dis_sites %>% select(ref_id, ref_pos) %>% unique() %>% nrow()))
    site_stats <- cbind(site_stats, dis_sites_types)
    
    maj_wide %<>% filter(!is.na(pid)) %>% select(-pid)
    
    #### WHAT TO DO with those disagreed_sites?
    discard_sites <- left_join(maj_wide, dis_sites %>% group_by(ref_id, ref_pos) %>% summarise(num_bt2dbs = n()) %>% ungroup(), by=c("refid" = "ref_id", "pos1" = "ref_pos")) %>% filter(! is.na(num_bt2dbs))
    
    site_stats$sites_discard <- discard_sites %>% select(refid, pos1) %>% unique() %>% nrow()
    
    maj_wide <- left_join(maj_wide, dis_sites %>% group_by(ref_id, ref_pos) %>% summarise(num_bt2dbs = n()) %>% ungroup(), by=c("refid" = "ref_id", "pos1" = "ref_pos")) %>% filter( is.na(num_bt2dbs)) %>% select(-num_bt2dbs)
    
    
    site_stats %>% write.table(stats_fp, sep="\t", row.names = F, quote = F)
    maj_wide %>% write.table(cal_maj_fp, sep="\t", row.names = F, quote = F)
  }
}



maj_pass_one_parallel <- function(r_tempdir, species_outdir, species_under_investigation, bt2_levels, aligned_sites_fp, numCores=16, min_allele_depth=2) {
  # 20201022: this doesn't work on EC2: ssh 
  #' Accumulate useful information
  #' Analysis flow: pileup => major_long => major_wide
  
  ##https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
  ##https://privefl.github.io/blog/a-guide-to-parallelism-in-r/
  
  cl <- parallel::makeCluster(numCores)
  registerDoParallel(numCores)
  print(paste(species_under_investigation, "start"))
  
  foreach(sim_cov=50:1) %dopar% {
    true_sites <- true_sites_from_aligned(aligned_sites_fp)
    print(paste(species_under_investigation, "--", sim_cov))
    cal_maj_fp <- file.path(r_tempdir, paste("maj.wide_", sim_cov, "X.tsv", sep=""))
    stats_fp <- file.path(species_outdir, paste("maj.stats_", sim_cov, "X.tsv", sep=""))
    gen_maj_per_simcov(sim_cov, bt2_levels, cal_maj_fp, stats_fp) 
  }
  
  parallel::stopCluster(cl)
}


