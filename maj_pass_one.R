## Chunyu Zhao 2020-10-21
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



pileup_to_major <- function(pileup_fp, min_allele_depth=2) {
  #' Read in MIDAS's pileup result and generate the major alleles for each covered sites
  #' callable allele: allele being covered by at least min_allele_depth reads
  
  if (! file.exists(pileup_fp)) {
    system(paste("lz4 -d -c", paste(pileup_fp, ".lz4", sep=""), " > ", pileup_fp, ep="")) 
  }
  
  
  ## Ignore alleles with less than min_allele_depth (2 reads)
  pile.long <- read_delim(pileup_fp, delim = "\t", col_types = cols()) %>% 
    gather(allele, counts, count_a:count_t) %>% filter(counts >= min_allele_depth)
  
  ## Compute number_of_allels per genomic site after the previous filter
  sites <- pile.long %>% group_by(ref_id, ref_pos) %>% dplyr::count() %>% ungroup() %>% dplyr::rename(allele_counts = n)
  sites %>% dplyr::count(allele_counts)
  
  
  major <- pile.long %>% group_by(ref_id, ref_pos) %>% filter(counts == max(counts)) %>% filter(n() == 1) %>% ungroup() %>%
    mutate(allele = toupper(sub("count_", "", allele)))
  
  
  major %<>% left_join(sites, by=c("ref_id", "ref_pos"))
  stopifnot(nrow(filter(major, is.na(allele_counts))) == 0)
  
  return(major)    
}



specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))



gen_maj_per_simcov <- function(sim_cov, bt2_levels, cal_maj_fp, stats_fp, min_allele_depth=2) {
  
  site_stats <- data.frame(bt2_db = bt2_levels)
  
  ######################################## pileup => maj.long
  maj_list <- list()
  for (bt2_db in bt2_levels) {
    sample_name <- paste("art_", bt2_db, "_", sim_cov, "X", sep="")
    pileup_fp <- file.path(datadir, "4_midas", bt2_db, "out", sample_name, "snps", paste(species_under_investigation, ".snps.tsv", sep=""))
    maj_list[[as.character(bt2_db)]] <- pileup_to_major(pileup_fp, min_allele_depth)
  }
  maj_long <- bind_rows(maj_list, .id = "bt2_db") %>% select(bt2_db, everything())
  
  
  ## How many sites are genotyped?
  cs <- maj_long %>% group_by(bt2_db) %>% dplyr::count() %>% ungroup() %>% dplyr::rename(total_genotyped_sites = n)
  if (nrow(cs) > 0) {
    site_stats <- left_join(site_stats, cs, by=c("bt2_db"))
  }
  
  
  # How many genotyped sites are not covered by the true_sites? [could from non-uniquely aligned regions] (repeats)
  maj_long %<>% left_join(true_sites, by=c("ref_id" = "refid", "ref_pos" = "pos1", "ref_allele" = "refallele"))
  cs <- maj_long %>% filter(is.na(sitetype)) %>% group_by(bt2_db) %>% summarise(sites_outside = n()) 
  if (nrow(cs) > 0) {
    site_stats <- left_join(site_stats, cs, by=c("bt2_db"))
  }
  
  maj_long %<>% filter(!is.na(sitetype))
  
  
  ## Among all the genotyped sites overlapped with true_sites, how many are not accurate?
    #### For `ref` sites, major_allele == refallele == altallele; For `snp` sites, major_allele == altallele
    #### Therefore, we use major_allele != altallele 
    #### Two possibilities: (1) true_sites were wrong (2) pileup errors
  incorrect_sites <- maj_long %>% filter(allele != altallele)
  cs <- incorrect_sites %>% dplyr::count(bt2_db, sitetype) %>% mutate(sitetype = paste("incorrect_sites_", sitetype, sep="")) %>% spread(sitetype, n, fill = 0)
  if (nrow(cs) > 0) {
    site_stats <- left_join(site_stats, cs, by=c("bt2_db"))
  }
  
  incorrect_sites %<>% group_by(ref_id, ref_pos, qryid, pos2) %>% summarise(num_bt2dbs = n()) %>% ungroup()
  
  
  ## For some sites that was disagreed by one or more bt2 but not all, we discard that sites for downstream analysis
  maj_long <- left_join(maj_long, incorrect_sites, by=c("ref_id", "ref_pos", "qryid", "pos2"))
  discard_sites <- maj_long  %>% filter(!is.na(num_bt2dbs)) 
  cs <- discard_sites %>% dplyr::count(bt2_db, sitetype) %>% mutate(sitetype = paste("discard_sites_", sitetype, sep="")) %>% spread(sitetype, n, fill = 0)
  if (nrow(cs) > 0) {
    site_stats <- left_join(site_stats, cs, by=c("bt2_db"))
  }
  
  
  tmpdir <- file.path(dirname(cal_maj_fp))
  if (nrow(discard_sites)) {
    discard_sites %>% write.table(file.path(tmpdir, paste("discard_sites_", sim_cov, "X.tsv", sep="")), sep="\t", quote = F, row.names = F)
  }
  
  
  maj_long %<>% filter(is.na(num_bt2dbs)) %>% select(-num_bt2dbs)
  stopifnot(nrow(maj_long %>% filter(allele != altallele)) == 0)
  
  ## Did we genotype the correct allele?
  n1 <- maj_long %>% filter(ref_allele == allele & sitetype != "ref") %>% nrow()
  if (n1 > 0) {
    maj_long %>% filter(ref_allele == allele & sitetype != "ref") %>%
      write.table(file.path(dirname(cal_maj_fp), paste("incor_genotyped_ref_sites_", sim_cov, "X.tsv", sep="")), sep="\t", quote = F, row.names = F)
    
    maj_long %<>% filter(ref_allele == allele & sitetype == "ref")
  }
  n2 <- maj_long %>% filter(ref_allele != allele & altallele == allele & sitetype != "snp") %>% nrow()
  if (n2 > 0) {
    maj_long %>% filter(ref_allele != allele & altallele == allele & sitetype != "snp") %>% 
      write.table(file.path(dirname(cal_maj_fp), paste("incor_genotyped_snp_sites_", sim_cov, "X.tsv", sep="")), sep="\t", quote = F, row.names = F)
    
    maj_long %<>% filter(ref_allele != allele & altallele == allele & sitetype == "snp") 
  } 
  
  
  cs <- maj_long %>% group_by(bt2_db, sitetype) %>% dplyr::count() %>% spread(sitetype, n, fill=0)
  if (nrow(cs) > 0) {
    site_stats <- left_join(site_stats, cs, by=c("bt2_db"))
  }
  
  ######################################## maj.long => maj.wide: potentially this can be ported into MIDAS
  ## AMONG all the detected true sites by at least one bt2db, how many are accurate?
  maj_wide <- maj_long %>% select(bt2_db, ref_id, ref_pos, ref_allele, counts) %>% spread(bt2_db, counts, fill = 0)
  
  if (FALSE) {
    add_freq <- maj_long %>% left_join(maj_wide %>% select(ref_id, ref_pos) %>% mutate(pass = TRUE), by=c("ref_id", "ref_pos")) %>% filter(pass == TRUE) %>% 
      mutate(allele_frequency = specify_decimal(counts/depth,3)) %>%
      select(bt2_db, ref_id, ref_pos, ref_allele, allele_frequency) %>%
      mutate(bt2_db = paste("allele_freq_", bt2_db, sep="")) %>%
      spread(bt2_db, allele_frequency, fill = 0)
    stopifnot(nrow(maj_wide) == nrow(add_freq))
  }
  
  add_depth <- maj_long %>% left_join(maj_wide %>% select(ref_id, ref_pos) %>% mutate(pass = TRUE), by=c("ref_id", "ref_pos")) %>% filter(pass == TRUE) %>% 
    select(bt2_db, ref_id, ref_pos, ref_allele, depth) %>%
    mutate(bt2_db = paste("depth_", bt2_db, sep="")) %>%
    spread(bt2_db, depth, fill = 0)
  stopifnot(nrow(maj_wide) == nrow(add_depth))
  
  
  maj_wide %<>% left_join(add_depth, by=c("ref_id", "ref_pos", "ref_allele"))
  maj_wide %<>% left_join(true_sites, by=c("ref_id"="refid", "ref_pos" = "pos1", "ref_allele" = "refallele")) %>% filter(!is.na(sitetype))
  

  site_stats %>% write.table(stats_fp, sep="\t", row.names = F, quote = F)
  maj_wide %>% write.table(cal_maj_fp, sep="\t", row.names = F, quote = F)
}






####################### main functions #######################
## Input
if (length(args) > 0 ) {
  datadir <- args[1]
  species_under_investigation <- args[2]
  aligned_sites_fp <- args[3]
  sim_cov <- args[4]
  min_allele_depth = 2
} else {
  species_under_investigation <- 103703
  datadir <- file.path("/Users/chunyu.zhao/Documents/20200729_snps_reads_mapping/20200917_prep_s3", species_under_investigation)
  sim_cov <- 20
  aligned_sites_fp <- file.path(datadir, "5_nucmer", paste(species_under_investigation, "_aligned_sites.tsv", sep=""))
  min_allele_depth = 2
}


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

true_sites <- true_sites_from_aligned(aligned_sites_fp) %>% filter(sitetype!= "mism")

## Compute
gen_maj_per_simcov(sim_cov, bt2_levels, cal_maj_fp, stats_fp, min_allele_depth=2) 

