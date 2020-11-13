## Chunyu Zhao 2020-10-21
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(magrittr)


specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))


gen_maj_analysis <- function(r_tempdir, r_outdir, species_under_investigation, bt2_levels, sim_cov) {
  ############## Compuate stats per maj per sim_cov 
  
  species_outdir <- file.path(r_outdir, species_under_investigation)
  
  
  site_presence_fp <- file.path(species_outdir, paste("sites_presence_", sim_cov, "X.tsv", sep=""))
  site_abuncorr_fp <- file.path(species_outdir, paste("sites_abuncorr_", sim_cov, "X.tsv", sep=""))
  site_abun_fp <- file.path(species_outdir, paste("sites_abundance_", sim_cov, "X.tsv", sep=""))
  site_af_fp <- file.path(species_outdir, paste("sites_allelefreq_", sim_cov, "X.tsv", sep=""))
  
  
  ######### Read in maj_wide from maj_pass_one
  cal_maj_fp <- file.path(r_tempdir, paste("maj.wide_", sim_cov, "X.tsv", sep=""))
  stopifnot(file.exists(cal_maj_fp))
  
  maj_wide <- read_delim(cal_maj_fp, delim = "\t", col_types = cols())
  
  
  ############### What to do with sites not in ani.100 but in other bt2_anis? What's the reason for this observation?
  #### does ani.100 miss out sites (FN) or the other bt2_indexes reports FP sites
  extra_sites <- maj_wide %>% filter(ani.100 == 0)
  extra_sites %<>% gather(ani, readcounts, bt2_levels) %>% filter(readcounts > 0) %>% group_by(ref_id, qryid, ref_pos, pos2, sitetype) %>% summarise(nz = n(), totalRC = sum(readcounts)) %>% ungroup() 
  extra_sites %>% group_by(sitetype) %>% summarise(counts = n(), mean_nz = mean(nz), med_nz = median(nz), med_rc = median(totalRC), mean_rc = mean(totalRC))
  extra_sites %>% dplyr::count(sitetype)
  ## WHAT TO DO with the extra_sites?
  
  
  ################ DOUBE CHECK ME
  maj_wide %<>% filter(ani.100 > 0) #<- need to ask Katie about this one
  
  
  ############### Part One: presence-absence: among the detected-true-sites by at least one bt2_index, what is the number of sites that have been missed/called across bt2_anis.
  list_of_presence <- lapply(bt2_levels, function(x) {
    maj_wide %>% select(sitetype, all_of(x)) %>% dplyr::rename(ani := !!x) %>% mutate(presence = ifelse(ani > 0, "called", "missed")) %>% dplyr::count(sitetype, presence) %>% dplyr::rename(!!x := n)
  })
  df_presence <- Reduce(function(x, y) full_join(x, y, by=c("sitetype", "presence")), list_of_presence)
  df_presence
  
  
  df_presence %>% mutate(sim_cov = sim_cov) %>% select(sim_cov, everything()) %>%
    write.table(site_presence_fp, sep="\t", quote = F, row.names = F)
  
  
  
  ############### Part Two: abundance
  .prefix <- "ref_"
  list_of_refabun <- lapply(bt2_levels[2:length(bt2_levels)], function(x, prefix=.prefix) {
    maj_wide %>% select(ref_id, qryid, ref_pos, pos2, sitetype, `ani.100`, all_of(x)) %>% dplyr::rename(ani := !!x) %>% mutate(relabun = ani / `ani.100`) %>% dplyr::rename("{.prefix}{x}" := relabun) %>% select(-one_of(c("ani", "ani.100")))
  })
  df_relabun <- Reduce(function(x, y) full_join(x, y, by=c("ref_id", "qryid", "ref_pos", "pos2", "sitetype")), list_of_refabun)
  
  
  stopifnot(nrow(maj_wide) == nrow(df_relabun))
  maj_wide %<>% left_join(df_relabun, by=c("ref_id", "qryid", "ref_pos", "pos2", "sitetype"))
  
  
  #### findInterval: left close, right open [ ) <= no more
  #### findInterval: right close, left open ( ] <= 2020-11-10
  maj_long <- maj_wide %>% gather(ani, rel, colnames(df_relabun)[grepl(paste(.prefix, "ani", sep=""), colnames(df_relabun))])
  
  
  sites_abun <- maj_long %>% 
    #mutate(binold = findInterval(rel, seq(0, max(maj_long$rel)+0.1, 0.1))) %>% # <- first verson
    mutate(bin = findInterval(rel, seq(0, max(maj_long$rel)+0.1, 0.1), all.inside = T, left.open = T)) %>%
    group_by(ani, sitetype) %>% dplyr::count(bin) %>% ungroup()
  
  
  sites_abun %>% mutate(sim_cov = sim_cov) %>% select(sim_cov, everything()) %>%
    write.table(site_abun_fp, sep="\t", quote = F, row.names = F)
  
  
  
  #### Check if the relative readcounts is related to ani.100 abosulte read depths
  sites_corr <- maj_long %>% group_by(sitetype, ani) %>% summarise(corr = cor(ani.100, rel)) %>% ungroup() %>% spread(sitetype, corr, fill = 0)
  
  sites_corr %>% mutate(sim_cov = sim_cov) %>% select(sim_cov, everything()) %>%
    write.table(site_abuncorr_fp, sep="\t", quote = F, row.names = F)
  
  #### The above correlation shows the coverage is not the deciding factor but probably genomic content.
  
  
  
  ############### Part Three: welcome to the world of allele frequency
  maj_wide %<>% select(ref_id, ref_pos, ref_allele, qryid, pos2, sitetype, bt2_levels, sapply(bt2_levels, function(x) paste("depth_", x, sep=""), simplify = T, USE.NAMES = F))
  
  rc <- maj_wide %>% gather(bt2_db, readcounts, bt2_levels) %>% select(ref_id, ref_pos, qryid, pos2, sitetype, bt2_db, readcounts)
  dp <- maj_wide %>% gather(bt2_db, depth, sapply(bt2_levels, function(x) paste("depth_", x, sep=""), simplify = T, USE.NAMES = F)) %>% select(ref_id, ref_pos, qryid, pos2, sitetype, bt2_db, depth) %>%
    mutate(bt2_db = gsub("depth_", "", bt2_db))
  stopifnot(nrow(rc) == nrow(dp))
  
  af <- left_join(rc, dp, by=c("ref_id", "ref_pos", "qryid", "pos2", "sitetype", "bt2_db"))
  af %<>% mutate(allele_frequency = specify_decimal(readcounts/depth,3))
  af %<>% mutate(bin = findInterval(allele_frequency, seq(0.0, 1.0, 0.1), all.inside = T, left.open = T))
  
  #### When read counts and depths are zero, it doesn't mean the allele frequency is zero
  af_stats <- af %>% group_by(bt2_db, sitetype) %>% dplyr::count(bin) %>% ungroup() %>% mutate(bin = ifelse(is.na(bin), 11, as.numeric(bin)))
    
  af_stats %>% mutate(sim_cov = sim_cov) %>% select(sim_cov, everything()) %>%
    write.table(site_af_fp, sep="\t", quote = F, row.names = F)
  
}



####################### main functions #######################
if (length(args) > 0 ) {
  datadir <- args[1]
  species_under_investigation <- args[2]
  sim_cov <- args[3]
} else {
  datadir <- "~/Desktop/100107"
  species_under_investigation <- "100107"
  sim_cov <- 10
}


r_tempdir <- file.path(datadir, "6_rtemp")
r_outdir <- file.path(datadir, "7_rdata")


bt2_levels <- data.frame(bt2_db = list.files(file.path(datadir, "4_midas"))) %>% 
  separate(bt2_db, sep="\\.", into = c("_", "bt2_ani"), remove = F) %>% mutate(bt2_ani = as.numeric(bt2_ani)) %>% 
  mutate(bt2_db = as.factor(bt2_db)) %>% mutate(bt2_db = fct_reorder(bt2_db, bt2_ani, .desc=TRUE)) %>% .$bt2_db %>% levels()


gen_maj_analysis(r_tempdir, r_outdir, species_under_investigation, bt2_levels, sim_cov)
