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


maj_pass_two <- function(r_tempdir, r_outdir, species_under_investigation, bt2_levels) {
  #' MERGE maj_wide => calculate stats
  
  species_outdir <- file.path(r_outdir, species_under_investigation)
  
  
  det_power_fp <- file.path(r_outdir, paste(species_under_investigation, "_detection_power.tsv", sep=""))
  site_presence_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_presence.tsv", sep=""))
  site_abuncorr_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_abuncorr.tsv", sep=""))
  site_abun_fp <- file.path(r_outdir, paste(species_under_investigation, "_sites_abundance.tsv", sep=""))
  
  
  list_of_detection_power <- list()
  site_presence_list <- list()
  site_abuncorr_list <- list()
  site_abun_list <- list()

  
  for (sim_cov in 50:1) {
    print(sim_cov)
    
    cal_maj_fp <- file.path(r_tempdir, paste("maj.wide_", sim_cov, "X.tsv", sep=""))
    stopifnot(file.exists(cal_maj_fp))
    
    stats_fp <- file.path(species_outdir, paste("maj.stats_", sim_cov, "X.tsv", sep=""))
    stopifnot(file.exists(stats_fp))
    
    
    hist_abun_fp <- file.path(species_outdir, paste("hist_detected_sites_relabun_", sim_cov, "X.pdf", sep=""))
    
    
    maj_wide <- read_delim(cal_maj_fp, delim = "\t", col_types = cols())
    list_of_detection_power[[as.character(sim_cov)]] <- read_delim(stats_fp, delim = "\t", col_types = cols())
    
    ############### What to do with sites not in ani.100 but in other bt2_anis? What's the reason for this observation?
    #### does ani.100 miss out sites (FN) or the other bt2_indexes reports FP sites
    extra_sites <- maj_wide %>% filter(ani.100 == 0)
    extra_sites %<>% gather(ani, readcounts, bt2_levels) %>% filter(readcounts > 0) %>% group_by(refid, qryid, pos1, pos2, sitetype) %>% summarise(nz = n(), totalRC = sum(readcounts)) %>% ungroup() 
    extra_sites %>% group_by(sitetype) %>% summarise(counts = n(), mean_nz = mean(nz), med_nz = median(nz), med_rc = median(totalRC), mean_rc = mean(totalRC))
    extra_sites %>% dplyr::count(sitetype)
    ## WHAT TO DO with the extra_sites?
    
    
    maj_wide %<>% filter(ani.100 > 0) #<- need to ask Katie about this one
    
    
    ############### Part One: presence-absence: among the detected-true-sites by at least one bt2_index, what is the number of sites that have been missed/called across bt2_anis.
    list_of_presence <- lapply(bt2_levels, function(x) {
      maj_wide %>% select(sitetype, all_of(x)) %>% dplyr::rename(ani := !!x) %>% mutate(presence = ifelse(ani > 0, "called", "missed")) %>% dplyr::count(sitetype, presence) %>% dplyr::rename(!!x := n)
    })
    df_presence <- Reduce(function(x, y) full_join(x, y, by=c("sitetype", "presence")), list_of_presence)
    df_presence
    
    site_presence_list[[as.character(sim_cov)]] <- df_presence
    
    
    
    ############### Part Two: abundance
    .prefix <- "ref_"
    list_of_refabun <- lapply(bt2_levels[2:length(bt2_levels)], function(x, prefix=.prefix) {
      maj_wide %>% select(refid, qryid, pos1, pos2, sitetype, `ani.100`, all_of(x)) %>% dplyr::rename(ani := !!x) %>% mutate(relabun = ani / `ani.100`) %>% dplyr::rename("{.prefix}{x}" := relabun) %>% select(-one_of(c("ani", "ani.100")))
    })
    df_relabun <- Reduce(function(x, y) full_join(x, y, by=c("refid", "qryid", "pos1", "pos2", "sitetype")), list_of_refabun)
    
    stopifnot(nrow(maj_wide) == nrow(df_relabun))
    maj_wide %<>% left_join(df_relabun, by=c("refid", "qryid", "pos1", "pos2", "sitetype"))
    
    
    #### Histogram of relative read counts compared to ani.100, for varying sitetypes
    #### K: this is a histogram of *sites*, and what's the y-axis.
    if (FALSE) {
      maj_wide %>% 
        gather(ani, rel, colnames(df_relabun)[grepl(.prefix, colnames(df_relabun))]) %>%
        ggplot(aes( x = rel)) + 
        geom_histogram() + 
        scale_y_log10() +
        theme_bw() + 
        facet_wrap(ani~sitetype, scales = "free", ncol = 3) + 
        geom_vline(xintercept = 1, color = "red") + 
        ggtitle(paste("For sites detected in ani.100, the distributuion of relative read counts")) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        ggsave(hist_abun_fp, width = 10, height = 12, useDingbats = F)
    }
    
    
    #### findInterval: left close, right open [ )
    #findInterval(c(1, 0, -0.01, 1.11, 0.92, 0.83), seq(0, 1.1, 0.1), all.inside = T)
    
    maj_long <- maj_wide %>% gather(ani, rel, colnames(df_relabun)[grepl(.prefix, colnames(df_relabun))])
    
    sites_abun <- maj_long %>% 
      mutate(bin = findInterval(rel, seq(0, max(maj_long$rel)+0.1, 0.1))) %>%
      group_by(ani, sitetype) %>% dplyr::count(bin) %>% ungroup()
    site_abun_list[[as.character(sim_cov)]] <- sites_abun
    
    #### Check if the relative readcounts is related to ani.100 abosulte read depths
    sites_corr <- maj_long %>% group_by(sitetype, ani) %>% summarise(corr = cor(ani.100, rel)) %>% ungroup() %>% spread(sitetype, corr, fill = 0)
    site_abuncorr_list[[as.character(sim_cov)]] <- sites_corr
    #### The above correlation shows the coverage is not the deciding factor but probably genomic content.
    
    
    #### TODO: maj_wides are too big to save, let me think about it ...
    #maj_sites_list[[as.character(sim_cov)]] <- maj_wide
  }
  
  
  bind_rows(list_of_detection_power, .id="sim_cov") %>% 
    write.table(det_power_fp, sep="\t", quote = F, row.names = F)
  
  
  bind_rows(site_abuncorr_list, .id = "sim_cov") %>%
    write.table(site_abuncorr_fp, sep="\t", quote = F, row.names = F)
  
  bind_rows(site_abun_list, .id = "sim_cov") %>%
    write.table(site_abun_fp, sep="\t", quote = F, row.names = F)
  
  
  bind_rows(site_presence_list, .id = "sim_cov") %>%
    write.table(site_presence_fp, sep="\t", quote = F, row.names = F)
  
}


####################### main functions #######################
datadir <- args[1]
species_under_investigation <- args[2]


r_tempdir <- file.path(datadir, "6_rtemp")
r_outdir <- file.path(datadir, "7_rdata")

bt2_levels <- data.frame(bt2_db = list.files(file.path(datadir, "4_midas"))) %>% 
  separate(bt2_db, sep="\\.", into = c("_", "bt2_ani"), remove = F) %>% mutate(bt2_ani = as.numeric(bt2_ani)) %>% 
  mutate(bt2_db = as.factor(bt2_db)) %>% mutate(bt2_db = fct_reorder(bt2_db, bt2_ani, .desc=TRUE)) %>% .$bt2_db %>% levels()


maj_pass_two(r_tempdir, r_outdir, species_under_investigation, bt2_levels)
