## Chunyu Zhao 2020-12-09
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)


library(tidyverse)
library(magrittr)

get_major_allele <- function(count_a, count_c, count_g, count_t) {
  acgt_counts <- c(count_a, count_c, count_g, count_t)
  names(acgt_counts) <- c("A", "C", "G", "T")

  alleles_above_cutoff <- sort(acgt_counts[acgt_counts >= min_allele_depth], decreasing = TRUE)
  allele_counts <- length(alleles_above_cutoff)

  if (allele_counts == 0) {
    # no genotyped
    return(paste(c(NA, 0), collapse = ","))
  }

  alleles_above_cutoff <- alleles_above_cutoff[1:min(allele_counts, 2)]

  major_allele <- names(alleles_above_cutoff[1])
  minor_allele <- names(alleles_above_cutoff[min(allele_counts, 2)]) # for fixed sites, same as major allele

  if (allele_counts > 1 && as.numeric(alleles_above_cutoff[major_allele]) == as.numeric(alleles_above_cutoff[minor_allele])) {
    # undetermined
    major_called <- c(major_allele, -1)
  } else {
    major_called <- c(major_allele, as.numeric(alleles_above_cutoff[major_allele]))
  }
  return(paste(major_called, collapse = ","))
}


gen_sites_summary <- function(pileup_fp, min_allele_depth) {

  if (! file.exists(pileup_fp)) {
    system(paste("lz4 -d -c", paste(pileup_fp, ".lz4", sep=""), " > ", pileup_fp, ep=""))
  }

  pileup <- read_delim(pileup_fp, delim = "\t", col_types = cols())

  pileup <- left_join(true_sites, pileup, by=c("refid" = "ref_id", "pos1" = "ref_pos", "refallele" = "ref_allele"))
  sites_summary <- true_sites %>% dplyr::count(sitetype) %>% dplyr::rename(total_sites_counts = n)


  un_aligned <- pileup %>% filter(is.na(depth)) %>% dplyr::count(sitetype) %>% dplyr::rename(sites_unaligned = n)
  sites_summary <- left_join(sites_summary, un_aligned, by=c("sitetype"))


  aligned <- pileup %>% filter(!is.na(depth)) %>%
    rowwise() %>%
    mutate(major_allele_counts = get_major_allele(count_a, count_c, count_g, count_t))
  aligned %<>% separate(major_allele_counts, into=c("major_allele", "read_counts"), sep=",")


  un_called <- aligned %>% filter(read_counts <= 0)
  sites_summary <- left_join(sites_summary, un_called %>% dplyr::count(sitetype) %>% dplyr::rename(sites_uncalled = n), by=c("sitetype"))


  called <- aligned %>% filter(read_counts > 0) %>% mutate(called_label = ifelse(altallele == major_allele, "called_correct", "called_incorrect"))
  sites_summary <- left_join(sites_summary, called %>% dplyr::count(sitetype, called_label) %>% spread(called_label, n, fill = 0), by=c("sitetype"))

  return(sites_summary)
}




####################### main functions #######################
## Input
datadir <- args[1]
species_under_investigation <- args[2]
true_sites_fp <- args[3]
sim_cov <- args[4]
do_filter <- args[5]


## Output
if (do_filter == 1) {
  midas_dir <- file.path(datadir, "4_midas")
  min_allele_depth = 2
  r_outdir <- file.path(datadir, "6_rdata")
  species_outdir <- file.path(r_outdir,species_under_investigation)
  o_summary_fp <- file.path(species_outdir, paste("sites_summary_", sim_cov, "X.tsv", sep=""))
} else {
  midas_dir <- file.path(datadir, "4_midas_nofilter")
  min_allele_depth = 1
  r_outdir <- file.path(datadir, "6_rdata_nofilter")
  species_outdir <- file.path(r_outdir,species_under_investigation)
  o_summary_fp <- file.path(species_outdir, paste("sites_summary_", sim_cov, "X.tsv", sep=""))
}


## Prepare
bt2_levels <- data.frame(bt2_db = list.files(file.path(datadir, "4_midas"))) %>%
  separate(bt2_db, sep="\\.", into = c("_", "bt2_ani"), remove = F) %>% mutate(bt2_ani = as.numeric(bt2_ani)) %>%
  mutate(bt2_db = as.factor(bt2_db)) %>% mutate(bt2_db = fct_reorder(bt2_db, bt2_ani, .desc=TRUE)) %>% .$bt2_db %>% levels()

true_sites <- read_delim(true_sites_fp, delim = "\t", col_types = cols())


stats_list <- list()
for (bt2_db in bt2_levels) {
  sample_name <- paste("art_", bt2_db, "_", sim_cov, "X", sep="")
  pileup_fp <- file.path(midas_dir, bt2_db, "out", sample_name, "snps", paste(species_under_investigation, ".snps.tsv", sep=""))
  stats_list[[as.character(bt2_db)]] <- gen_sites_summary(pileup_fp, min_allele_depth)
}

bind_rows(stats_list, .id = "bt2_db") %>% select(bt2_db, everything()) %>% write.table(o_summary_fp, sep="\t", quote = F, row.names = F)
