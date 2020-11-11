## Chunyu Zhao 2020-10-21
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)


library(tidyverse)
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



datadir <- args[1]
species_under_investigation <- args[2]


read_coords <- function(coordsfile) {
  #' Read in the show-coords resulted file
  #'
  lines <- scan(coordsfile, 'a', sep='\n', skip=5, quiet=T)
  lines <- gsub("\\s+", "\t", str_trim(gsub("\\|", "", lines)))
  coords <- data.frame(do.call('rbind', strsplit(lines, "\t", fixed=T)), stringsAsFactors = FALSE) %>%
    set_colnames(c("s1", "e1", "s2", "e2", "len1", "len2", "pident", "len_ref", "len_query", "cov1", "cov2", "ref", "query"))
  coords[,1:11] <- sapply(coords[,1:11],as.numeric)
  coords %<>% select(ref, query, everything())
  coords
}



read_genome <- function(fasta_name) {
  #' Read in fasta file and return the genome seq in a data frame
  #' Return a dataframe with refid-refpos-refallele
  #'
  require("Biostrings")
  fastaFile <- readDNAStringSet(fasta_name)
  contig_ids = names(fastaFile)
  contig_seqs = paste(fastaFile)

  seqs <- bind_rows(lapply(1:length(contig_ids), function (x) tibble(refid = contig_ids[[x]], refallele = unlist(strsplit(contig_seqs[[x]], "", fixed=T))) %>%
                             rowid_to_column(var="refpos"))) %>% select(refid, refpos, refallele)
  seqs
}



read_nucmer_alnfile <- function(aln_file) {
  #' Read in the show-aln pairwise alignment file
  cmd_wc <- paste("wc -l ", aln_file, " | awk '{ print $1 }'", sep="")
  nrows <- system(command=cmd_wc, intern=TRUE)
  nrows <- as.numeric(nrows)


  alns <- list()
  bi = 1

  connection <- file(aln_file, open="r")
  counter <- 1
  while (counter < nrows) {
    line <- scan(file=connection, 'a', nlines=1, quiet=TRUE)
    counter <- counter + 1

    if (sum(grepl("Alignments", line)) & line[2] == "Alignments") {
      # this only be only ONCE
      ref = line[4]
      qry = line[6]
    }


    if (sum(grepl("BEGIN", line)) & line[2] == "BEGIN") {
      # start a new alignment block
      #print(paste("START alignment block", bi))
      curr_align <- data.frame(s1 = line[6], e1 = line[8], s2 = line[11], e2 = line[13], fmt1 = line[5], fmt2 = line[10], stringsAsFactors = F)
      curr_ref = ""
      curr_qry = ""
      curr_alnlen = as.numeric(line[8]) - as.numeric(line[6]) + 1

      while (TRUE) {
        if (nchar(curr_ref) - sum(str_detect(strsplit(curr_ref, "")[[1]], "\\.")) == curr_alnlen) {
          # save the old results
          curr_align$curr_ref <- curr_ref
          curr_align$curr_qry <- curr_qry
          alns[[bi]] <- curr_align
          #print(paste("END alignment block", bi))
          bi = bi + 1
          break
        }

        line <- scan(file=connection, 'a', nlines=3, quiet=TRUE)
        curr_ref <- paste(curr_ref, line[2], sep = "")
        curr_qry <- paste(curr_qry, line[4], sep = "")
        counter = counter + 3
      }
    }
  }
  line <- scan(file=connection, 'a', nlines=1, quiet=TRUE)
  stopifnot(str_count(line[1], "=") == 60)
  close(connection)

  bind_rows(alns) %>%
    mutate(s1 = as.numeric(as.character(s1)), s2 = as.numeric(as.character(s2)), e1 = as.numeric(as.character(e1)), e2 = as.numeric(as.character(e2)))
}



read_snps <- function(snps_file) {
  showsnps <- read_delim(snps_file, delim="\t", col_names = c("pos1", "sub_ref", "sub_qry", "pos2", "buff", "dist", "frm1", "frm2", "ref", "query"))
  
  ## When two genomes' ANI is 100, the following would happen.
  if (nrow(showsnps) == 0) {
    showsnps <- data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("pos1", "sub_ref", "sub_qry", "pos2", "buff", "dist", "frm1", "frm2", "ref", "query"))))
  }
  
  return(showsnps)
}



gen_tp_from_aln <- function(curr_snps, curr_block, refgenome, qrygenome) {
  #' Generate pairwise alignment for one alignment block from nucmer's show-coords, show-alns, and show-snps
  #'

  ## expand the whole genome pairwise alignment
  alnref = strsplit(curr_block$curr_ref, "")[[1]]
  alnqry = strsplit(curr_block$curr_qry, "")[[1]]
  aln_df = data.frame(refid = curr_block$ref, qryid = curr_block$query, subref = alnref, subqry = alnqry, stringsAsFactors = F) %>% mutate(rn = row_number())

  addpos1 <- aln_df %>% filter(subref != '.') %>% mutate(pos1 = curr_block$s1 : curr_block$e1)
  addpos2 <- aln_df %>% filter(subqry != '.') %>% mutate(pos2 = curr_block$s2 : curr_block$e2)

  aln_df <- left_join(aln_df, select(addpos1, rn, pos1), by="rn") %>% left_join(select(addpos2, rn, pos2), by="rn")
  aln_df %<>% mutate(pos1 = ifelse(is.na(pos1), -1, as.numeric(pos1))) %>% mutate(pos2 = ifelse(is.na(pos2), -1, as.numeric(pos2))) %>% select(-rn)


  ## the total site counts of homo-sites + snp-sites for both sequences are the same.
  ## len1 = homo + snps + alnqry=='-'
  ## len2 = homo + snps + alnref=='-'
  sitesref = curr_block$len1 - sum(alnqry == '.')
  sitesqry = curr_block$len2 - sum(alnref == '.')
  stopifnot(sitesref == sitesqry)


  aln_df %<>% left_join(refgenome, by=c("refid" = "refid", "pos1" = "refpos"))
  stopifnot(nrow(filter(aln_df, subref !='.' & toupper(subref) != refallele)) == 0)


  if (curr_block$fmt2 == "-1") {
    aln_df %<>% left_join(qrygenome, by=c("qryid" = "qryid", "pos2" = "qrypos")) %>% mutate(altallele = rc_allele) %>% select(-one_of(c("rc_allele", "qryallele")))
  } else {
    aln_df %<>% left_join(qrygenome, by=c("qryid" = "qryid", "pos2" = "qrypos")) %>% mutate(altallele = qryallele) %>% select(-one_of(c("rc_allele", "qryallele")))
  }
  stopifnot(nrow(filter(aln_df, subqry !='.' & toupper(subqry) != altallele)) == 0)


  aln_df %<>% filter(subref != '.' & subqry != '.') %>% select(-one_of(c("subref", "subqry"))) %>%
    mutate(sitetype = ifelse(refallele == altallele, "ref", "mism"))

  
  if (nrow(curr_snps) > 0) {
    curr_snps %<>% filter(sub_ref != "." & sub_qry != '.') %>% select(ref, query, pos1, pos2, sub_ref, sub_qry) %>% mutate(issnp = "yes")
    aln_df <- left_join(aln_df, curr_snps, by=c("refid" = "ref", "qryid" = "query", "pos1", "pos2", "refallele" = "sub_ref", "altallele" = "sub_qry"))
    aln_df %<>% mutate(issnp = ifelse(is.na(issnp), "no", issnp)) %>% mutate(sitetype = ifelse(issnp == "yes", "snp", sitetype)) %>% select(-issnp)
  }
  
  aln_df %>% dplyr::count(sitetype)

  return(aln_df)
}


visualize_coords <- function(showcoords, plot_fp) {
  showcoords %>%
    group_by(ref, query, len_ref, len_query) %>% summarise(total_aln_len_1 = sum(len1), total_aln_len_2 = sum(len2)) %>%
    mutate(aln_cov_1 = total_aln_len_1 / len_ref, aln_cov_2 = total_aln_len_2 / len_query) %>%
    ggplot(aes(x = ref, y = query)) +
    geom_tile(aes(fill = aln_cov_1)) +
    theme_bw() +
    scale_fill_material("red") +
    #scale_fill_gradient(low = "white", high = "red") +
    coord_equal() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste("showcoords")) + theme(plot.title = element_text(hjust = 0.5)) +
    ggsave(plot_fp, width = 6, height = 6, useDingbats = F)
}



aligned_sites_from_alns <- function(datadir, species_under_investigation) {
  #' Generate aligned sties between ref and qry from the whole genome alignments
  numcer_dir <- file.path(datadir, "5_nucmer")

  ## input files
  snps_file <- file.path(numcer_dir, paste(species_under_investigation, ".snps", sep=""))
  coords_file <- file.path(numcer_dir, paste(species_under_investigation, ".coords", sep=""))

  ## output files
  plot_fp <- file.path(numcer_dir, paste(species_under_investigation, ".heatmap_coords_refqry.pdf", sep=""))
  tp_sites_fp <- file.path(numcer_dir, paste(species_under_investigation, "_aligned_sites.tsv", sep=""))

   ## 2020-10-29: when two genomes ANI is 100, the showsnps would be empty.
  showsnps <- read_snps(snps_file)
  showcoords <- read_coords(coords_file)
  visualize_coords(showcoords, plot_fp)


  refgenome <- read_genome(file.path(datadir, "uhgg_rep.fna"))
  qrygenome <- read_genome(file.path(datadir, "genome_under_investigation.fna")) %>%
    set_colnames(c("qryid", "qrypos", "qryallele")) %>% mutate(rc_allele = chartr("ATGC","TACG",qryallele))

  ## Let's start
  aligned_sites <- list()
  for (r in unique(showcoords$ref)) {
    print(r)
    ralnDir <- file.path(numcer_dir, r)

    allqs = unique(filter(showcoords, ref == r) %>% .$ query)
    chrom_aln <- list()
    for (q in allqs) {
      rq_alnfile <- file.path(ralnDir, paste(q, ".aln", sep=""))
      stopifnot(file.exists(rq_alnfile))

      showalns <- read_nucmer_alnfile(rq_alnfile)

      blocks <- filter(showcoords, ref == r & query == q) %>% arrange(s1)
      blocks <- left_join(blocks, showalns, by=c("s1", "e1", "s2", "e2"))

      curr_snps <- showsnps %>% filter(ref == r & query == q)

      lob <- blocks %>%
        mutate(rn = row_number()) %>%
        group_by(rn) %>%
        group_map( ~ gen_tp_from_aln(curr_snps, .x, refgenome, qrygenome))

      chrom_aln[[q]] <- bind_rows(lob)
    }
    aligned_sites[[r]] <- bind_rows(chrom_aln)
  }

  write.table(bind_rows(aligned_sites), tp_sites_fp, sep="\t", quote=F, row.names = F)
}



aligned_sites_from_alns(datadir, species_under_investigation)
