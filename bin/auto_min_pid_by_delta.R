## Chunyu Zhao 2020-10-21
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)


library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(gridExtra)
library(magrittr)


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


auto_min_pid_by_delta <- function(coordsfile, plot_fp, max_pid_delta = 0.01, option = 2) {
  ## Parameter defining the maximum identity gap between identity of each aligned block and whole-genome ANI, all alignments with identity less than ANI * (1 - delta) will be purged, [0, 1]

  # coords also include indels, need to filter them out.
  showcoords <- read_coords(coordsfile)

  avg_pid = mean(showcoords$pident)

  min_pid_by_delta = avg_pid * (1-max_pid_delta) ## 0.99 average
  min_pid_by_95quantile = quantile(showcoords$pident, 0.05, 0.95)[[1]]

  p1 <- showcoords %>% group_by(ref, query) %>% summarise(avg_pident = mean(pident), total_len = sum(len1)) %>%
    ggplot(aes(x = total_len, y = avg_pident)) + geom_point(aes(size = total_len), shape = 1) +
    geom_hline(yintercept = min_pid_by_delta, color = "pink") +
    geom_hline(yintercept = min_pid_by_95quantile, color = "red", linetype="longdash") +
    scale_y_log10()
  p1 <- p1 + geom_point(data=showcoords, aes(x = len1, y = pident), color = "lightblue", shape = 20)

  p2 <- showcoords %>% ggplot(aes(x = pident)) + geom_histogram() +
    geom_vline(xintercept = min_pid_by_delta, color = "pink") +
    geom_vline(xintercept = min_pid_by_95quantile, color = "red", linetype="longdash")


  ggsave(plot_fp, arrangeGrob(p1, p2, ncol = 2), width = 10, height = 5)

  if (option == 2) {
    return(round(min_pid_by_95quantile))
  } else {
    return(round(min_pid_by_delta))
  }

}


min_coords_pid <- auto_min_pid_by_delta(args[1], args[2])
cat(min_coords_pid)
