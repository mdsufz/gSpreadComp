#!/usr/bin/env Rscript

#Test path

setwd("//wsl.localhost/Ubuntu/home/kasmanas/mSpreadComp")

#Select Samples based on number of MAGs recovered
#Describe MAGs and Save figures

#### Load libs and inputs

library("optparse")
library("dplyr")
library("tidyr")
library("ggplot2")
library("data.table")

library("patchwork")
library("viridis")
library("gridExtra")
library("ggforce")


option_list = list(
  make_option(c("--mags_data"), type="character", default=NULL, 
              help="Merged Genome information", metavar="character"),
  make_option(c("--selected_lib"), type="character", default=NULL, 
              help="Libraaries selected based on number of MAGs", metavar="character"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="gene merged result for all genomes file path", metavar="character"),
  make_option(c("--norm_gene_prev"), type="character", default=NULL, 
              help="Normalized gene prevalence per library", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formated files", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# out.path <- opt$out
# mags_data_df <- data.table::fread(opt$mags_data)
# selected_lib <- data.table::fread(opt$selected_lib)
# gene_df <- data.table::fread(opt$gene)
# norm_gene_prev_df <- data.table::fread(opt$norm_gene_prev)

#### TEST INPUT ####

mags_data_df <- data.table::fread("test_output/genome_quality_norm/genome_data_merged.csv")
selected_lib <- data.table::fread("test_output/genome_quality_norm/selected_samples.csv")
gene_df <- data.table::fread("test_data/deeparg_df_format_mSpread.csv")
norm_gene_prev_df <- data.table::fread("test_output/genome_quality_norm/gene_prevalence_per_library.csv")


#### Process initial load data ####

selected_lib <- as.character(selected_lib$x)

mags_data_df <- mags_data_df %>%
  filter(Library %in% selected_lib)

gene_df <- gene_df %>%
  merge.data.frame(.,
                   mags_data_df[,c("Library", "Genome")],
                   by = "Genome") %>%
  filter(Library %in% selected_lib)

norm_gene_prev_df <- norm_gene_prev_df %>%
  filter(Library %in% selected_lib)


















