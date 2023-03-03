#!/usr/bin/env Rscript

#Test path

#setwd("//wsl.localhost/Ubuntu/home/kasmanas/mSpreadComp")

#Select Samples based on number of MAGs recovered
#Describe MAGs and Save figures

#### Load libs and inputs

library("optparse")
library("dplyr")
library("tidyr")
library("ggplot2")
library("data.table")
library("viridis")
library("pheatmap")
library("forcats")

#library("patchwork")
#library("gridExtra")
#library("ggforce")

#library(stringr)


option_list = list(
  make_option(c("--mags_data"), type="character", default=NULL, 
              help="Merged Genome information", metavar="character"),
  make_option(c("--selected_lib"), type="character", default=NULL, 
              help="Libraries selected based on number of MAGs", metavar="character"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="gene merged result for all genomes", metavar="character"),
  make_option(c("--norm_gene_prev"), type="character", default=NULL, 
              help="Normalized gene prevalence per library", metavar="character"),
  make_option(c("--spread_taxa"), type="character", default=NULL, 
              help="Taxa level to use calculate gene spread [default: Phylum]", metavar="character"),
  make_option(c("--target_gene_col"), type="character", default=NULL, 
              help="Name of the Target Gene column [default: Gene_id]", metavar="character"),
  make_option(c("--patho_db"), type="character", default=NULL, 
              help="Dataset with the human bacterial pathogens taxonomic classification", metavar="character"),
  make_option(c("--vf"), type="character", default=NULL, 
              help="Virulence factors blastx alignment formated dataset", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formated files", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out
mags_data_df <- data.table::fread(opt$mags_data)
selected_lib <- data.table::fread(opt$selected_lib)
gene_df <- data.table::fread(opt$gene)
norm_gene_prev_df <- data.table::fread(opt$norm_gene_prev)
tax_level <- opt$spread_taxa
target_gene <- opt$target_gene_col
patho_db <-data.table::fread(opt$patho_db)
vf <- data.table::fread(opt$vf)

#target_gene <- "Gene_id"
target_gene <- "Gene_class"

#### TEST INPUT ####

mags_data_df <- data.table::fread("test_out/03_mspread_analysis_out/genome_quality_norm/genome_data_merged.csv")
selected_lib <- data.table::fread("test_out/03_mspread_analysis_out/genome_quality_norm/selected_samples.csv")
gene_df <- data.table::fread("test_data/deeparg_df_format_mSpread.csv")
norm_gene_prev_df <- data.table::fread("test_out/03_mspread_analysis_out/genome_quality_norm/gene_prevalence_per_library.csv")

vf_df <- data.table::fread("test_data/victors_df_format_mSpread.csv")
patho_db <- data.table::fread("installation/dependencies/patho_ncbi_20230222_taxa.csv")

#tax_level <- "Phylum"
#out.path <- "test_output"

#### Process initial load data ####

selected_lib <- as.character(selected_lib$x)

mags_data_df <- mags_data_df %>%
  filter(Library %in% selected_lib) %>%
  as.data.frame(.)

gene_df <- gene_df %>%
  merge.data.frame(.,
                   mags_data_df[,c("Library", "Genome", "Target")],
                   by = "Genome") %>%
  filter(Library %in% selected_lib)

norm_gene_prev_df <- norm_gene_prev_df %>%
  filter(Library %in% selected_lib)

vf_df_t <- vf_df %>%
  filter(Genome %in% as.vector(c(unique(mags_data_df$Genome)))) %>%
  as.data.frame(.)


### PARTICULAR FILTERING -> REMOVE AFTER ! ###

# only consider ARGs with 35% or higher percent identity
gene_df <- gene_df[gene_df$identity >= 35, ]

# load best bins to assign taxonomy to clusters
target_classes <- sort(c(unique(gene_df$Target)))

# remove underscore
mags_data_df[, tax_level] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, tax_level])
