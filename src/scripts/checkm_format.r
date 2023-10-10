#!/usr/bin/env Rscript

#How to call it
#Rscript checkm_format.r --checkm checkm_all_mags.tsv -o output_folder

#Preprocess script for gSpreadComp

#Expect input for the same genomes:
#CheckM
suppressPackageStartupMessages({
library("optparse")
library("dplyr")
library("tidyr")})

option_list = list(
  make_option(c("--checkm"), type="character", default=NULL, 
              help="CheckM quality result file path", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formatted tables", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out
checkm_df <- data.table::fread(opt$checkm)

### CheckM
checkm_df_format <- checkm_df %>%
  select(`Bin Id`, Completeness, Contamination, `Strain heterogeneity`) %>%
  rename(Genome = `Bin Id`)

##CheckM
write.csv(x = checkm_df_format, file = paste0(out.path, "/checkm_df_format_gSpread.csv"), row.names = F)
