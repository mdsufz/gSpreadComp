#!/usr/bin/env Rscript

#Format GTDBtk output for gSpreadComp

#How to call it
#Rscript gtdbtk_format.r --gtdb gtdbtk_all_mags.tsv --checkm checkm_all_mags.tsv --deeparg 01_deeparg_output_combined_edit.csv --meta 01_diet_metadata_std_all_edit.csv -o 02_format_input_dfs/

#Preprocess script for mSpreadComp

#Expect input for the same genomes:
#GTDB-tk
#CheckM
#DeepARG complete
#Sample Metadata

library("optparse")
library("dplyr")
library("tidyr")

option_list = list(
  make_option(c("--gtdb"), type="character", default=NULL, 
              help="GTDB-tk Taxa classification file path", metavar="character"),
  make_option(c("--checkm"), type="character", default=NULL, 
              help="CheckM quality result file path", metavar="character"),
  make_option(c("--deeparg"), type="character", default=NULL, 
              help="deeparg merged result for all genomes file path", metavar="character"),
  make_option(c("--meta"), type="character", default=NULL, 
              help="Library metadata file path", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formated tables", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out
gtdb_df <- data.table::fread(opt$gtdb)
checkm_df <- data.table::fread(opt$checkm)
deeparg_df <- data.table::fread(opt$deeparg)
meta_df <- data.table::fread(opt$meta)

#Testing with samples inside
# out.path <- opt$out
# 
# gtdb_df <- data.table::fread("01_input_data/gtdbtk_all_mags.tsv")
# checkm_df <- data.table::fread("01_input_data/checkm_all_mags.tsv")
# deeparg_df <- data.table::fread("01_input_data/01_deeparg_output_combined_edit.csv")
# meta_df <- data.table::fread("01_input_data/01_diet_metadata_std_all_edit.csv")


#Format tables to select only relevant columns form input to mSpreadComp pipeline

### GTDB
gtdb_df_format <- gtdb_df %>%
  select(user_genome, classification) %>%
  rename(Genome = user_genome)

# remove "s__" etc. (e.g. "p__Firmicutes" to "Firmicutes")
gtdb_df_format$Domain <- stringr::str_extract(gtdb_df_format$classification, "(?<=d__).*?(?=;)")
gtdb_df_format$Phylum <- stringr::str_extract(gtdb_df_format$classification, "(?<=p__).*?(?=;)")
gtdb_df_format$Class <- stringr::str_extract(gtdb_df_format$classification, "(?<=c__).*?(?=;)")
gtdb_df_format$Order <- stringr::str_extract(gtdb_df_format$classification, "(?<=o__).*?(?=;)")
gtdb_df_format$Family <- stringr::str_extract(gtdb_df_format$classification, "(?<=f__).*?(?=;)")
gtdb_df_format$Genus <- stringr::str_extract(gtdb_df_format$classification, "(?<=g__).*?(?=;)")
gtdb_df_format$Species <- stringr::str_extract(gtdb_df_format$classification, "(?<=s__).*")

# replace empty with NA
gtdb_df_format$Phylum[gtdb_df_format$Phylum == ""] <- NA
gtdb_df_format$Class[gtdb_df_format$Class == ""] <- NA
gtdb_df_format$Order[gtdb_df_format$Order == ""] <- NA
gtdb_df_format$Family[gtdb_df_format$Family == ""] <- NA
gtdb_df_format$Genus[gtdb_df_format$Genus == ""] <- NA
gtdb_df_format$Species[gtdb_df_format$Species == ""] <- NA

#Remove classification col
gtdb_df_format <- gtdb_df_format %>% select(-classification)
