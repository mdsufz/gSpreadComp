#!/usr/bin/env Rscript

#Format GTDBtk output for gSpreadComp

#How to call it
#Rscript gtdbtk_format.r --gtdb gtdbtk_result.tsv -o output_folder

#Expect input for the same genomes:
#GTDB-tk

library("optparse")
library("dplyr")
library("tidyr")

option_list = list(
  make_option(c("--gtdb"), type="character", default=NULL, 
              help="GTDB-tk Taxa classification file path", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formatted tables", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out
gtdb_df <- data.table::fread(opt$gtdb)

#Format tables to select only relevant columns for input to gSpreadComp pipeline

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

#Save formatted tables
##GTDB
write.csv(x = gtdb_df_format, file = paste0(out.path, "/gtdb_df_format_gSpread.csv"), row.names = F)
