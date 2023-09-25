#!/usr/bin/env Rscript

#Format Victors or VFDB output for gSpreadComp

#### Load libs and inputs

library("optparse")
library("dplyr")
library("tidyr")

option_list = list(
  make_option(c("--vf"), type="character", default=NULL, 
              help="Path to the created blastx merged result file", metavar="character"),
  make_option(c("--vf_header"), type="character", default=NULL, 
              help="Path to the file with the headers from the selected vf database", metavar="character"),
  make_option(c("--e_value"), type="numeric", default=1e-50, 
              help="Evalue used in blastx", metavar="numeric"),
  make_option(c("--vf_prefix"), type="character", default="both", 
              help="Prefix to be used, based on the VF database name", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              help="Output directory path to save formatted tables", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out

vf.path <- opt$vf
vf.header.path <- opt$vf_header
e_value <- opt$e_value
vf_prefix <- opt$vf_prefix

##################################

#1) Load BlastX results for diet bin mapped to VFDB2022 and VictorsDB

vf_blastx_df <- data.table::fread(vf.path, sep = "\t")


vf_headers <- data.table::fread(vf.header.path, header = F, sep = NULL)
vf_headers$V1 <- sub(".", "", vf_headers$V1)

names(vf_blastx_df) <- c('qseqid',
                           'sseqid',
                           'pident',
                           'length',
                           'mismatch',
                           'gapopen',
                           'qstart',
                           'qend',
                           'sstart',
                           'send',
                           'evalue',
                           'bitscore',
                           'qcovs',
                           'qcovhsp',
                           'mag')


#2) Filter out BlastX results 

vf_blastx_df_filt <- vf_blastx_df %>%
  filter(evalue < e_value)

#3) Load Virulence Factors headers with annotation
####
vf_headers_split <- tidyr::separate(vf_headers,
                                      V1,
                                      into = c("header_id","header_class"),sep = " ",
                                      remove = FALSE, extra = "merge") %>%
  select(-V1) %>%
  filter(!duplicated(header_id))

#4) Merge Annotation to the filtered dataset

vf_merge <- merge.data.frame(x = vf_blastx_df_filt,
                 y = vf_headers_split,
                 by.x = "sseqid",
                 by.y = "header_id")

vf_merge$mag <- gsub(pattern = ".out",
                                   replacement = "",
                                   x = vf_merge$mag)

vf_merge$mag <- gsub(pattern = paste0(vf_prefix, "_"), 
                     replacement = "", 
                     x = vf_merge$mag)

#5) Save merge dfs to file 

vf_merge <- vf_merge %>%
  mutate(Sequence_id = paste0(mag, "_", qseqid)) %>%
  select(mag, Sequence_id, sseqid, header_class, evalue, bitscore) %>%
  rename(Genome = mag, VF_found = sseqid, VF_class = header_class)

# Save the modified dataframe to a file with a dynamic name based on vf_prefix
data.table::fwrite(vf_merge, file = paste0(out.path, "/", vf_prefix, "_format_gSpread.csv"))

#6) Count the number of Unique Virulence factors in each bin

vf_count <- vf_merge %>%
  select(VF_found, Genome, VF_class) %>%
  unique(.) %>%
  count(Genome) %>%
  as.data.frame(.)

#6.1) Save count tables to file
data.table::fwrite(vf_count, file = paste0(out.path, "/", vf_prefix, "_per_genome_unique_count.csv"))
