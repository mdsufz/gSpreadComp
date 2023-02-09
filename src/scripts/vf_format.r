#!/usr/bin/env Rscript

#Format Victors and VFDB output for mSpreadComp

#### Load libs and inputs

library("optparse")
library("dplyr")
library("tidyr")

option_list = list(
  make_option(c("--victors"), type="character", default=NULL, 
              help="Path to the created Victors merged result file", metavar="character"),
  make_option(c("--vfdb"), type="character", default=NULL, 
              help="Path to the created VFDB merged result file", metavar="character"),
  make_option(c("--vic_header"), type="character", default=NULL, 
              help="Path to the file with the headers from Victors database", metavar="character"),
  make_option(c("--vfdb_header"), type="character", default=NULL, 
              help="Path to the file with the headers from VFDB database", metavar="character"),			  
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formated tables", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out

vic.path <- opt$victors
vic.header.path <- opt$vic_header

vfdb.path <- opt$vfdb
vfdb.header.path <- opt$vfdb_header


##################################

#1) Load BlastX results for diet bin mapped to VFDB2022 and VictorsDB

vfdb_blastx_df <- data.table::fread(vfdb.path, sep = "\t")

victors_blastx_df <- data.table::fread(vic.path, sep = "\t")


vfdb_headers <- data.table::fread(vfdb.header.path, header = F, sep = NULL)
vfdb_headers$V1 <- sub(".", "", vfdb_headers$V1)

victors_headers <- data.table::fread(vic.header.path, header = F, sep = NULL)
victors_headers$V1 <- sub(".", "", victors_headers$V1)

names(vfdb_blastx_df) <- c('qseqid',
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

names(victors_blastx_df) <- c('qseqid',
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
vfdb_blastx_df_filt <- vfdb_blastx_df %>%
  filter(evalue < 1e-50)

victors_blastx_df_filt <- victors_blastx_df %>%
  filter(evalue < 1e-50)

#3) Load Virulence Factors headers with annotation

vfdb_headers_split <- tidyr::separate(vfdb_headers,
                               V1,
                               into = c("header_id","header_class"),sep = " ",
                               remove = FALSE,extra = "merge") %>%
  select(-V1) %>%
  filter(!duplicated(header_id))

####
victors_headers_split <- tidyr::separate(victors_headers,
                                      V1,
                                      into = c("header_id","header_class"),sep = " ",
                                      remove = FALSE, extra = "merge") %>%
  select(-V1) %>%
  filter(!duplicated(header_id))

#4) Merge Annotation to filtered dataset

vfdb_merge <- merge.data.frame(x = vfdb_blastx_df_filt,
                                  y = vfdb_headers_split,
                                  by.x = "sseqid",
                                  by.y = "header_id")

vfdb_merge$mag <- gsub(pattern = ".out",
                          replacement = "",
                          x = vfdb_merge$mag)
						  
						  
vfdb_merge$mag <- gsub(pattern = "vfdb_",
                                   replacement = "",
                                   x = vfdb_merge$mag)



####

victors_merge <- merge.data.frame(x = victors_blastx_df_filt,
                 y = victors_headers_split,
                 by.x = "sseqid",
                 by.y = "header_id")

victors_merge$mag <- gsub(pattern = ".out",
                                   replacement = "",
                                   x = victors_merge$mag)

victors_merge$mag <- gsub(pattern = "victors_",
                                   replacement = "",
                                   x = victors_merge$mag)

#5) Save merge dfs to file 

vfdb_merge <- vfdb_merge %>%
  mutate(Sequence_id = paste0(mag, "_", qseqid)) %>%
  select(mag, Sequence_id, sseqid, header_class, evalue, bitscore) %>%
  rename(Genome = mag, VFDB_VF_found = sseqid, VFDB_VF_class = header_class)

victors_merge <- victors_merge %>%
  mutate(Sequence_id = paste0(mag, "_", qseqid)) %>%
  select(mag, Sequence_id, sseqid, header_class, evalue, bitscore) %>%
  rename(Genome = mag, Victor_VF_found = sseqid, Victor_VF_class = header_class)


data.table::fwrite(victors_merge, file=paste0(out.path, "/victors_e_50_ann_format.csv"))
data.table::fwrite(vfdb_merge, file=paste0(out.path, "/vfdb_e_50_ann_format.csv"))

#6) Count number of Unique Virulence factors in each bin

vfdb_count <- vfdb_merge %>%
  select(VFDB_VF_found, Genome, VFDB_VF_class) %>%
  unique(.) %>%
  count(Genome) %>%
  as.data.frame(.)


victors_count <- victors_merge %>%
  select(Victor_VF_found, Genome, Victor_VF_class) %>%
  unique(.) %>%
  count(Genome) %>%
  as.data.frame(.)

#6.1) Save count tables to file
data.table::fwrite(vfdb_count, file=paste0(out.path, "/vfdb_per_genome_unique_count.csv"))
data.table::fwrite(victors_count, file=paste0(out.path, "/victors_per_genome_unique_count.csv"))






