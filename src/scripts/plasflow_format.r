#!/usr/bin/env Rscript

#Format plasflow output for gSpreadComp

#### Load libs and inputs
suppressPackageStartupMessages({
library("optparse")
library("dplyr")
library("tidyr")})

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to the created Plasflow results directories", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formatted tables", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out
wd.path <- opt$input

#1)Load Plasflow results systematically

setwd(wd.path)

files <- dir(recursive=TRUE,
             full.names=TRUE,
             pattern="\\_plasflow_out.tsv$")


#Start dataframe with one example

plasflow_bind_results <- read.delim(file = files[1])
bin_name <- sub(".fa.*", "", files[1])

plasflow_bind_results$bin_id <- bin_name

plasflow_bind_results <- plasflow_bind_results %>%
  select(bin_id, label ,contig_name, contig_length, id)

no_plasflow_found_bins <- c()

total_iterations <- length(files)

for (i in 2:length(files)) {


  # Compute the percentage of completion
  percent_complete <- (i / total_iterations) * 100
  # Print the percentage of completion
  cat(sprintf("\rProgress: %.0f%%", percent_complete))
  
  args_found_df <- read.delim(file = files[i])
  
  bin_name <- sub(".*_dir/", "", files[i])
  bin_name <- sub(".fa.*", "", bin_name)
  
  if (nrow(args_found_df) > 0) {
    
    args_found_df$bin_id <- bin_name
    
    args_found_df <- args_found_df %>%
      select(bin_id, label, contig_name, contig_length, id)
    
    plasflow_bind_results <- rbind.data.frame(plasflow_bind_results,
                                             args_found_df)
  } else {
    no_plasflow_found_bins <- c(no_plasflow_found_bins, bin_name)
    
  }
  
  # Flush the output to show the progress bar
  flush.console()
  
}

# End the progress bar with a new line
cat("\n")

#Format
plasflow_bind_results_format <- plasflow_bind_results %>%
  separate(label, into = c("sequence_type", "tax")) %>%
  mutate(Genome = gsub("./", "", bin_id)) %>%
  mutate(Sequence_id = paste0(Genome, "_", contig_name)) %>%
  select(Genome, sequence_type, Sequence_id)

#2) Reorder plasflow_bind_results_format and save to file

write.csv(plasflow_bind_results_format, file = paste0(out.path, "/plasflow_combined_format_gSpread.csv"),
          row.names = F)

write.csv(no_plasflow_found_bins, file = paste0(out.path, "/genomes_with_no_found_plasflow.csv"),
          row.names = F)
