#!/usr/bin/env Rscript

#Format DeepARG output for mSpreadComp

#### Load libs and inputs
suppressPackageStartupMessages({
library("optparse")
library("dplyr")
library("tidyr")})

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to the created DeepARG results directories", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formatted tables", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out
wd.path <- opt$input

#1)Load deeparg results systematically

setwd(wd.path)

files <- dir(recursive=TRUE,
             full.names=TRUE,
             pattern="\\.mapping.ARG$")


#Start dataframe with one example

deeparg_bind_results <- read.delim(file = files[1])
bin_name <- sub(".fa.*", "", files[1])

deeparg_bind_results$bin_id <- bin_name

no_deeparg_found_bins <- c()

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
    
    deeparg_bind_results <- rbind.data.frame(deeparg_bind_results,
                                             args_found_df)
  } else {
    no_deeparg_found_bins <- c(no_deeparg_found_bins, bin_name)
    
  }
  
  # Flush the output to show the progress bar
  flush.console()
  
}

# End the progress bar with a new line
cat("\n")

#Format
deeparg_bind_results_format <- deeparg_bind_results %>%
  mutate(Genome = gsub("./", "", bin_id)) %>%
  mutate(Gene_sequence_location = paste0(Genome, "_", read_id)) %>%
  select(Genome, X.ARG, predicted_ARG.class, Gene_sequence_location, probability, identity, alignment.evalue) %>%
  rename(Gene_id = X.ARG, Gene_class = predicted_ARG.class)

#2) Reorder deeparg_bind_results_format and save to file

write.csv(deeparg_bind_results_format, file = paste0(out.path, "/deeparg_df_format_gSpread.csv"),
          row.names = F)

write.csv(deeparg_bind_results, file = paste0(out.path, "/deeparg_df_combined_raw.csv"),
          row.names = F)

write.csv(no_deeparg_found_bins, file = paste0(out.path, "/genomes_with_no_found_deeparg.csv"),
          row.names = F)



