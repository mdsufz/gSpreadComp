#!/usr/bin/env Rscript

#Test path

#setwd("//wsl.localhost/Ubuntu/home/kasmanas/mSpreadComp")


#### Load libs and inputs

library("data.table")
library("dplyr")
library("aricode") # NMI
library("optparse")
library("ggplot2")
library("tidyr")
library("gridExtra")


library("viridis")
library("pheatmap")
library("forcats")


# packages
library(stringr)

#library("patchwork")
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

#target_gene <- "Gene_id"
target_gene <- "Gene_class"

#### TEST INPUT ####

mags_data_df <- as.data.frame(data.table::fread("test_output/03_mspread_analysis_out/genome_quality_norm/genome_data_merged.csv"))

mags_data_df[, 'Domain'] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, 'Domain'])
mags_data_df[, 'Phylum'] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, 'Phylum'])
mags_data_df[, 'Class'] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, 'Class'])
mags_data_df[, 'Order'] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, 'Order'])
mags_data_df[, 'Family'] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, 'Family'])
mags_data_df[, 'Genus'] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, 'Genus'])
mags_data_df[, 'Species'] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, 'Species'])


#selected_lib <- data.table::fread("test_output/genome_quality_norm/selected_samples.csv")
gene_df <- data.table::fread("test_data/deeparg_df_format_mSpread.csv")
#norm_gene_prev_df <- data.table::fread("test_output/genome_quality_norm/gene_prevalence_per_library.csv")

tax_level <- "Phylum"
out.path <- "test_output"



# calculate NMI and chi-squared p-value per taxonomic category per Gene
target_genes <- unique(gene_df %>%
                         as.data.frame(.) %>%
                         select(all_of(target_gene)))

target_genes_t <- as.character(target_genes[, target_gene])

# set up data frame to collect results
nmi.results <- rbind(
  setNames(do.call("rbind",
                   replicate(length(target_genes_t),
                             as.data.frame(unique(cbind("Phylum", mags_data_df[, 'Phylum']))),
                             simplify = F)),
           c("tax.level", "taxon")),
  setNames(do.call("rbind",
                   replicate(length(target_genes_t),
                             as.data.frame(unique(cbind("Class", mags_data_df[, 'Class']))),
                             simplify = F)),
           c("tax.level", "taxon")),
  setNames(do.call("rbind",
                   replicate(length(target_genes_t),
                             as.data.frame(unique(cbind("Order", mags_data_df[, 'Order']))),
                             simplify = F)),
           c("tax.level", "taxon")),
  setNames(do.call("rbind",
                   replicate(length(target_genes_t),
                             as.data.frame(unique(cbind("Family", mags_data_df[, 'Family']))),
                             simplify = F)),
           c("tax.level", "taxon")),
  setNames(do.call("rbind",
                   replicate(length(target_genes_t),
                             as.data.frame(unique(cbind("Genus", mags_data_df[, 'Genus']))),
                             simplify = F)),
           c("tax.level", "taxon")),
  make.row.names = F)

# make sure that there are unique combinations of taxons and arg classes
nmi.results[, target_gene] <- c(
  rep(target_genes_t, each = length(unique(mags_data_df[, 'Phylum']))),
  rep(target_genes_t, each = length(unique(mags_data_df[, 'Class']))),
  rep(target_genes_t, each = length(unique(mags_data_df[, 'Order']))),
  rep(target_genes_t, each = length(unique(mags_data_df[, 'Family']))),
  rep(target_genes_t, each = length(unique(mags_data_df[, 'Genus']))))



# make sure that there are unique combinations of taxons and target genes
# tax_level <- as.character(nmi.results[, "tax.level"])
# taxon <- as.character(nmi.results[, "taxon"])
# 
# 
# nmi.results <- expand.grid(tax_level, taxon, target_genes_t)
# colnames(nmi.results) <- c("tax.level", "taxon", target_gene)

# NA values for calculated numbers
nmi.results$nmi <- NA
nmi.results$p.value <- NA
nmi.results$n <- NA

# combine all data
myData <- merge(mags_data_df,
                gene_df,
                by = "Genome",
                all.x = T)

myData <- as.data.frame(myData)

# calculate NMI and p-value
print("NMI and p-value calculation started!")
# Initialize the progress bar
pb <- txtProgressBar(min = 0, max = nrow(nmi.results), style = 3)

for (idx in 1:nrow(nmi.results)) {
  # extract relevant data
  #extr.data <- myData[which((myData[, colnames(myData) == nmi.results$tax.level[idx]] == nmi.results$taxon[idx]) & 
  #                            (myData[, target_gene] == nmi.results[idx, target_gene])), 
  #                    c(target_gene, colnames(myData)[which(colnames(myData) == nmi.results$tax.level[idx]) + 1])]
  
  
  test <-  myData[which((myData[, colnames(myData) == nmi.results$tax.level[idx]] == nmi.results$taxon[idx]) & 
                          (myData[, target_gene] == nmi.results[idx, target_gene])), ]
  
  
  
  # remove NA data
  extr.data <- extr.data[complete.cases(extr.data),]
  
  # if the both variables have different level values more than two
  if (nrow(extr.data) >= 1 & length(unique(extr.data[, 1])) > 1 & length(unique(extr.data[, 2]))> 1 ) {
    # compute NMI for all genes in the ARG class in all species found in the taxa
    nmi.results$nmi[idx] <- NMI(extr.data[, 1], extr.data[, 2])
    nmi.results$n[idx] <- nrow(extr.data)
    nmi.results$p.value[idx] <- chisq.test(extr.data[, 1], extr.data[, 2])$p.value
  }
  
  # Update the progress bar
  setTxtProgressBar(pb, idx)
}

# Close the progress bar
close(pb)

# extract the significant results
nmi.report <- nmi.results[which((nmi.results$nmi > 0.5) & (nmi.results$p.value < 0.05)), ]

# save tables
fwrite(nmi.results, "tables\\nmi_results.csv", sep = ",", row.names = F)
fwrite(nmi.report, "tables\\nmi_significant.csv", sep = ",", row.names = F)

# plot the significant NMI ARGs
p <- list()
for (i in 1:nrow(nmi.report)) {
  # extract data
  extr.data <- myData[which((myData[which(colnames(myData) == nmi.report$tax.level[i])] ==  nmi.report$taxon[i]) & 
                              (myData[, target_gene] == nmi.report$arg.class[i])), ]
  
  # dplyr needs a specific column name, can't use a variable
  extr.data$taxons <- extr.data[, which(colnames(extr.data) == nmi.report$tax.level[i]) + 1]
  
  # data to plot
  data.plot <- as.data.frame(
    extr.data %>% group_by(taxons, arg) %>% summarize(n = n_distinct(bin))
  )
  
  # plot
  #
  p[[i]] <- ggplot(data.plot, aes(x = taxons, y = n, fill = arg)) +
    geom_col(position="stack") +
    theme_bw() + 
    theme(axis.text = element_text(size = 12),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    labs(x = NULL, y = "Number of MAGs", title = paste0(nmi.report$taxon[i], ":", nmi.report$arg.class[i])) +
    coord_flip()
}
pdf(file = paste0("figures\\nmi_significant.pdf"),
    height = 16.5, width = 23.5)
do.call(grid.arrange, p)
dev.off()














