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
# tax_level <- opt$spread_taxa

#target_gene <- "Gene_id"
target_gene <- "Gene_class"

#### TEST INPUT ####

mags_data_df <- data.table::fread("test_output/genome_quality_norm/genome_data_merged.csv")
selected_lib <- data.table::fread("test_output/genome_quality_norm/selected_samples.csv")
gene_df <- data.table::fread("test_data/deeparg_df_format_mSpread.csv")
norm_gene_prev_df <- data.table::fread("test_output/genome_quality_norm/gene_prevalence_per_library.csv")

tax_level <- "Phylum"

out.path <- "test_output"

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



### PARTICULAR FILTERING -> REMOVE AFTER ! ###

# only consider ARGs with 35% or higher percent identity
gene_df <- gene_df[gene_df$identity >= 35, ]

# load best bins to assign taxonomy to clusters
target_classes <- sort(c(unique(gene_df$Target)))

# remove underscore
mags_data_df[, tax_level] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, tax_level])


for (t in target_classes) {
  #We consider Species NA as the same if they have the same Genus in the same sample, otherwise is a different gOTU
  mags_data_fil <- mags_data_df %>%
    filter(Target == t) %>%
    mutate(unique_group = if_else(condition = is.na(Species), 
            true = paste(sep = "_", Genome, Genus), false = Species))
    
  
  # put all data together in a table
  gene.tax <- mags_data_fil %>%
    select(Genome, all_of(tax_level), unique_group)
  
  gene.tax <- merge(gene.tax,
                         gene_df,
                         by = "Genome")
  
  # summarize data
  #Count NA as unique species
  gene.tax.cnt <- 
    merge(
      as.data.frame(gene.tax %>% 
                      group_by(across(all_of(tax_level)), across(all_of(target_gene)), Target) %>% 
                      summarize(
                        gotu.gene.count = n_distinct(unique_group, na.rm = F))
      ),
      as.data.frame(gene.tax %>% 
                      group_by(across(all_of(tax_level))) %>% 
                      summarize(
                        gotu.count = n_distinct(unique_group, na.rm = F))
      ),
      by = tax_level
      
    )
  gene.tax.cnt$prevalence <- gene.tax.cnt$gotu.gene.count/gene.tax.cnt$gotu.count
  gene.tax.cnt <- gene.tax.cnt[!is.na(gene.tax.cnt[, target_gene]),]
  
  # Build matrix
  
  gene.tax.cnt.matrix <- gene.tax.cnt
  
  ac.tx <- matrix(0, 
                  nrow = length(unique(gene.tax.cnt.matrix[, target_gene])), 
                  ncol = length(unique(gene.tax.cnt.matrix[, tax_level])))
  rownames(ac.tx) <- unique(gene.tax.cnt.matrix[, target_gene])
  colnames(ac.tx) <- unique(gene.tax.cnt.matrix[, tax_level])
  
  for (i in 1:nrow(ac.tx)) {
    for (j in 1:ncol(ac.tx)) {
      prev <- gene.tax.cnt.matrix$prevalence[
        which((gene.tax.cnt.matrix[, target_gene] == rownames(ac.tx)[i]) & (gene.tax.cnt.matrix[, tax_level] == colnames(ac.tx)[j]))]
      if (length(prev) > 0) {
        ac.tx[i, j] <- prev
      }
    }
  }
  
  # color scale breaks
  breaksList = seq(0, 1, by = 0.01)
  # count number of mags in tax
  cnt.taxa <- as.data.frame(
    gene.tax %>%
      group_by(across(all_of(tax_level))) %>%
      summarize(unique_gotus = n_distinct(unique_group))
  )
  
  # Modify ordering of the clusters using clustering callback option
  # callback = function(hc, mat){
  #    sv = svd(t(mat))$v[,1]
  #    dend = reorder(as.dendrogram(hc), wts = sv)
  #    as.hclust(dend)
  #  }
  
  
  # only use if taxon has 5 gOTUs or more
  ac.tx.red <- ac.tx[, which(colnames(ac.tx) %in%
                               cnt.taxa[cnt.taxa$unique_gotus >= 5,
                                        tax_level])]
  
  
  # calculate weighted average prevalence (WAP)
  #We consider Species NA as the same if they have the same Genus in the same sample,
  #otherwise is a different gOTU
  
  unique_gotus <- length(unique(mags_data_fil$unique_group))
  wap <- rowSums(ac.tx.red %*% diag(cnt.taxa$unique_gotus[match(colnames(ac.tx.red),
                                                                cnt.taxa[, tax_level])])/unique_gotus)
  
  #Reorder rows (gene) from ac.tx.red
  wap <- sort(wap, decreasing = T)
  ac.tx.red <- ac.tx.red[names(wap),]
  
  pdf(file = paste0(out.path,"/prevalence_gene_", tax_level, "_", t, ".pdf"),
      height = 11.69, width = 8.27)
  
  print(
    
  pheatmap(ac.tx.red, 
           color = magma(100),
           breaks = breaksList,
           labels_col = paste0(colnames(ac.tx.red), " (", cnt.taxa$unique_gotus[match(colnames(ac.tx.red), cnt.taxa[, tax_level])], ")"),
           labels_row = paste0(rownames(ac.tx.red), " (", round(wap, 2), ")"),
           cutree_rows = 4,
           fontsize = 16,
           clustering_method = "ward.D2",
           cluster_cols = T,
           cluster_rows = F,
           main = t)
  
  )
  
  dev.off()
  
  #Save prev matrix and WAP vector to file

  write.csv(x = gene.tax.cnt, file = paste0(out.path, "/gene_prevalence_per_gOTU_", tax_level, "_", t, ".csv"), row.names = F)
  
  write.csv(x = wap, file = paste0(out.path, "/WAP_gene_", tax_level, "_", t, ".csv"), row.names = T)
  
  
}







