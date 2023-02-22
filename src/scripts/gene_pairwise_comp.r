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

#Set Default parameters
#target_gene <- "Gene_id"
target_gene <- "Gene_class"
tax_level <- "Phylum"


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
              help="Name of the Target Gene column to perform the pairwise comparisons [default: Gene_id]", metavar="character"),
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


#### TEST INPUT ####

#mags_data_df <- data.table::fread("test_output/03_mspread_analysis_gene_class_out/genome_quality_norm/genome_data_merged.csv")
#selected_lib <- data.table::fread("test_output/03_mspread_analysis_gene_class_out/genome_quality_norm/selected_samples.csv")
#gene_df <- data.table::fread("test_data/deeparg_df_format_mSpread.csv")
#norm_gene_prev_df <- data.table::fread("test_output/03_mspread_analysis_gene_class_out/genome_quality_norm/gene_prevalence_per_library.csv")
#out.path <- "test_output"

#### Process initial load data ####

selected_lib <- as.character(selected_lib$x)

mags_data_df <- mags_data_df %>%
  filter(Library %in% selected_lib) %>%
  as.data.frame(.)

#REMOVE AFTER DONE!
# remove underscore
mags_data_df[, tax_level] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, tax_level])

gene_df <- gene_df %>%
  merge.data.frame(.,
                   mags_data_df[,c("Library", "Genome", "Target")],
                   by = "Genome") %>%
  filter(Library %in% selected_lib)

norm_gene_prev_df <- norm_gene_prev_df %>%
  filter(Library %in% selected_lib)

target_gene_sum_prev <- norm_gene_prev_df %>%
  group_by(across(all_of(target_gene))) %>%
  summarise(prev_sum = sum(gene.genome.prev)) %>%
  tibble::deframe() %>%
  sort()

### PARTICULAR FILTERING -> REMOVE AFTER ! ###

# only consider ARGs with 35% or higher percent identity
gene_df <- gene_df[gene_df$identity >= 35, ]

# load best bins to assign taxonomy to clusters
target_classes <- sort(c(unique(gene_df$Target)))


remov.arg_class <- c("oxazolidinone", "aminoglycoside:aminocoumarin", "rifamycin",
                     "nucleoside", "tetracenomycin_C", "fosfomycin",
                     "polyamine:peptide")

norm_gene_prev_df <- norm_gene_prev_df %>%
  filter(!(get(target_gene) %in% remov.arg_class))


#######
#Plotting

pdf(file = paste0(out.path,"/box_plot_gene_prevalence_per_library", ".pdf"),
    width = 20, height = 11.69)

ggplot(norm_gene_prev_df,
       aes(x=reorder(get(target_gene), gene.genome.prev), y=gene.genome.prev, fill=Target)) + 
  geom_boxplot() +
  xlab("Target Gene") + ylab("Normalized Target Gene Prevalence per Library") +
  theme(axis.text.x = element_text(angle = 45, size=15, hjust=1)) +
  theme(legend.text = element_text(size=15))   

dev.off()


#####

#### Perform Pairwise analyses
dat <- norm_gene_prev_df

dat <- dat %>% 
  arrange(Target) %>%
  select(Target, all_of(target_gene), gene.genome.prev)

mean_gene_target <- dat %>%
  group_by(across(all_of(target_gene)), Target) %>%
  summarise(
    n = n(),
    mean = mean(gene.genome.prev),
    sd = sd(gene.genome.prev)
  ) %>%
  ungroup()

#Do the test

stat.test <- dat %>%
  group_by(across(all_of(target_gene))) %>%
  rstatix::t_test(gene.genome.prev ~ Target, p.adjust.method = "bonferroni")

# Remove unnecessary columns and display the outputs
stat.test.select <- stat.test %>%
  select(-.y., -statistic, -df) %>%
  filter(!(p.adj.signif == "ns"))

# Create the plot
# myplot <- ggpubr::ggboxplot(
#   dat, x = "Target", y = "gene.genome.prev",
#   fill = "Target", palette = "npg", legend = "none"
# ) +
#   facet_wrap(~get(target_gene))

# # Add statistical test p-values
stat.test.select <- stat.test.select %>% rstatix::add_xy_position(x = "Target")
# myplot + ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif")

#Create individual plots
#palette = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")

graphs <- dat %>%
  group_by(across(all_of(target_gene))) %>%
  rstatix::doo(
    ~ggpubr::ggboxplot(
      data =., x = "Target", y = "gene.genome.prev",
      fill = "Target", legend = "none"
    ), 
    result = "plots"
  )
#graphs

variables <- as.character(as.data.frame(graphs)[, target_gene])
plots <- graphs$plots %>% purrr::set_names(variables)
for(variable in variables){
  
  stat.test.i <- filter(stat.test.select, get(target_gene) == variable) 
  graph.i <- plots[[variable]] + 
    labs(title = variable) +
    ggpubr::stat_pvalue_manual(stat.test.i, label = "p.adj.signif")
  
  pdf(file = paste0(out.path,"/pairwise_per_gene_boxplots/boxplot_gene_prev_per_library_", variable, ".pdf"),
      width = 12, height = 8)
  
  print(graph.i)
  
  dev.off()
}


#Save Tables to file

##Mean and SD of Target gene prev per Target metadata
write.csv(x = mean_gene_target, file = paste0(out.path, "/gene_prev_mean_per_target.csv"), row.names = F)


##t-test analysis of all pairwise comparisons
write.csv(x = stat.test, file = paste0(out.path, "/gene_per_target_pairwise_comp.csv"), row.names = F)

















