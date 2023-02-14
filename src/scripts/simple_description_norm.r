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
library("patchwork")
library("data.table")
library("viridis")
library("gridExtra")
library("ggforce")


#target_gene <- "Gene_id"
#target_gene <- "Gene_class"

#Young samples to remove -> REMOVE ON THE FINAL VERSION

young_samples <- c("EA_SRX3062589", "EA_SRX3062590", "EA_SRX3062593", "EA_SRX3062617", "EA_SRX3062647", 
                   "EA_SRX3062648", "EA_SRX3062649", "EA_SRX3062650", "EA_SRX3062651", "EA_SRX3062652", 
                   "EA_SRX3062653", "EA_SRX3062654", "EA_SRX3062655", "EA_SRX3062656", "EA_SRX3062657", 
                   "EA_SRX3062658", "EA_SRX3062660", "EA_SRX3062661", "EA_SRX3062662", "EA_SRX3062663", 
                   "EA_SRX3062664")

####

option_list = list(
  make_option(c("--gtdb"), type="character", default=NULL, 
              help="GTDB-tk Taxa classification file path", metavar="character"),
  make_option(c("--checkm"), type="character", default=NULL, 
              help="CheckM quality result file path", metavar="character"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="gene merged result for all genomes file path", metavar="character"),
  make_option(c("--meta"), type="character", default=NULL, 
              help="Library metadata file path", metavar="character"),
  make_option(c("--nmag_filter"), type="integer", default=0, 
              help="Minimum number of mags per Library. [default = 0]", metavar="character"),
  make_option(c("--target_gene_col"), type="character", default=NULL, 
              help="Name of the Target Gene column [default: Gene_id]", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output directory path to save formated tables", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.path <- opt$out
gtdb_df <- data.table::fread(opt$gtdb)
checkm_df <- data.table::fread(opt$checkm)
gene_df <- data.table::fread(opt$gene)
meta_df <- data.table::fread(opt$meta)
n_mag <- opt$nmag_filter
target_gene <- opt$target_gene_col


#1) Load metadata, mags taxa, and mags quality

############ TEST ##############

# checkm_df <- data.table::fread(file = "test_data/checkm_df_format_mSpread.csv")
# 
# #####
# gtdb_df <- data.table::fread(file = "test_data/gtdb_df_format_mSpread.csv")
# 
# #####
# gene_df <- data.table::fread(file = "test_data/deeparg_df_format_mSpread.csv")
# 
# #####
# meta_df <- data.table::fread(file = "test_data/meta_df_format_mSpread.csv")

meta_df <- meta_df %>%
  filter(!(Library %in% young_samples))

####2) Merger
mags_data <- merge.data.frame(x = checkm_df,
                              y = meta_df,
                              by = "Genome")

mags_data <- merge.data.frame(x = mags_data,
                              y = gtdb_df,
                              by = "Genome")


target_classes <- sort(as.vector(unique(mags_data$Target)))

#3) create plot of mag quality
mags_data <- as.data.frame(mags_data)
mags_data$quality <- "Low"
mags_data$quality[(mags_data$Completeness >= 50) & (mags_data$Contamination < 10)] <- "Medium"
mags_data$quality[(mags_data$Completeness > 90) & (mags_data$Contamination < 5)] <- "High"
mags_data$quality.score <-mags_data$Completeness - 5*mags_data$Contamination


# summarize HQ and MQ MAGs

mag.quality <- as.data.frame(mags_data %>%
                                    group_by(quality, Target) %>%
                                    summarize(count = n_distinct(Genome)) %>% 
                                    mutate(percent = round(prop.table(count) * 100, 1)))

mag.quality <- as.data.frame(mags_data %>%
                               group_by(quality) %>%
                               summarize(count = n_distinct(Genome)) %>% 
                               mutate(percent = round(prop.table(count) * 100, 1)))

mag.scatter <- mags_data %>%  ggplot(., aes(x = Contamination, y = Completeness, color = Target)) +
  geom_point(size = 2, show.legend = TRUE) +
  facet_wrap(~Target) +
  #scale_color_viridis(option = "viridis") +
  scale_y_continuous(limits = c(50, 100)) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Contamination (%)") + ylab("Completeness (%)") +
  theme(legend.position="none")
#annotation_custom(tableGrob(mag.quality, rows = NULL, theme = ttheme_minimal()), 
#                  xmin=5, xmax=10, ymin=50, ymax=70)

mag.scatter.quality <- mags_data %>%  ggplot(., aes(x = Contamination, y = Completeness, color = quality.score)) +
  geom_point(size = 2, show.legend = TRUE) +
  facet_wrap(~Target) +
  scale_color_viridis(option = "viridis") +
  scale_y_continuous(limits = c(50, 100)) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Contamination (%)") + ylab("Completeness (%)") 


quality.t.p <- tableGrob(mag.quality, rows = NULL, theme = ttheme_minimal())
grid.arrange(quality.t.p)

mag.completeness.hist <- ggplot(mags_data, aes(x = Completeness, fill = cut(Completeness, 100))) +
  facet_wrap(~Target) +
  geom_histogram(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis(option = "viridis", discrete = T) +
  theme_void() +
  scale_x_continuous(limits = c(50, 100)) +
  coord_flip() +
  theme(strip.text = element_blank(), strip.background = element_blank())

#mag.completeness.hist

mag.contamination.hist <- ggplot(mags_data, aes(x = Contamination, fill = cut(Contamination, 100))) + 
  geom_histogram(bins = 200, show.legend = FALSE) +
  facet_wrap(~Target) +
  scale_fill_viridis(option = "viridis", discrete = T, direction = -1) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_void() +
  theme(strip.text = element_blank(), strip.background = element_blank())

#mag.contamination.hist

#2) Filter out samples based on number of MAGs
#Calculate total number of bins recovered per sample
mags.per.sample <- meta_df %>%
  select(Library, Genome) %>%
  unique(.) %>%
  group_by(Library) %>%
  mutate(t_mags = n()) %>%
  select(-Genome) %>%
  unique(.) %>%
  merge.data.frame(.,
                   meta_df,
                   by = "Library") %>%
  select(Library, Target, t_mags) %>%
  unique(.)


#Total number of samples per diet type
table(mags.per.sample$Target)

#Calculate mags per sample mean per class
mags.per.sample %>%
  group_by(Target) %>%
  summarise_at(vars(t_mags), list(name = mean))

#Keep samples only with 10 mags or more and for the ancient 5 or more
mags.per.sample.fil <- mags.per.sample %>%
  filter(t_mags > n_mag)

mags.per.sample.fil %>%
  group_by(Target) %>%
  summarise_at(vars(t_mags), list(name = mean))

#Total number of samples per diet type after filtering
table(mags.per.sample.fil$Target)

#Selected samples
samples.fil <- mags.per.sample.fil$Library


# Normalize prevalence
mags.per.sample.gene.prev <- gene_df %>%
  merge.data.frame(.,
                   meta_df,
                   by = "Genome") %>%
  select(Library, all_of(target_gene), Genome) %>%
  unique(.) %>%
  group_by(Library, across(all_of(target_gene))) %>%
  mutate(present.gene = n()) %>%
  select(-Genome) %>%
  unique(.) %>%
  merge.data.frame(.,
                   mags.per.sample,
                   by = "Library") %>%
  mutate(gene.genome.prev = present.gene/t_mags)

##########
# count number of unique args per mag
avg_unique_genes <- tibble()

for (t in target_classes) {
  
  mag.gene.cnt <- as.data.frame(
    gene_df %>%
      merge.data.frame(.,
                       mags_data[,c("Genome", "Target")]) %>%
      filter(Target == t) %>%
      group_by(Genome) %>% 
      summarize(unique.gene = n_distinct(across(all_of(target_gene)), na.rm = T))
  )
  
  # create histogram with number of args per mag
  pdf(file = paste0(out.path, "/target_gene_per_genome_hist_", t, "_plot.pdf"),
      height = 8.27, width = 11.69)
  
  print(
  ggplot(mag.gene.cnt, aes(x = unique.gene)) + 
    geom_histogram(binwidth = 1, color = "darkcyan", fill = "darkturquoise") + 
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      text = element_text(size = 20)
    ) +
    labs(x = "Number of Unique Target Genes per MAG", y = "Number of MAGs") +
    ggtitle(t)
  )
  
  dev.off()
  
  # average number of Target Gene per gOTU MAGs
  avg_t <- mean(mag.gene.cnt$unique.gene)
  sd_t <- sd(mag.gene.cnt$unique.gene)
  
  avg_unique_genes <- rbind.data.frame(avg_unique_genes, c(t, avg_t, sd_t))

}
names(avg_unique_genes) <- c("Target", "Avg_unique_target_gene", "Standard_deviation")


#Save Figures and Objects to files

#files

##Merged Genome info dataset
write.csv(x = mags_data, file = paste0(out.path, "/genome_data_merged.csv"), row.names = F)

##Selected samples after filtering
write.csv(x = samples.fil, file = paste0(out.path, "/selected_samples.csv"), row.names = F)

##Prev gene table
write.csv(x = mags.per.sample.gene.prev, file = paste0(out.path, "/gene_prevalence_per_library.csv"), row.names = F)

##Avg unique gene per target table
write.csv(x = avg_unique_genes, file = paste0(out.path, "/avg_unique_gene_per_target.csv"), row.names = F)

#Figures
###
for (i in 1:length(target_classes)) {
  
  pdf(paste0(out.path, "/quality_scatter_", target_classes[i], "_plot.pdf"),
      width = 7, height = 5)
  
  print(
    mag.contamination.hist +
      facet_wrap_paginate(~Target, ncol = 1, nrow = 1, page = i) +
      plot_spacer() +
      mag.scatter +
      facet_wrap_paginate(~Target, ncol = 1, nrow = 1, page = i) +
      mag.completeness.hist +
      facet_wrap_paginate(~Target, ncol = 1, nrow = 1, page = i) +
      plot_layout(
        ncol = 2, 
        nrow = 2, 
        widths = c(4, 1),
        heights = c(1, 4)
      )
  )
  
  dev.off()
}

###
#Filter out samples that recovered less than x MAGs per
pdf(paste0(out.path, "/before_nGenomes_violin", "_plot.pdf"),
    width = 7, height = 5)

mags.per.sample %>%
  ggplot(aes(x=Target, y=t_mags, fill=Target)) +
  geom_violin() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  #theme_ipsum() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("MAGs per sample") +
  xlab("")+
  labs(y = "Number of MAGs in a Library")

dev.off()

###
pdf(paste0(out.path, "/after_nGenomes_violin", "_plot.pdf"),
    width = 7, height = 5)

mags.per.sample.fil %>%
  ggplot(aes(x=Target, y=t_mags, fill=Target)) +
  geom_violin() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  #theme_ipsum() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("MAGs per sample filtered") +
  xlab("")+
  labs(y = "Number of MAGs in a Library")

dev.off()







		   
