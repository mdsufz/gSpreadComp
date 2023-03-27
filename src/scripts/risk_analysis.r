#!/usr/bin/env Rscript

#Test path
setwd("/mnt/Chapter3_Diet/mspreadcom_r_dev/")
#setwd("//wsl.localhost/Ubuntu/home/kasmanas/mSpreadComp")

#### Load libs and inputs


library("optparse")
library("dplyr")
library("ggplot2")
library("data.table")
library("ggpubr")
library("MCDA")


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
  make_option(c("--patho_db"), type="character", default=NULL, 
              help="Dataset with the human bacterial pathogens taxonomic classification", metavar="character"),
  make_option(c("--vf"), type="character", default=NULL, 
              help="Virulence factors blastx alignment formated dataset", metavar="character"),
  make_option(c("--plasmid"), type="character", default=NULL, 
              help="Plasmid prediction dataset", metavar="character"),
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
patho_db <-data.table::fread(opt$patho_db)
vf_df <- data.table::fread(opt$vf)
plas_df <- data.table::fread(opt$plasmid)

#target_gene <- "Gene_id"
#target_gene <- "Gene_class"

#### TEST INPUT ####

#mags_data_df <- data.table::fread("03_mspread_analysis_out/genome_quality_norm/genome_data_merged.csv")
#selected_lib <- data.table::fread("03_mspread_analysis_out/genome_quality_norm/selected_samples.csv")
#gene_df <- data.table::fread("test_data/deeparg_df_format_mSpread.csv")
#gene_df <- data.table::fread("02_format_input_dfs/deeparg_df_format_mSpread.csv")
#norm_gene_prev_df <- data.table::fread("03_mspread_analysis_out/genome_quality_norm/gene_prevalence_per_library.csv")

#vf_df <- data.table::fread("02_format_input_dfs/victors_df_format_mSpread.csv")
#patho_db <- data.table::fread("/mnt/mSpreadComp_test/mSpreadComp/installation/dependencies/patho_ncbi_20230222_taxa.csv")
#plas_df <- data.table::fread("02_format_input_dfs/plasflow_df_format_mSpread.csv")

#tax_level <- "Phylum"
#out.path <- "03_mspread_analysis_out"

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

vf_df <- vf_df %>%
  filter(Genome %in% as.vector(c(unique(mags_data_df$Genome)))) %>%
  as.data.frame(.)

# load best bins to assign taxonomy to clusters
target_classes <- sort(c(unique(gene_df$Target)))

### PARTICULAR FILTERING -> REMOVE AFTER ! ###

# only consider ARGs with 35% or higher percent identity
gene_df <- gene_df[gene_df$identity >= 35, ]

# remove underscore
mags_data_df[, "Phylum"] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, "Phylum"])
mags_data_df[, "Family"] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, "Family"])
mags_data_df[, "Genus"] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, "Genus"])
mags_data_df[, "Species"] <- gsub("_([[:alpha:]]{1,2})", "", mags_data_df[, "Species"])

#2) Map NCBI Pathogens to recovered MAGs
mags_data_df$pathogen_potential <- "unclassified"

pp <- apply(X = mags_data_df, MARGIN = 1, FUN = function(x){
  pp <- "unclassified"
  
  if (x["Family"] %in% patho_db$Family) {
    pp <- "Low"
  }
  if (x["Genus"] %in% patho_db$Genus) {
    pp <- "Medium"
  }
  if (!is.na(x["Species"]) & x["Species"] %in% patho_db$Species) {
    pp <- "High"
  }
  return(pp)
  
})

mags_data_df$pathogen_potential <- pp

mags_data_df$pathogen_potential <- factor(mags_data_df$pathogen_potential,
                                          levels=c("High","Medium","Low", "unclassified"))

#Pathogens distribution per Target

path_abs_freq <- mags_data_df %>%
  group_by(Target, pathogen_potential) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))

pdf(paste0(out.path, "/pathogens_results", "barplot_patho_target_proportion_plot.pdf"),
    width = 7, height = 5)

path_abs_freq %>%
  filter(pathogen_potential != 'unclassified') %>%
  ggplot(aes(x = Target, y = prop, fill = pathogen_potential)) +
  geom_col(position = "dodge") + 
  theme(text = element_text(size = 15)) +
  labs(y = "Proportion of pathogenic Species")

dev.off()

#Count Unique VF per MAG
vf_count <- vf_df %>%
  select(Victor_VF_found, Genome, Victor_VF_class) %>%
  unique(.) %>%
  count(Genome) %>%
  as.data.frame(.) %>%
  rename(unique_vfs = n)

mags_data_df <- merge.data.frame(mags_data_df,
                                 vf_count,
                                 by = "Genome")

mean_vfs_target <- mags_data_df %>%
  group_by(across(all_of(tax_level)), Target) %>%
  summarise(
    n = n(),
    mean_unique_VFs = mean(unique_vfs),
    sd_unique_VFs = sd(unique_vfs)
  ) %>%
  ungroup()


#Do the test
#Select only tax level present in all Target
select.tax <- mean_vfs_target %>%
  group_by(across(all_of(tax_level))) %>%
  summarise(n = n_distinct(Target)) %>%
  filter(n == length(unique(mags_data_df$Target))) %>%
  ungroup() %>%
  as.data.frame(.)

select.tax <- as.character(select.tax[,1])

#### VF per Pathogen ####
mags_data_df %>%
  #filter(quality == "High") %>%
  group_by(pathogen_potential) %>%
  mutate(avg_vfs = mean(unique_vfs)) %>%
  mutate(sd_vfs = sd(unique_vfs)) %>%
  select(pathogen_potential, avg_vfs, sd_vfs) %>%
  unique(.)


pdf(paste0(out.path, "/pathogens_results", "boxplot_patho_target_uniqueVFs_plot.pdf"),
    width = 7, height = 5)

mags_data_df %>% 
  group_by(pathogen_potential) %>% 
  ggplot(aes(x = pathogen_potential,
             y = unique_vfs, fill = Target)) + 
  geom_boxplot() +
  #facet_wrap(~Target) +
  theme_bw() +
  labs(title = "", x = "Pathogen Potential", y = "Number of Unique VFs")

dev.off()
### VF : VFs per Taxa per Diet ###

pdf(paste0(out.path, "/pathogens_results", "boxplot_Phylum_target_uniqueVFs_plot.pdf"),
    width = 7, height = 5)

mags_data_df %>%
  #filter(quality == "High") %>%
  filter(pathogen_potential != "High") %>%
  ggplot(aes(x = get(tax_level), y = unique_vfs, fill = Target)) + 
  geom_boxplot() +
  #facet_wrap(~Target) +
  theme_bw() +
  labs(title = paste0("Virulence Factors (VF) per ", tax_level, "\n Removing Highly pathogenic"),
       x = tax_level, y = "Number of Unique VFs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text = element_text(size = 15))

dev.off()

pdf(paste0(out.path, "/pathogens_results", "boxplot_Phylum_filtered_target_uniqueVFs_plot.pdf"),
    width = 7, height = 5)

mags_data_df %>%
  #filter(quality == "High") %>%
  filter(pathogen_potential != "High") %>%
  filter(get(tax_level) %in% select.tax) %>%
  ggplot(aes(x = get(tax_level), y = unique_vfs, fill = Target)) + 
  geom_boxplot() +
  #facet_wrap(~Target) +
  theme_bw() +
  labs(title = paste0("Virulence Factors (VF) per ", tax_level, "\n Removing Highly pathogenic - Filtered"),
       x = tax_level, y = "Number of Unique VFs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(text = element_text(size = 15))

dev.off()


pdf(paste0(out.path, "/pathogens_results", "boxplot_Phylum_Patho_target_uniqueVFs_plot.pdf"),
    width = 7, height = 5)

mags_data_df %>%
  filter(pathogen_potential == "High") %>%
  ggplot(aes(x = get(tax_level), y = unique_vfs, fill = Target)) + 
  geom_boxplot() +
  #facet_wrap(~Target) +
  theme_bw() +
  labs(title = paste0("Virulence Factors (VF) per ", tax_level, "\n Only Highly pathogenic"),
       x = tax_level, y = "Number of Unique VFs")+
  theme(text = element_text(size = 15))

dev.off()

##
#### Perform Pairwise analyses on selected target tax level
dat <- mags_data_df %>%
  filter(pathogen_potential != "High")

dat <- dat %>% 
  arrange(Target) %>%
  select(Target, all_of(tax_level), unique_vfs) %>%
  group_by(across(all_of(tax_level)), Target) %>%
  mutate(n = n()) %>%
  ungroup(.) %>%
  filter(n >= 3) %>%
  select(-n)

mean_gene_target <- dat %>%
  group_by(across(all_of(tax_level)), Target) %>%
  summarise(
    n = n(),
    mean = mean(unique_vfs),
    sd = sd(unique_vfs)
  ) %>%
  ungroup()


#Do the test
#Select only Phylum present in all Target
select.tax <- mean_gene_target %>%
  group_by(across(all_of(tax_level))) %>%
  summarise(n = n_distinct(Target)) %>%
  filter(n == length(unique(dat$Target))) %>%
  as.data.frame(.)

select.tax <- as.character(select.tax[,1])

dat <- dat %>%
  filter(get(tax_level) %in% select.tax)

#TO DO LIKE THIS ALL PAIRWISE MUST BE POSSIBLE !
stat.test <- dat %>%
  group_by(across(all_of(tax_level))) %>%
  rstatix::t_test(unique_vfs ~ Target, p.adjust.method = "bonferroni")

# Remove unnecessary columns and display the outputs
stat.test.select <- stat.test %>%
  select(-.y., -statistic, -df) %>%
  filter(!(p.adj.signif == "ns"))

# # Add statistical test p-values
stat.test.select <- stat.test.select %>% rstatix::add_xy_position(x = "Target")
# myplot + ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif")

#Create individual plots
#palette = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")

graphs <- dat %>%
  group_by(across(all_of(tax_level))) %>%
  rstatix::doo(
    ~ggpubr::ggboxplot(
      data =., x = "Target", y = "unique_vfs",
      fill = "Target", legend = "none"
    ), 
    result = "plots"
  )
#graphs

variables <- as.character(as.data.frame(graphs)[, tax_level])
plots <- graphs$plots %>% purrr::set_names(variables)
for(variable in variables){
  
  stat.test.i <- filter(stat.test.select, get(tax_level) == variable) 
  graph.i <- plots[[variable]] + 
    labs(title = variable) +
    ggpubr::stat_pvalue_manual(stat.test.i, label = "p.adj.signif")
  
  pdf(file = paste0(out.path,"/pathogens_results/pairwise_VF_per_tax_boxplots/boxplot_tax_VF_", variable, ".pdf"),
      width = 12, height = 8)
  
  print(graph.i)
  
  dev.off()
}



#Plasmid gene and pathogens
plas_df <- plas_df %>%
  filter(Genome %in% as.vector(c(unique(mags_data_df$Genome)))) %>%
  as.data.frame(.)

all_df <- merge.data.frame(x = gene_df,
                           y = plas_df %>% select(-Genome),
                           by.x = "Gene_sequence_location",
                           by.y = "Sequence_id") %>%
  rename(Sequence_id = Gene_sequence_location)

all_df <- merge.data.frame(all_df,
                           mags_data_df %>% select(-c(Target, Library)),
                           by = "Genome")

all_df <- all_df %>%
  group_by(Genome) %>%
  mutate(unique_gene_target = n_distinct(Gene_id))

####
all_df.vfs <- merge.data.frame(vf_df,
                               y = plas_df %>% select(-Genome),
                               by.x ="Sequence_id" ,
                               by.y = "Sequence_id")

all_df.vfs <- merge.data.frame(all_df.vfs,
                               y = mags_data_df,
                               by.x ="Genome" ,
                               by.y = "Genome")


all_df.vfs <- merge.data.frame(all_df.vfs,
                               y = all_df,
                               by = c("Genome", "Sequence_id", "Library", "Target", "sequence_type", "Completeness",
                                      "Contamination", "Strain heterogeneity", "Domain", "Phylum", "Class",
                                      "Order", "Family", "Genus", "Species", "quality", "quality.score",
                                      "pathogen_potential", "unique_vfs"),
                               all = TRUE)




all_df.vfs <- all_df.vfs %>%
  group_by(Genome) %>%
  mutate(unique_gene_target = n_distinct(Gene_id, na.rm = T))



vfs_per_seq_type_count <- all_df.vfs %>%
  group_by(Genome, sequence_type) %>%
  summarise(n=n_distinct(Victor_VF_found, na.rm = T)) %>%
  pivot_wider(names_from = sequence_type, values_from = n) %>%
  replace(is.na(.), 0) %>%
  rename(unique_vf_in_chromosome = chromosome,
         unique_vf_in_plasmid = plasmid,
         unique_vf_in_unclassified = unclassified)


gene_per_seq_type_count <- all_df %>%
  group_by(Genome, sequence_type) %>%
  summarise(n=n_distinct(Gene_id)) %>%
  pivot_wider(names_from = sequence_type, values_from = n) %>%
  replace(is.na(.), 0) %>%
  rename(unique_gene_in_chromosome = chromosome,
         unique_gene_in_plasmid = plasmid,
         unique_gene_in_unclassified = unclassified)


all_df.vfs <- merge.data.frame(all_df.vfs,
                               gene_per_seq_type_count,
                               by = "Genome")

all_df.vfs <- merge.data.frame(all_df.vfs,
                               vfs_per_seq_type_count,
                               by = "Genome")



all_df <- merge.data.frame(all_df,
                           vfs_per_seq_type_count,
                           by = "Genome")

all_df <- merge.data.frame(all_df,
                           gene_per_seq_type_count,
                           by = "Genome")

##################
mean_target_gene_per_patho <- all_df %>%
  select(Genome, Library, Target, unique_gene_target, unique_vfs, sequence_type, pathogen_potential) %>%
  unique(.) %>%
  group_by(Target, sequence_type, pathogen_potential) %>%
  summarise(mean_unique_gene_target = mean(unique_gene_target),
            sd_unique_gene_target = sd(unique_gene_target),
            n = n())


# mean_target_gene_per_patho %>%
#   filter(n>1) %>%
#   ggplot(aes(x = pathogen_potential, y = mean_unique_gene_target, fill = Target)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_errorbar(aes(ymin = mean_unique_gene_target - sd_unique_gene_target,
#                     ymax = mean_unique_gene_target + sd_unique_gene_target),
#                 width = 0.4, position = position_dodge(0.9)) +
#   facet_wrap(~sequence_type)


all_df$pathogen_potential <- as.factor(all_df$pathogen_potential)

pdf(paste0(out.path, "/pathogens_results", "boxplot_Patho_target_uniqueGene_plot.pdf"),
    width = 7, height = 5)

all_df %>%
  select(Genome, pathogen_potential, unique_gene_target, Target) %>%
  unique(.) %>%
  ggplot(aes(x = pathogen_potential, y = unique_gene_target, fill = Target)) + 
  geom_boxplot() +
  #facet_wrap(~Target) +
  theme_bw() +
  labs(title = "Target Gene in Pathogens", x = "Pathogen Potential", y = "Total Unique Target Gene")

dev.off()

stat.test.gene.patho <- all_df %>%
  select(pathogen_potential, unique_gene_target, Genome) %>%
  unique(.) %>%
  as.data.frame(.) %>%
  #group_by(pathogen_potential) %>%
  rstatix::t_test(unique_gene_target ~ pathogen_potential, p.adjust.method = "bonferroni")

stat.test.gene.patho <- stat.test.gene.patho %>%
  filter(p.adj.signif != 'ns') %>%
  rstatix::add_xy_position(x = "pathogen_potential")


pdf(paste0(out.path, "/pathogens_results", "boxplot_Patho_uniqueGene_plot.pdf"),
    width = 7, height = 5)

ggboxplot(all_df %>%
            select(pathogen_potential, unique_gene_target, Genome) %>%
            unique(.) %>%
            as.data.frame(.), x = "pathogen_potential", y = "unique_gene_target",
          fill = "pathogen_potential",
          xlab = "Pathogen Potential", ylab = "Total Unique Target Gene") +
  stat_pvalue_manual(stat.test.gene.patho) +
  theme(legend.position="none")

dev.off()


#######################
#### Calculate Risk####
#######################

performanceTable <- all_df %>%
  select(Genome, pathogen_potential, unique_vfs,
         unique_gene_target, unique_gene_in_chromosome, unique_gene_in_plasmid, unique_gene_in_unclassified,
         unique_vf_in_chromosome, unique_vf_in_plasmid, unique_vf_in_unclassified) %>%
  unique(.)

row.names(performanceTable) <- performanceTable$Genome
performanceTable <- performanceTable %>% select(-Genome)

mean_per_patho <- performanceTable %>%
  group_by(pathogen_potential) %>%
  summarise(across(everything(), mean)) %>%
  as.data.frame(.)

row.names(mean_per_patho) <- mean_per_patho$pathogen_potential
mean_per_patho <- mean_per_patho %>% 
  select(-pathogen_potential) 

mean_per_patho_norm <- as.data.frame(lapply(mean_per_patho, function(x) x/sum(x)))

#weights were defined based on the highest proportion of every attribute in the pathogens species
weights_patho <- as.numeric(mean_per_patho_norm[1,]) #The first row is the pathogens High
weights_patho <- weights_patho / sum(weights_patho)

criteriaMinMax <- c("max", "max", "max", "max", "max", "max", "max", "max")

#positiveIdealSolutions <- c(0.179573776, 0.171636015, 0.159499658, 0.087302767)
#negativeIdealSolutions <- c(0.212610118, 0.124958799, 0.131352659, 0.085797547)

performanceTable <- performanceTable %>%
  select(colnames(mean_per_patho))

names(weights_patho) <- colnames(performanceTable)
names(criteriaMinMax) <- colnames(performanceTable)
#names(positiveIdealSolutions) <- colnames(performanceTable)
#names(negativeIdealSolutions) <- colnames(performanceTable)

risk_criteria <- TOPSIS(performanceTable = performanceTable,
                        criteriaWeights = weights_patho,
                        criteriaMinMax = criteriaMinMax)

risk_criteria <- as.data.frame(risk_criteria)

risk_criteria$Genome <- row.names(risk_criteria)


all_df <- merge.data.frame(all_df,
                           risk_criteria,
                           by = "Genome")


all_df.vfs <- merge.data.frame(all_df.vfs,
                               risk_criteria,
                               by = "Genome")



mags_data_df <- all_df %>%
  select(c("Genome", "Completeness", "Contamination", "Strain heterogeneity", "quality", "quality.score",
           "Library", "Target",
           "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species",
           "pathogen_potential","risk_criteria", "unique_vfs", "unique_gene_target",
           "unique_gene_in_chromosome", "unique_gene_in_plasmid",
           "unique_gene_in_unclassified",
           "unique_vf_in_chromosome", "unique_vf_in_plasmid", 
           "unique_vf_in_unclassified")) %>%
  unique(.)

#Plots for the calculated risk measure

pdf(paste0(out.path, "/pathogens_results", "boxplot_Patho_uniqueGene_per_seqtype_plot.pdf"),
    width = 7, height = 5)

all_df %>%
  select(pathogen_potential, unique_gene_target, Genome, Target, sequence_type) %>%
  unique(.) %>%
  ggplot(aes(x = pathogen_potential, y = unique_gene_target, fill = Target)) + 
  geom_boxplot() +
  facet_wrap(~sequence_type) +
  theme_bw() +
  labs(title = "Target Gene based on Sequence type", x = "Pathogen Potential", y = "Unique Gene Targets")

dev.off()


all_df %>%
  select(pathogen_potential, risk_criteria, Genome, Target, Phylum) %>%
  unique(.) %>%
  filter(pathogen_potential != "High") %>%
  filter(Phylum %in% select.tax) %>%
  ggplot(aes(x = Phylum, y = risk_criteria, fill = Target)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Risk Measure per Phylum", x = "Phylum", y = "Risk Measure")


all_df %>%
  select(pathogen_potential, risk_criteria, Genome, Target, Phylum) %>%
  unique(.) %>%
  filter(pathogen_potential != "High") %>%
  #filter(Phylum %in% select.tax) %>%
  ggplot(aes(x = Target, y = risk_criteria, fill = Target)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Risk Measure per Target", x = "Target", y = "Risk Measure")


##### Find the HGT in plasmids ####
cycle_sample_finder <- function(mags_info_table,
                                sample_name_col = "Library",
                                taxa_level_col = "Family",
                                target_gene_col = "Gene_id"){
  
  #5.x) Find cycles in the form of: Sample -> Taxa1 -> ARG -> Taxa2 -> Sample
  cycles.df <- data.frame(matrix(nrow = 0, ncol = 4))
  
  names(cycles.df) <-  c(sample_name_col,
                         paste0(taxa_level_col,"1"),
                         target_gene_col,
                         paste0(taxa_level_col,"2"))
  
  
  n_iter <- length(unique(mags_info_table[, sample_name_col])) # Number of iterations of the loop
  
  counter <- 0
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       #width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for (s in unique(mags_info_table[, sample_name_col])) {
    
    i <- which(unique(mags_info_table[, sample_name_col]) == s)
    #select the taxa connected to sample s
    tax <- mags_info_table %>% 
      filter(.[[sample_name_col]] == s) %>%
      select(all_of(taxa_level_col)) %>%
      unique(.)
    
    tax <- as.vector(tax[,1])
    
    #Filter if tax > 2
    if (length(tax) > 1) {
      
      #Loop through tax and find connections to ARG
      for (t1 in 1:(length(tax)-1)) {
        
        #select gene connected to t1 in sample s
        gene.t1 <- mags_info_table %>%
          filter(.[[taxa_level_col]] == tax[t1]) %>%
          filter(.[[sample_name_col]] == s) %>%
          select(all_of(target_gene_col)) %>%
          unique(.)
        
        gene.t1 <- as.vector(gene.t1[,target_gene_col])
        
        for (t2 in (t1+1):length(tax)) {
          
          #select gene connected to t2 in sample s
          gene.t2 <- mags_info_table %>%
            filter(.[[taxa_level_col]] == tax[t2]) %>%
            filter(.[[sample_name_col]] == s) %>%
            select(all_of(target_gene_col)) %>%
            unique(.)
          
          gene.t2 <- as.vector(gene.t2[,target_gene_col])
          
          #Check if t1 and t2 are connected to the same arg in sample s
          gene.both <- gene.t1[gene.t1 %in% gene.t2]
          counter <- counter + 1
          
          #Create cycle if gene.both exists
          if (length(gene.both)) {
            
            for (gene in gene.both) {
              
              row <- c(s, tax[t1],
                       gene, tax[t2])
              cycles.df <- rbind.data.frame(cycles.df, row)
              counter <- counter + 1
              
            }
            
          }
          
        }
        
      }
      
    }
    
    # Sets the progress bar to the current state
    setTxtProgressBar(pb, i)
    
  }
  names(cycles.df) <- c(sample_name_col, paste0(taxa_level_col,"1"),
                        target_gene_col, paste0(taxa_level_col,"2"))
  
  close(pb) # Close the connection
  total_count <- counter
  return(list(cycles.df, total_count))
  
}

#### Find VFs involved in plasmid HGT ####
print("Finding Virulence factors involved in plasmid-mediated HGT:")

cycles_found_list_vf <- list()
all_cycles_found_df_vf <- tibble()
tax_level = "Family"
for (t in target_classes) {
  
  print(paste0("VFs plasmid-HGT calculation for ", t, ":"))
  
  cycles_found <- cycle_sample_finder(mags_info_table = as.data.frame(all_df.vfs %>%
                                                                        filter(Target == t) %>%
                                                                        filter(sequence_type == "plasmid") %>%
                                                                        drop_na(Victor_VF_found)),
                                      sample_name_col = "Library",
                                      taxa_level_col = tax_level,
                                      target_gene_col = "Victor_VF_found")
  
  total_plas_df <- cycles_found[[1]]
  
  total_paths <- cycles_found[[2]]
  
  perc_cycles <- nrow(cycles_found[[1]]) / cycles_found[[2]]
  
  cycles_found_list_vf[[t]][["percentage_of_plasmids_carring_gene"]] <- perc_cycles
  cycles_found_list_vf[[t]][["plasmid_HGT_df"]] <- total_plas_df
  
  total_plas_df$Target <- t
  
  if(nrow(total_plas_df) > 0){
    
    all_cycles_found_df_vf <- rbind.data.frame(all_cycles_found_df_vf, total_plas_df)
    
  }
}


#### Find Target Gene involved in plasmid HGT ####
cycles_found_list <- list()
all_cycles_found_df <- tibble()
tax_level = "Family"
for (t in target_classes) {
  
  print(paste0("Target Gene plasmid-HGT calculation for ", t, ":"))
  
  cycles_found <- cycle_sample_finder(mags_info_table = as.data.frame(all_df %>%
                                                                        filter(Target == t) %>%
                                                                        filter(sequence_type == "plasmid") %>%
                                                                        drop_na(Gene_id)),
                                      sample_name_col = "Library",
                                      taxa_level_col = tax_level,
                                      target_gene_col = "Gene_id")
  
  total_plas_df <- cycles_found[[1]]
  
  total_paths <- cycles_found[[2]]
  print(total_paths)
  
  perc_cycles <- nrow(cycles_found[[1]]) / cycles_found[[2]]
  
  cycles_found_list[[t]][["percentage_of_plasmids_carring_gene"]] <- perc_cycles
  cycles_found_list[[t]][["plasmid_HGT_df"]] <- total_plas_df
  
  total_plas_df$Target <- t
  
  if(nrow(total_plas_df) > 0){
    
    all_cycles_found_df <- rbind.data.frame(all_cycles_found_df, total_plas_df)
    
  }
  
}

#Creating Network files per Target

for (t in target_classes) {
  print(paste0("Creating network Visualization files for ", t))
  dir.create(paste0(out.path, "/network_vis_files/", t))
  all_df.vfs.tmp <- all_df.vfs %>%
    filter(Target == t) %>%
    unique(.) %>%
    as.data.frame(.)
  
  #File generated in pathogen script
  nodes.genomes <- all_df.vfs.tmp %>%
    select(Genome, risk_criteria, Target, Library, 
           Domain, Phylum, Class,
           Order, Family, Genus, Species,
           pathogen_potential) %>%
    unique(.) %>%
    as.data.frame(.) %>%
    rename(ID = Genome)
  
  
  nodes.vfs <- all_df.vfs.tmp %>%
    select(Victor_VF_found) %>%
    unique(.) %>%
    mutate(risk_criteria = min(nodes.genomes$risk_criteria), Target = "VF", Library= "VF", 
           Domain= "VF", Phylum= "VF", Class= "VF",
           Order= "VF", Family= "VF", Genus= "VF", Species= "VF",
           pathogen_potential= "VF") %>%
    as.data.frame(.) %>%
    rename(ID = Victor_VF_found)
  
  nodes.gene <- all_df.vfs.tmp %>%
    select(Gene_id) %>%
    unique(.) %>%
    mutate(risk_criteria = min(nodes.genomes$risk_criteria), Target = "Gene", Library= "Gene", 
           Domain= "Gene", Phylum= "Gene", Class= "Gene",
           Order= "Gene", Family= "Gene", Genus= "Gene", Species= "Gene",
           pathogen_potential= "Gene") %>%
    as.data.frame(.) %>%
    rename(ID = Gene_id)
  
  nodes.gene.genomes <- rbind.data.frame(nodes.genomes, nodes.gene)
  #Create edges
  
  edge.genome.vfs <- table(all_df.vfs.tmp$Genome,
                           all_df.vfs.tmp$Victor_VF_found) %>%
    as.data.frame(.) %>%
    rename(Source = Var1, Target = Var2) %>%
    filter(Freq > 0)
  
  
  edge.genome.gene <- table(all_df.vfs.tmp$Genome,
                            all_df.vfs.tmp$Gene_id) %>%
    as.data.frame(.) %>%
    rename(Source = Var1, Target = Var2) %>%
    filter(Freq > 0)
  
  
  edge.gene.vf <- table(all_df.vfs.tmp$Gene_id,
                        all_df.vfs.tmp$Victor_VF_found) %>%
    as.data.frame(.) %>%
    rename(Source = Var1, Target = Var2) %>%
    filter(Freq > 0)
  
  
  edge.total <- rbind.data.frame(edge.genome.vfs,
                                 edge.genome.gene,
                                 edge.gene.vf)
  
  
  edge.total.fil <- edge.total %>%
    filter(Freq > 1)
  
  
  ##Save to  csv
  
  write_csv(nodes.genomes, file = paste0(out.path, "/network_vis_files/", t ,"_nodes_genomes.csv"), row.names = F)
  
  write_csv(nodes.gene, file = paste0(out.path, "/network_vis_files/",t ,"_nodes_genes.csv"), row.names = F)
  
  write_csv(nodes.vfs, file = paste0(out.path,"/network_vis_files/",t,"_nodes_vfs.csv"), row.names = F)
  
  write_csv(nodes.gene.genomes, file = paste0(out.path,"/network_vis_files/",t,"_nodes_gene_genome.csv"), row.names = F)
  
  write_csv(edge.total, file = paste0(out.path,"/network_vis_files/",t ,"_edge_all.csv"), row.names = F)
  
  write_csv(edge.total.fil, file = paste0(out.path,"/network_vis_files/",t ,"_edge_all_fil.csv"), row.names = F)
  
  write_csv(edge.genome.vfs, file = paste0(out.path,"/network_vis_files/",t ,"_edge_genome_vfs.csv"), row.names = F)
  
  write_csv(edge.genome.gene, file = paste0(out.path,"/network_vis_files/",t ,"_edge_genome_gene.csv"), row.names = F)
  
  write_csv(edge.gene.vf, file = paste0(out.path,"/network_vis_files/",t ,"_edge_gene_vf.csv"), row.names = F)
  
}


#Save tables to files


## path_abs_freq
write.csv(x = path_abs_freq,
          file = paste0(out.path, "/pathogens_results/", "/target_pathogens_proportion.csv"), row.names = F)

## mean_vfs_target
write.csv(x = mean_vfs_target, 
          file = paste0(out.path, "/pathogens_results/", "/target_tax_VFs.csv"), row.names = F)

## select.tax 
write.csv(x = select.tax,
          file = paste0(out.path, "/common_tax_target.csv"), row.names = F)

## stat.test.select 
write.csv(x = stat.test.select,
          file = paste0(out.path, "/pathogens_results/", "/stat_diff_tax_vf_target.csv"), row.names = F)

## stat.test.gene.patho
write.csv(x = stat.test.gene.patho,
          file = paste0(out.path, "/pathogens_results/", "/stat_diff_patho_potential.csv"), row.names = F)

## all_cycles_found_df_vf
write.csv(x = all_cycles_found_df_vf,
          file = paste0(out.path, "/hgt_events_results/", "/all_hgt_found_vf.csv"), row.names = F)

## all_cycles_found_df
write.csv(x = all_cycles_found_df,
          file = paste0(out.path, "/hgt_events_results/", "/all_hgt_found_gene.csv"), row.names = F)

## all_df.vfs
write.csv(x = all_df.vfs,
          file = paste0(out.path, "/magos_complete_annotation.csv"), row.names = F)

## mags_data_df 
write.csv(x = mags_data_df,
          file = paste0(out.path, "/mags_summary_results.csv"), row.names = F)





