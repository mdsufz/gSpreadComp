#!/usr/bin/env Rscript

#Select Samples based on the number of MAGs recovered
#Describe MAGs and Save figures

#### Load libs and inputs
suppressPackageStartupMessages({
library("optparse")
library("dplyr")
library("tidyr")
library("ggplot2")
library("data.table")
library("viridis")
library("pheatmap")
library("forcats")
library("rstatix")
library("tidyverse")
library("purrr")
library("progress")})


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
              help="output directory path to save formatted files", metavar="character")
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

# Only consider ARGs with 35% or higher percent identity
gene_df <- gene_df[gene_df$identity >= 35, ]

# Load best bins to assign taxonomy to clusters
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
#Create a function to remove special char
sanitize_string <- function(string) {
  # Use gsub to replace non-alphanumeric characters (excluding spaces) with an empty string
  sanitized <- gsub("[^[:alnum:][:space:]]", "", string)
  return(sanitized)
}

########## Speedup function ####
perform_t_tests <- function(df, target_gene) {
  
  compute_pairwise <- function(gene_class, sub_df, pb) {
    targets <- unique(sub_df$Target)
    
    if (length(targets) < 2) {
      pb$tick(length(targets))
      map_df(targets, ~ setNames(data.frame(
        gene_class = gene_class,
        .y. = NA,
        group1 = NA,
        group2 = NA,
        n1 = NA,
        n2 = NA,
        statistic = NA,
        df = NA,
        p = NA,
        notes = paste("Not enough groups for comparison within this", target_gene)
      ), c(target_gene, ".y.", "group1", "group2", "n1", "n2", "statistic", "df", "p", "notes")))
    } else {
      pb$tick(choose(length(targets), 2))
      combn(targets, 2, simplify = FALSE) %>%
        map_df(~ {
          data_pair <- sub_df %>% filter(Target %in% .x)
          data_pair$Target <- droplevels(data_pair$Target)
          
          if (all(table(data_pair$Target) >= 2)) {
            data_pair %>%
              group_by(!!sym(target_gene)) %>%
              rstatix::t_test(
                gene.genome.prev ~ Target, 
                p.adjust.method = "bonferroni"
              ) %>%
              mutate(notes = NA)
          } else {
            setNames(data.frame(
              gene_class = gene_class,
              .y. = NA,
              group1 = .x[1],
              group2 = .x[2],
              n1 = NA,
              n2 = NA,
              statistic = NA,
              df = NA,
              p = NA,
              notes = "Not enough observations in each group for comparison"
            ), c(target_gene, ".y.", "group1", "group2", "n1", "n2", "statistic", "df", "p", "notes"))
          }
        })
    }
  }
  
  # Calculate total steps for the progress bar
  total_steps <- df %>%
    group_by(!!sym(target_gene)) %>%
    summarise(n_combinations = sum(ifelse(n_distinct(Target) < 2, n_distinct(Target), choose(n_distinct(Target), 2)))) %>%
    pull(n_combinations) %>%
    sum()
  
  # Create a progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent Elapsed: :elapsed",
    total = total_steps,
    width = 60
  )
  
  results_df <- df %>%
    group_by(!!sym(target_gene)) %>%
    group_map(~ compute_pairwise(.x[[target_gene]][1], .x, pb), .keep = TRUE) %>%
    bind_rows()
  
  results_df$p_adjusted <- p.adjust(results_df$p, method = "bonferroni")
  
  results_df <- results_df %>%
    mutate(
      p.adj.signif = case_when(
        is.na(p_adjusted) ~ NA_character_,
        p_adjusted < 0.001 ~ '***',
        p_adjusted >= 0.001 & p_adjusted < 0.01 ~ '**',
        p_adjusted >= 0.01 & p_adjusted < 0.05 ~ '*',
        p_adjusted >= 0.05 & p_adjusted < 0.1 ~ '.',
        p_adjusted >= 0.1 ~ 'ns'
      )
    )
  
  return(results_df)
}

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

# Reorder the levels of the Target factor based on alphabetical order
dat$Target <- factor(dat$Target, levels = sort(unique(dat$Target)))

results_df <- perform_t_tests(dat, target_gene)

#### Create Graphs ####

graphs <- dat %>%
  group_by(across(all_of(target_gene))) %>%
  rstatix::doo(
    ~ggpubr::ggboxplot(
      data =., x = "Target", y = "gene.genome.prev",
      fill = "Target", legend = "none"
    ), 
    result = "plots"
  )

######
# Identify all unique target groups
all_targets <- sort(unique(dat$Target))

# Generate a consistent color palette for all target groups
color_palette <- scales::hue_pal()(length(all_targets))
names(color_palette) <- all_targets

# Annotate the plots using the provided significance levels in results_df
graphs$plots <- map2(graphs$plots, graphs[[target_gene]], function(p, current_gene_class_value) {
  
  # 1. Extract Relevant Annotations
  annotations_for_this_gene_class <- results_df %>% 
    filter(!!sym(target_gene) == current_gene_class_value) %>% 
    select(group1, group2, p.adj.signif)
  
  # Initialize a base y-position
  base_y_position <- max(p$data$gene.genome.prev, na.rm=TRUE) + 0.05
  
  # 2. Iterate Over Each Annotation
  for (i in 1:nrow(annotations_for_this_gene_class)) {
    ann <- annotations_for_this_gene_class[i, ]
    
    # Skip "ns" annotations or NA values
    if (is.na(ann$p.adj.signif) || ann$p.adj.signif == "ns") {
      next
    }
    
    # 3. Add Annotations to the Plot
    p <- p + geom_signif(
      comparisons = list(c(ann$group1, ann$group2)), 
      test = NULL, 
      annotations = ann$p.adj.signif, 
      y_position = base_y_position  # Use the current base y-position
    )
    
    # Increment the base y-position for the next annotation
    base_y_position <- base_y_position + 0.05
  }
  
  # Use a consistent color palette and add the title
  p <- p + 
    scale_fill_manual(values = color_palette) + 
    ggtitle(current_gene_class_value)
  
  pdf_file_path <- paste0(out.path, "/pairwise_per_gene_boxplots/boxplot_gene_prev_per_library_", 
                          sanitize_string(current_gene_class_value), ".pdf")
  pdf(file = pdf_file_path, width = 12, height = 8)
  print(p)
  dev.off()
  
  return(p)
})

#Save Tables to file

##Mean and SD of Target gene prev per Target metadata
write.csv(x = mean_gene_target, file = paste0(out.path, "/gene_prev_mean_per_target.csv"), row.names = F)

##t-test analysis of all pairwise comparisons
write.csv(x = results_df, file = paste0(out.path, "/gene_per_target_pairwise_comp.csv"), row.names = F)
