#!/bin/bash

help_message () {
	echo ""
	echo "Usage: gspreadcomp gspread [options] -o output_dir"
	echo "Options:"
	echo ""
	echo "	--checkm	STR		Path to the formatted Quality estimation dataframe"
	echo "	--gene		STR		Path to the formatted target Gene dataframe to calculate the spread"
	echo "	--gtdbtk	STR		Path to the formatted Taxonomy assignment dataframe"
	echo "	--meta		STR		Path to the formatted Sample's Metadata dataframe"
	echo "	--vf		STR		Path to the formatted Virulence Factors assignment dataframe"
	echo "	--plasmid	STR		Path to the formatted Plasmid prediction dataframe"
	echo "	--nmag		INT		Minimum number of Genomes per Library accepted [default=0]"
	echo "	--spread_taxa	STR		Taxonomic level to check gene spread [default=Phylum]"
	echo "	--target_gene_col	STR		Name of the column from the gene dataset with the Gene_ids to analyse [default=Gene_id]"
	echo "	-t			INT		number of threads"
	echo "	-o			STR		output directory"
	echo "	-h --help			print this message"
	echo ""
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
out="false"; checkm="false"; gene="false"; gtdbtk="false"; meta="false"; vf="false"; plasmid="false"; nmag=0; t=1 ;spread_taxa="Phylum";target_gene_col="Gene_id"

# load in params
OPTS=`getopt -o ht:o: --long help,checkm:,gene:,gtdbtk:,meta:,vf:,nmag:,spread_taxa:,target_gene_col:,plasmid: -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                --checkm) checkm=$2; shift 2;;
				--gene) gene=$2; shift 2;;
				--gtdbtk) gtdbtk=$2; shift 2;;
				--meta) meta=$2; shift 2;;
				--vf) vf=$2; shift 2;;
				--plasmid) plasmid=$2; shift 2;;
				--nmag) nmag=$2; shift 2;;
				--spread_taxa) spread_taxa=$2; shift 2;;
				--target_gene_col) target_gene_col=$2; shift 2;;
				-t) threads=$2; shift 2;;
                -o) out=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done


######################################################################################################o#
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# loading conda environment
echo '-------> START MODULE gspread'
#conda activate gspreadcomp_env
config_path="$(which config)"
database="${config_path/config/database}"
source $config_path
source $database

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

### Normalize and Describe sample only Mags vs Target meta ####


### Perform Complete-Contamination analysis ####

echo "Genome Quality description and Gene normalization started!"
out=`realpath $out`
initial_processing_path=$out/genome_quality_norm
mkdir $initial_processing_path

Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/simple_description_norm.r --gtdb $gtdbtk \
 --checkm $checkm \
 --gene $gene \
 --meta $meta \
 --nmag_filter $nmag \
 --target_gene_col $target_gene_col \
 --out $initial_processing_path

echo "Genome Quality description and Gene normalization finished!"


### Perform Beta dispersion analysis ####



### Perform pairwise analysis irrespective o taxa ####

echo "Gene per Target pairwise comparison started!"
out=`realpath $out`
gene_pair_path=$out/gene_pairwise_comp_results
mkdir $gene_pair_path

mkdir $gene_pair_path/pairwise_per_gene_boxplots

Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/gene_pairwise_comp.r --mags_data $initial_processing_path/genome_data_merged.csv \
 --selected_lib $initial_processing_path/selected_samples.csv \
 --gene $gene \
 --norm_gene_prev $initial_processing_path/gene_prevalence_per_library.csv \
 --spread_taxa $spread_taxa \
 --target_gene_col $target_gene_col \
 --out $gene_pair_path

echo "Gene per Target pairwise comparison finished!"


### Perform Gene prevelance spread per gOTU ####

echo "Gene prevelance spread per gOTU started!"
out=`realpath $out`
gene_spread_path=$out/gene_spread_results
mkdir $gene_spread_path

Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/gene_gotu_spread.r --mags_data $initial_processing_path/genome_data_merged.csv \
 --selected_lib $initial_processing_path/selected_samples.csv \
 --gene $gene \
 --norm_gene_prev $initial_processing_path/gene_prevalence_per_library.csv \
 --spread_taxa $spread_taxa \
 --target_gene_col $target_gene_col \
 --out $gene_spread_path

echo "Gene prevelance spread per gOTU finished!"

### Perform NMI figures ####



### Perform Pathogens - Target Gene analysis ####

patho_bac_db=$mSPREAD_DEPENDENCIES_PATH/patho_ncbi_20230222_taxa.csv

echo "Plasmid-HGT and Risk analysis started!"
out=`realpath $out`
pathogens_path=$out/pathogens_results
hgt_path=$out/hgt_events_results
vis_files=$out/network_vis_files

mkdir $pathogens_path
mkdir $hgt_path
mkdir $vis_files
mkdir $pathogens_path/pairwise_VF_per_tax_boxplots

Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/risk_analysis.r --mags_data $initial_processing_path/genome_data_merged.csv \
 --selected_lib $initial_processing_path/selected_samples.csv \
 --gene $gene \
 --norm_gene_prev $initial_processing_path/gene_prevalence_per_library.csv \
 --spread_taxa $spread_taxa \
 --target_gene_col $target_gene_col \
 --patho_db $patho_bac_db \
 --vf $vf \
 --plasmid $plasmid \
 --out $out

echo "Plasmid-HGT and Risk analysis finished!"

### Report Generation ###

echo "Report generation started!"
out=`realpath $out`

rmd_file_path=$mSPREAD_CONDA_ENVIRONMENT_PATH/bin/gspread_report_generator.Rmd

Rscript -e "rmarkdown::render('$rmd_file_path', output_file='$out/gSpread_report.html', params=list(resource_path='$out', taxa='$spread_taxa', gene_col='$target_gene_col'))"


echo "Report generation finished!"
