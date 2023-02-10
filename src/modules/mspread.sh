#!/bin/bash

help_message () {
	echo ""
	echo "Usage: mspreadcomp mspread [options] -o output_dir"
	echo "Options:"
	echo ""
	echo "	--checkm	STR		Path to the formatted Quality estimation dataframe"
	echo "	--gene		STR		Path to the formatted target Gene dataframe to calculate the spread"
	echo "	--gtdbtk	STR		Path to the formatted Taxonomy assignment dataframe"
	echo "	--meta		STR		Path to the formatted Sample's Metadata dataframe"
	echo "	--vf		STR		Path to the formatted Virulence Factors assignment dataframe"
	echo "	--nmag		INT		Minimum number of Genomes per Library accepted [default=0]"
	echo "	-t			INT		number of threads"
	echo "	-o			STR		output directory"
	echo "	-h --help			print this message"
	echo ""
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
out="false"; checkm="false"; gene="false"; gtdbtk="false"; meta="false"; vf="false";nmag=0

# load in params
OPTS=`getopt -o ht:o: --long help,checkm:,gene:,gtdbtk:,meta:,vf:,nmag: -- "$@"`
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
				--nmag) nmag=$2; shift 2;;
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
echo '-------> START MODULE mSpread'
conda activate mSpreadComp_env
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

Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/simple_description_norm.r --gtdb $gtdbtk --checkm $checkm --gene $gene --meta $meta --nmag_filter $nmag --out $initial_processing_path

echo "Genome Quality description and Gene normalization finished!"


### Perform Beta dispersion analysis ####



### Perform pairwise analysis irrespective o taxa ####



### Perform ARG prevelance analysis and spread per gOTU ####



### Perform NMI figures ####




### Perform Plasmid description ####



### Perform Pathogens description ####



### Perform Network Analysis ####




