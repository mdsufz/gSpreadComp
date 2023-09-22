#!/bin/bash

help_message () {
	echo ""
	echo "Usage: mspreadcomp quality [options] --genome_dir genome_folder -o output_dir"
	echo "Options:"
	echo ""
	echo "	--genome_dir STR	folder with the genomes to estimate quality (in fasta format)"
	echo "	--extension STR		fasta file extension (e.g. fa or fasta) [default: fa]"
	echo "	-o STR				output directory"
	echo "	-t INT      		number of threads [default: 1]"
	echo "	-h --help			print this message"
	echo ""
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
threads=1; out="false"; genome_dir="false"; extension=fa

# load in params
OPTS=`getopt -o ht:o: --long help,genome_dir:,extension: -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                --genome_dir) genome_dir=$2; shift 2;;
				--extension) extension=$2; shift 2;;
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
echo '------- START MODULE Quality Estimation'
conda activate mSpreadComp_env
config_path="$(which config)"
database="${config_path/config/database}"
source $config_path
source $database

# load checkm env
conda activate "$mSPREAD_DEPENDENCIES_ENVS_PATH"/checkm_env

#Set CheckM DB
CHECKM_DB="$DATABASES_LOCATION"checkm

echo ${CHECKM_DB} | checkm data setRoot ${CHECKM_DB}


# check if all parameters are entered
if [ "$out" = "false" ] || [ "$genome_dir" = "false" ]; then 
	help_message; exit 1
fi

if [ -z "$CHECKM_DB" ]; then 
	echo "No CHECKM database found."
	echo "Please make sure you installed the CHECKM database and configured its path"
	echo "You can follow the instructions on the mSpreadComp Github page"
	help_message; exit 1
fi

echo "Your CHECKM_DB database is at $CHECKM_DB"


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

#Run
if [ -f "$out"/outputcheckm.tsv ];
then echo "-> Quality check already done. Please check: "$out"/outputcheckm.tsv"
else
checkm lineage_wf -t $threads --reduced_tree --tab_table -x $extension -f "$out"/outputcheckm.tsv $genome_dir $out
fi

conda deactivate

#Format output
echo "Formatting Taxonomy CheckM results:"

Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/checkm_format.r --checkm "$out"/outputcheckm.tsv -o $out
