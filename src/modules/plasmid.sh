#!/bin/bash

help_message () {
	echo ""
	echo "Usage: mspreadcomp plasmid [options] --genome_dir genome_folder -o output_dir"
	echo "Options:"
	echo ""
	echo "	--genome_dir STR	folder with the genomes to be classified (in fasta format)"
	echo "	--extension STR		fasta file extension (e.g. fa or fasta) [default: fa]"
	echo "	--threshold NUM		threshold for probability filtering [default: 0.7]"
	echo "	-o STR				output directory"
	echo "	-h --help			print this message"
	echo ""
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
out="false"; genome_dir="false"; extension=fa; threshold=0.7

# load in params
OPTS=`getopt -o ho: --long help,genome_dir:,extension:,threshold: -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                --genome_dir) genome_dir=$2; shift 2;;
				--extension) extension=$2; shift 2;;
				--threshold) threshold=$2; shift 2;;
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
echo '------- START MODULE Plasmid Prediction'
conda activate mSpreadComp_env
config_path="$(which config)"
database="${config_path/config/database}"
source $config_path
source $database

# load checkm env
conda activate "$mSPREAD_DEPENDENCIES_ENVS_PATH"/plasflow_env

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$genome_dir" = "false" ]; then 
	help_message; exit 1
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

#Run PlasFlow

genomes_path=`realpath $genome_dir`

for g in $genomes_path/*.$extension; do
	genome=${g##*/}
	mkdir $out/$genome;
	PlasFlow.py --input $g --output $out/$genome/"$genome"_plasflow_out.tsv --threshold $threshold;
done

conda deactivate

#Format output
echo "Formatting Plasflow results:"

Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/plasflow_format.r -i $out -o $out

