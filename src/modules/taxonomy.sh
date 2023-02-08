#!/bin/bash

help_message () {
	echo ""
	echo "Usage: mspreadcomp taxonomy [options] --genome_dir genome_folder -o output_dir"
	echo "Options:"
	echo ""
	echo "	--genome_dir STR	folder with the bins to be classified (in fasta format)"
	echo "	--extension STR		fasta file extension (e.g. fa or fasta) [default: fa]"
	echo "	-o STR				output directory"
	echo "	-t INT      		number of threads"
	echo ""
	echo "";}

########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
threads=1; out="false"; genome_dir="false"; extension=fa

# load in params
OPTS=`getopt -o ht:o: --long help -- "$@"`
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
echo '------- START MODULE TAXONOMY'
conda activate mSpreadComp_env
config_path="$(which config)"
database="${config_path/config/database}"
source $config_path
source $database

# load gtdbtk env
conda activate "$mSPREAD_DEPENDENCIES_ENVS_PATH"/gtdbtk_env

GTDBTK_DATA_PATH=$(realpath "$DATABASES_LOCATION"gtdbtk/release*)

export GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$genome_dir" = "false" ]; then 
	help_message; exit 1
fi

if [ -z "$GTDBTK_DATA_PATH" ]; then 
	echo "No GTDBtk database found."
	echo "Please make sure you installed the GTDBtk database and configured its path"
	echo "You can follow the instructions on the mSpreadComp Github page"
	help_message; exit 1
fi

echo "Your GTDBtk database is at $GTDBTK_DATA_PATH"


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

#Run
if [ -f "$out"/gtdbtk.bac120.summary.tsv ] || [ -f "$out"/gtdbtk.ar53.summary.tsv ];
then echo "";
else
gtdbtk classify_wf --extension $extension --cpus "$threads" --genome_dir $genome_dir --out_dir $out
fi

#Create merged results file
if [ -f "$out"/gtdbtk.bac120.summary.tsv ] || [ -f "$out"/gtdbtk.ar53.summary.tsv ]; 
then awk 'FNR==1 && NR!=1 {next;}{print}' "$out"/gtdbtk.*summ*.tsv > "$out"/gtdbtk_result.tsv;
echo "GTDBtk results generated!"
else
echo "Error: GTDBtk summary files not found"
fi

conda deactivate