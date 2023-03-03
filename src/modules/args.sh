#!/bin/bash

help_message () {
	echo ""
	echo "Usage: mspreadcomp  [options] --genome_dir genome_folder -o output_dir"
	echo "Options:"
	echo ""
	echo "	--genome_dir STR	folder with the genomes to be classified (in fasta format)"
	echo "	--extension STR		fasta file extension (e.g. fa or fasta) [default: fa]"
	echo "	--min_prob NUM		Minimum probability cutoff for DeepARG [Default: 0.8]"
	echo "	--arg_alignment_identity NUM		Identity cutoff for sequence alignment for DeepARG [Default: 35]"
	echo "	--arg_alignment_evalue NUM		Evalue cutoff for DeepARG [Default: 1e-10]"
	echo "	--arg_alignment_overlap NUM		Alignment read overlap for DeepARG [Default: 0.8]"
	echo "	--arg_num_alignments_per_entry NUM		Diamond, minimum number of alignments per entry [Default: 1000]"
	echo "	-o STR				output directory"
	echo "	-h --help			print this message"
	echo ""
	echo "";}
	

########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
out="false"; genome_dir="false"; extension=fa; min_prob=0.8; arg_alignment_identity=35 ;arg_alignment_evalue=1e-10; arg_alignment_overlap=0.8; arg_num_alignments_per_entry=1000;

# load in params
OPTS=`getopt -o ho: --long help,genome_dir:,extension:,min_prob:arg_alignment_identity:,arg_alignment_evalue:,arg_alignment_overlap:,arg_num_alignments_per_entry: -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                --genome_dir) genome_dir=$2; shift 2;;
				--extension) extension=$2; shift 2;;
				--min_prob) min_prob=$2; shift 2;;
				--arg_alignment_identity) arg_alignment_identity=$2; shift 2;;
				--arg_alignment_evalue) arg_alignment_evalue=$2; shift 2;;
				--arg_alignment_overlap) arg_alignment_overlap=$2; shift 2;;
				--arg_num_alignments_per_entry) arg_num_alignments_per_entry=$2; shift 2;;
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
echo '------- START MODULE ARGs Prediction'
conda activate mSpreadComp_env
config_path="$(which config)"
database="${config_path/config/database}"
source $config_path
source $database

# load deeparg env
conda activate "$mSPREAD_DEPENDENCIES_ENVS_PATH"/deeparg_env



# Load databases paths
DEEPARG_DB_PATH="$DATABASES_LOCATION"deeparg/database/v2/features.fasta

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$genome_dir" = "false" ]; then 
	help_message; exit 1
fi


if [ -z "$DEEPARG_DB_PATH" ]; then 
	echo "No DeepARG database found."
	echo "Please make sure you installed the DeepARG database and configured its path"
	echo "You can follow the instructions on the mSpreadComp Github page"
	echo "In addition, you can check if the the DeepARG host repository is online"
	help_message; exit 1
fi

echo "Your DeepARG database is at $DEEPARG_DB_PATH"
########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

#Run DeepARG

DEEPARG_DB_PATH="$DATABASES_LOCATION"deeparg/

genomes_path=`realpath $genome_dir`

for g in $genomes_path/*.$extension; do
	genome=${g##*/}
	mkdir $out/$genome;
	deeparg predict --model LS \
	-i $g \
	-o $out/$genome/"$genome"_deeparg_out \
	-d $DEEPARG_DB_PATH \
	--type nucl \
	--min-prob $min_prob \
	--arg-alignment-identity $arg_alignment_identity \
	--arg-alignment-evalue $arg_alignment_evalue \
	--arg-alignment-overlap $arg_alignment_overlap \
	--arg-num-alignments-per-entry $arg_num_alignments_per_entry \
	--model-version "v2"
done

conda deactivate

#Format output
echo "Formatting ARGs annotation results: NOT FINISHED!!!!"

#Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/deeparg_format.r -i $out -o $out

