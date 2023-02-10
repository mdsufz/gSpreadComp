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
	echo "	-t			INT		number of threads"
	echo "	-o			STR		output directory"
	echo "	-h --help			print this message"
	echo ""
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
out="false"; checkm="false"; gene="false"; gtdbtk="false"; meta="false"; vf="false"

# load in params
OPTS=`getopt -o ht:o: --long help,checkm:,gene:,gtdbtk:,meta:,vf: -- "$@"`
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

# load blast env
conda activate "$mSPREAD_DEPENDENCIES_ENVS_PATH"/blast_env

# Load databases paths
VICTORS_DB_PATH="$DATABASES_LOCATION"victors/gen_downloads_protein.php
VFDB_DB_PATH="$DATABASES_LOCATION"vfdb/VFDB_setB_pro.fas

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$genome_dir" = "false" ]; then 
	help_message; exit 1
fi

if [ -z "$VICTORS_DB_PATH" ]; then 
	echo "No Victors database found."
	echo "Please make sure you installed the Victors database and configured its path"
	echo "You can follow the instructions on the mSpreadComp Github page"
	help_message; exit 1
fi

if [ -z "$VFDB_DB_PATH" ]; then 
	echo "No VFDB database found."
	echo "Please make sure you installed the VFDB database and configured its path"
	echo "You can follow the instructions on the mSpreadComp Github page"
	help_message; exit 1
fi

echo "Your Victors database is at $VICTORS_DB_PATH"
echo "Your VFDB database is at $VFDB_DB_PATH"

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

### Normalize and Describe sample only Mags vs Target meta ####



### Perform Complete-Contamination analysis ####



### Perform Beta dispersion analysis ####



### Perform pairwise analysis irrespective o taxa ####



### Perform ARG prevelance analysis and spread per gOTU ####



### Perform NMI figures ####




### Perform Plasmid description ####



### Perform Pathogens description ####



### Perform Network Analysis ####




