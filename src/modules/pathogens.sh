#!/bin/bash

help_message () {
	echo ""
	echo "Usage: mspreadcomp pathogens [options] --genome_dir genome_folder -o output_dir"
	echo "Options:"
	echo ""
	echo "	--genome_dir STR	folder with the genomes to be aligned againt Virulence factors (in fasta format)"
	echo "	--extension STR		fasta file extension (e.g. fa or fasta) [default: fa]"
	echo "	--evalue NUM		evalue, expect value, threshold as defined by NCBI-BLAST [default: 1e-50]"
	echo "	--vf STR		select the virulence factors database to be used (e.g. victors, vfdb or both) [default: both]"
	echo "	-t INT      		number of threads"
	echo "	-o STR				output directory"
	echo "	-h --help			print this message"
	echo ""
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
out="false"; genome_dir="false"; extension=fa; evalue=1e-50; t=1; vf=both

# load in params
OPTS=`getopt -o ht:o: --long help,genome_dir:,extension:,evalue:,vf: -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                --genome_dir) genome_dir=$2; shift 2;;
		--extension) extension=$2; shift 2;;
		--evalue) evalue=$2; shift 2;;
		-t) threads=$2; shift 2;;
		--vf) vf=$2; shift 2;;	
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
echo '-------> START MODULE PATHOGENS'
conda activate mSpreadComp_env
config_path="$(which config)"
database="${config_path/config/database}"
source $config_path
source $database

# load blast env
conda activate "$mSPREAD_DEPENDENCIES_ENVS_PATH"/blast_env


# Check if the user has given valid values for --database
if [[ "$vf" != "victors" && "$vf" != "vfdb" && "$vf" != "both" ]]; then
    echo "Error: Invalid value for --vf. Use 'victors', 'vfdb', or 'both'."
    help_message
    exit 1
fi

# check if all parameters are entered
# Get the real paths and print them for troubleshooting
if [ "$out" != "false" ]; then 
    out=$(realpath "$out" 2>/dev/null) 
    if [ $? -ne 0 ]; then 
        echo "Error: Invalid output directory path provided: $out"
        exit 1
    fi
    echo "Using output directory: $out"
fi

if [ "$genome_dir" != "false" ]; then 
    genome_dir=$(realpath "$genome_dir" 2>/dev/null) 
    if [ $? -ne 0 ]; then 
        echo "Error: Invalid genome directory path provided: $genome_dir"
        exit 1
    fi
    echo "Using genome directory: $genome_dir"
fi

# check if all parameters are entered
if [ "$out" = "" ] || [ "$genome_dir" = "" ]; then 
    echo "Error: Output directory and Genome directory must be specified."
    help_message; exit 1
fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

if [[ "$vf" == "victors" || "$vf" == "both" ]]; then
    # Commands and logic related to VICTORS_DB_PATH
    
    echo "Running pathogens module for Victors DB"
    # Load databases paths
    VICTORS_DB_PATH="$DATABASES_LOCATION"victors/gen_downloads_protein.php
    
    if [ -z "$VICTORS_DB_PATH" ]; then 
	echo "No Victors database found."
	echo "Please make sure you installed the Victors database and configured its path"
	echo "You can follow the instructions on the mSpreadComp Github page"
	help_message; exit 1
	fi
	
    echo "Your Victors database is at $VICTORS_DB_PATH"
    
    num_files=$(ls $genome_dir/*.$extension | wc -l)
    counter=1

    for g in $genome_dir/*.$extension; do
	genome=${g##*/}
	mkdir $out/$genome 2>/dev/null;
	
	blastx -db $VICTORS_DB_PATH \
	-num_threads $threads \
	-query $g \
	-out $out/$genome/victors_${genome/.fa/.out} \
	-evalue $evalue \
	-outfmt "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
		
	#Add file name column
	#Victors
	name="$(basename -- $out/$genome/victors*.out)"
	sed -i "s/$/\t$name/" $out/$genome/victors*.out

	# Update the progress bar
	percent=$((100 * $counter / $num_files))
	printf "\rProgress: [%3d%%] File: %s" $percent $genome
	((counter++))

    done
    echo ""
    
    for f in $out/*/victors*.out; do
        cat $f >> $out/victors_merged.out
    done
    
    #Save DB headers
    cat $VICTORS_DB_PATH | grep ">" > $out/victors_db_headers.txt
    
    echo "Victors - Done!"

    
fi

if [[ "$vf" == "vfdb" || "$vf" == "both" ]]; then
    # Commands and logic related to VFDB_DB_PATH
    
    echo "Running pathogens module for VFDB"
    
    # Load databases paths
    VFDB_DB_PATH="$DATABASES_LOCATION"vfdb/VFDB_setB_pro.fas
    
    if [ -z "$VFDB_DB_PATH" ]; then 
	echo "No VFDB database found."
	echo "Please make sure you installed the VFDB database and configured its path"
	echo "You can follow the instructions on the mSpreadComp Github page"
	help_message; exit 1
	fi
    echo "Your VFDB database is at $VFDB_DB_PATH"
    
    #begin
    num_files=$(ls $genome_dir/*.$extension | wc -l)
    counter=1
    
    for g in $genome_dir/*.$extension; do
	genome=${g##*/}
	mkdir $out/$genome 2>/dev/null;
	
	blastx -db $VFDB_DB_PATH \
	-num_threads $threads \
	-query $g \
	-out $out/$genome/vfdb_${genome/.fa/.out} \
	-evalue $evalue \
	-outfmt "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
	
	#Add file name column
	
	#VFDB
	name="$(basename -- $out/$genome/vfdb*.out)"
	sed -i "s/$/\t$name/" $out/$genome/vfdb*.out
	
	# Update the progress bar
	percent=$((100 * $counter / $num_files))
	printf "\rProgress: [%3d%%] File: %s" $percent $genome
	((counter++))
    
    done
    echo ""
    
    for f in $out/*/vfdb*.out; do
        cat $f >> $out/vfdb_merged.out
    done
    
    #Save DBs headers
    cat $VFDB_DB_PATH | grep ">" >  $out/vfdb_headers.txt
    
    echo "VFDB - Done!"    
    
fi


conda deactivate

#Format output
echo "Formatting Virulence Factors results:"

if [[ "$vf" == "victors" || "$vf" == "both" ]]; then
    
    Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/vf_format.r --vf $out/victors_merged.out --vf_header $out/victors_db_headers.txt --e_value $evalue --vf_prefix "victors" -o $out
    
fi

if [[ "$vf" == "vfdb" || "$vf" == "both" ]]; then
    
    Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/vf_format.r --vf $out/vfdb_merged.out --vf_header $out/vfdb_headers.txt --e_value $evalue --vf_prefix "vfdb" -o $out
    
fi
