#!/bin/bash

help_message () {
	echo ""
	echo "Usage: mspreadcomp pathogens [options] --genome_dir genome_folder -o output_dir"
	echo "Options:"
	echo ""
	echo "	--genome_dir STR	folder with the genomes to be aligned againt Virulence factors (in fasta format)"
	echo "	--extension STR		fasta file extension (e.g. fa or fasta) [default: fa]"
	echo "	--evalue NUM		evalue, expect value, threshold as defined by NCBI-BLAST [default: 1e-50]"
	echo "	-t INT      		number of threads"
	echo "	-o STR				output directory"
	echo "	-h --help			print this message"
	echo ""
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Set defaults
out="false"; genome_dir="false"; extension=fa; evalue=1e-50

# load in params
OPTS=`getopt -o ht:o: --long help,genome_dir:,extension:,threshold: -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                --genome_dir) genome_dir=$2; shift 2;;
				--extension) extension=$2; shift 2;;
				--evalue) evalue=$2; shift 2;;
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
echo '-------> START MODULE PATHOGENS'
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

#Run VF count

genomes_path=`realpath $genome_dir`

num_files=$(ls $genomes_path/*.$extension | wc -l)
counter=1

for g in $genomes_path/*.$extension; do
	genome=${g##*/}
	mkdir $out/$genome;
	
	
	blastx -db $VICTORS_DB_PATH \
	-num_threads $threads \
	-query $g \
	-out $out/$genome/victors_${genome/.fa/.out} \
	-evalue $evalue \
	-outfmt "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
	
	blastx -db $VFDB_DB_PATH \
	-num_threads $threads \
	-query $g \
	-out $out/$genome/vfdb_${genome/.fa/.out} \
	-evalue $evalue \
	-outfmt "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
	
	#Add file name column
	#Victors
	name="$(basename -- $out/$genome/victors*.out)"
    sed -i "s/$/\t$name/" $out/$genome/victors*.out
	
	#VFDB
	name="$(basename -- $out/$genome/vfdb*.out)"
    sed -i "s/$/\t$name/" $out/$genome/vfdb*.out
	
	# Update the progress bar
	percent=$((100 * $counter / $num_files))
	printf "\rProgress: [%3d%%] File: %s" $percent $genome
	((counter++))

done

echo ""
echo "VF - Done!"

conda deactivate

#Format output
echo "Formatting Virulence Factors results:"

for f in $out/*/victors*.out; do
    cat $f >> $out/victors_merged.out
done

for f in $out/*/vfdb*.out; do
    cat $f >> $out/vfdb_merged.out
done

#Save DBs headers
#Save DB headers

cat $VICTORS_DB_PATH | grep ">" > $out/victors_db_headers.txt
cat $VFDB_DB_PATH | grep ">" >  $out/vfdb_headers.txt


Rscript $mSPREAD_CONDA_ENVIRONMENT_PATH/bin/vf_format.r --victors $out/victors_merged.out --vfdb $out/vfdb_merged.out --vic_header $out/victors_db_headers.txt --vfdb_header $out/vfdb_headers.txt -o $out



