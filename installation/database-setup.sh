#!/bin/bash

VERSION=1.0

help_message () {
        echo""
        echo "gspreadcomp database script v=$VERSION"
        echo "Usage: bash -i database-setup.sh --dbs [module] -o output_folder_for_dbs"
		echo "USE THE SAME DATABASE LOCATION OUTPUT FOLDER FOR ALL DATABASES USED WITH gspreadcomp"
        echo ""
        echo "  --dbs all			download and install the required and optional databases [default]"
        echo "  --dbs required              	download and install the required databases (Victors and VFDB) for gspreadcomp"
        echo "  --dbs optional              	download and install all the optional (ARGs, GTDB-tk, CheckM) databases for gspreadcomp"
        echo "  --dbs args			download and install the required and the ARGs databases."
        echo "  -o path/folder/to/save/dbs	output folder where you want to save the downloaded databases"
        echo "  --help | -h			show this help message"
        echo "  --version | -v			show database install script version"
        echo "";}
  

active_module="all"

while true; do
	case "$1" in
		--dbs) active_module=$2; shift 2;;
		-o) database_location=$2; shift 2;;
		-h | --help) help_message; exit 1; shift 1;;
		-v | --version) echo "$VERSION"; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done


mkdir -p "$database_location"
conda activate gspreadcomp_env
config_file="$(which config)"
source "$config_file"

echo DATABASES_LOCATION="$database_location" > ${config_file/config/database}


if [ "$active_module" = "all" ]; then
    ############################################### ALL ###############################################
	
    ### CheckM
    mkdir -p "$database_location"/checkm
    cd  "$database_location"/checkm
    if [ ! -f selected_marker_sets.tsv ]; then
    echo 'installing checkm database ...'
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xf checkm_data_2015_01_16.tar.gz
    rm -fr checkm_data_2015_01_16.tar.gz

    CHECKM_DB="$database_location"/checkm
    echo CHECKM_DB="$CHECKM_DB" >> "$config_path"
    else echo "-> your CheckM database is ready"
    fi
	cd -


    ### GTDB-tk
    mkdir -p  "$database_location"/"gtdbtk"
    cd "$database_location"/"gtdbtk"
    if [ ! -d release*  ]; then

    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    # Extract the contents (assuming a tar archive)
    tar xzf gtdbtk_data.*
    # Remove the downloaded file
    rm -fr gtdbtk_data.*

    else echo "-> your GTDBtk database is ready"
    fi
	cd -
	
### DeepARG
    mkdir -p  "$database_location"/"deeparg"
    cd "$database_location"/"deeparg"
    if [ ! -f model/v2/metadata_LS.pkl  ]; then
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env
	
	deeparg download_data -o "$database_location"/"deeparg"
	
	conda deactivate

   else echo "-> your DeepARG database is ready"
    fi
	cd -
	
	### Victors
    mkdir -p  "$database_location"/"victors"
    cd "$database_location"/"victors"
    if [ ! -f gen_downloads_protein.php.phr ]; then
	
    wget https://phidias.us/victors/downloads/gen_downloads_protein.php
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	makeblastdb -in gen_downloads_protein.php -dbtype prot
	
	conda deactivate

    else echo "-> your Victors database is ready"
    fi
	cd -
	
	### VFDB
    mkdir -p  "$database_location"/"vfdb"
    cd "$database_location"/"vfdb"
    if [ ! -f VFDB_setB_pro.fas.phr ]; then
	
    wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
	gzip -d VFDB_setB_pro.fas
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	makeblastdb -in VFDB_setB_pro.fas -dbtype prot
	
	conda deactivate

    else echo "-> your vfdb database is ready"
    fi
	cd -


elif [ "$active_module" = "required" ]; then

	############################################### REQUIRED ###############################################

	### Victors
    mkdir -p  "$database_location"/"victors"
    cd "$database_location"/"victors"
    if [ ! -f gen_downloads_protein.php.phr ]; then
	
    wget https://phidias.us/victors/downloads/gen_downloads_protein.php
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	makeblastdb -in gen_downloads_protein.php -dbtype prot
	
	conda deactivate

    else echo "-> your Victors database is ready"
    fi
	cd -
	
	### VFDB
    mkdir -p  "$database_location"/"vfdb"
    cd "$database_location"/"vfdb"
    if [ ! -f VFDB_setB_pro.fas.phr ]; then
	
    wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
	gzip -d VFDB_setB_pro.fas
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	makeblastdb -in VFDB_setB_pro.fas -dbtype prot
	
	conda deactivate

    else echo "-> your vfdb database is ready"
    fi
	cd -

	
 
elif [ "$active_module" = "optional" ]; then

	############################################### OPTIONAL ###############################################

	### CheckM
    mkdir -p "$database_location"/checkm
    cd  "$database_location"/checkm
    if [ ! -f selected_marker_sets.tsv ]; then
    echo 'installing checkm database ...'
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xf checkm_data_2015_01_16.tar.gz
    rm -fr checkm_data_2015_01_16.tar.gz

    CHECKM_DB="$database_location"/checkm
    echo CHECKM_DB="$CHECKM_DB" >> "$config_path"
    else echo "-> your CheckM database is ready"
    fi
	cd -


    ### GTDB-tk
    mkdir -p  "$database_location"/"gtdbtk"
    cd "$database_location"/"gtdbtk"
    if [ ! -d release*  ]; then

    wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
    # Extract the contents (assuming a tar archive)
    tar xzf gtdbtk_data.*
    # Remove the downloaded file
    rm -fr gtdbtk_data.*


    else echo "-> your GTDBtk database is ready"
    fi
	cd -
	
	### DeepARG
    mkdir -p  "$database_location"/"deeparg"
    cd "$database_location"/"deeparg"
    if [ ! -f model/v2/metadata_LS.pkl  ]; then
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env
	
	deeparg download_data -o "$database_location"/"deeparg"
	
	conda deactivate

   else echo "-> your DeepARG database is ready"
    fi
	cd -
	
 
elif [ "$active_module" = "args" ]; then

	############################################### REQUIRED AND ARGs ###############################################
	
	### DeepARG
    mkdir -p  "$database_location"/"deeparg"
    cd "$database_location"/"deeparg"
    if [ ! -f model/v2/metadata_LS.pkl  ]; then
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env
	
	deeparg download_data -o "$database_location"/"deeparg"
	
	conda deactivate

   else echo "-> your DeepARG database is ready"
    fi
	cd -
	
	### Victors
    mkdir -p  "$database_location"/"victors"
    cd "$database_location"/"victors"
    if [ ! -f gen_downloads_protein.php.phr ]; then
	
    wget https://phidias.us/victors/downloads/gen_downloads_protein.php
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	makeblastdb -in gen_downloads_protein.php -dbtype prot
	
	conda deactivate

    else echo "-> your Victors database is ready"
    fi
	cd -
	
	### VFDB
    mkdir -p  "$database_location"/"vfdb"
    cd "$database_location"/"vfdb"
    if [ ! -f VFDB_setB_pro.fas.phr ]; then
	
    wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
	gzip -d VFDB_setB_pro.fas
	
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	makeblastdb -in VFDB_setB_pro.fas -dbtype prot
	
	conda deactivate

    else echo "-> your vfdb database is ready"
    fi
	cd -
	
 
else
        comm "Please select a proper parameter."
        help_message
        exit 1
        
fi
