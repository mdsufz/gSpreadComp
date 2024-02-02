#### LITTLE HELP ####
## All conda environments is created inside
## dependencies/conda/envs directory.
## You can manually activate the environment
## using the following command
## $ conda activate path/to/gspreadcomp_env/dependencies/conda/envs/environment_name

#### LITTLE HELP ####
## All modules installed here will be named as "name_env"
## Ex.: gtdbtk_env

################# INSTALLATION'S PRE-CONFIGURATION #################

## Here is where pre-configuration will 
## be defined. In general, in config.sh we collect 
## some important paths

## Importing config.sh
#echo dirname "$dirname"
#echo '----'

# Check if conda is installed
if command -v conda > /dev/null; then
  echo "Conda is installed."
else
  echo "Conda is not installed. Please install Conda and try again."
  exit 1
fi

##Required to be installed:
# mamba in base environment

# Check if mamba is installed
if conda list | grep -q mamba; then
  echo "Mamba is already installed. Version:"
  mamba --version
else
  echo "Mamba is not installed. Please install it before proceeding"
  echo "For instance, you can use:"
  echo "conda install -c conda-forge mamba"
  echo "while you have your base environment activated"
  exit 1
fi

#### Installation started ####

#Run this from the repository directory
mSPREAD_CONDA_ENVIRONMENT_PATH=$CONDA_PREFIX/envs/gspreadcomp_env

#Save the gspreadcomp_env path to the config file
echo mSPREAD_CONDA_ENVIRONMENT_PATH="$mSPREAD_CONDA_ENVIRONMENT_PATH" > $(dirname $0)/temp
cat $(dirname $0)/temp  $(dirname $0)/.config_std > $(dirname $0)/config
rm -f $(dirname $0)/temp 

source $(dirname $0)/config
source $(dirname $0)/installation_utils.sh

## Moving the installation scripts to gspreadcomp_env's environment
# check if mudoger_env already exists
verify_if_main_env_exists "$mSPREAD_CONDA_ENVIRONMENT_PATH"
if [ $main_there == "yes" ]  # if yes, skip installation
then  
echo "Skip gspreadcomp conda env installation. Moving forward..."
else 
echo "Installing gspreadcomp conda env..."  # if no, move forward

start_pre_configuration

fi


echo "your gspreadcomp's path is $mSPREAD_CONDA_ENVIRONMENT_PATH"

yes | cp -rf $DEPENDENCIES_SCRIPTS_PATH $mSPREAD_CONDA_ENVIRONMENT_PATH

yes | cp -rf $INSTALLATION_SCRIPTS_PATH $mSPREAD_DEPENDENCIES_PATH


source installation/config
source installation/installation_utils.sh

chmod +x $mSPREAD_WORK_SCRIPTS/*
chmod +x $mSPREAD_WORK_MODULES/*
chmod +x $mSPREAD_MASTER_SCRIPT 
chmod +x $(dirname $0)/config

yes | cp -rf $mSPREAD_WORK_SCRIPTS/*  $mSPREAD_CONDA_ENVIRONMENT_PATH/bin
yes | cp -rf $mSPREAD_WORK_MODULES/*  $mSPREAD_CONDA_ENVIRONMENT_PATH/bin
yes | cp -rf $mSPREAD_MASTER_SCRIPT   $mSPREAD_CONDA_ENVIRONMENT_PATH/bin
yes | cp -rf $(dirname $0)/config  $mSPREAD_CONDA_ENVIRONMENT_PATH/bin

################# CHOOSING WICH MODULE TO INSTALL #################

## Giving the user the option of which modules to install

echo -e "\n### WELCOME TO gspreadcomp! ###\n"
echo "Do you want to install the complete gspreadcomp pipeline?"
echo "- Main: gspreadcomp metagenome analysis workflow"
echo "- Accessory: Antimicrobial genes annotation with DeepARG [Current inactive]"
echo "- Accessory: Taxonomical assigning with GTDB-tk"
echo "- Accessory: Genome Quality estimation with CheckM"
echo "[Y/N]"

while :
do
	read choose
	if [ $choose = y -o $choose = Y ];
	then
        echo "------> Installing complete workflow"
		install_main=$choose
		install_args=$choose
		install_taxa=$choose
		install_qual=$choose
		break
	elif [ $choose = n -o $choose = N ]
	then
		echo "Choose what you want to install:"
		break
	else
		echo "Command not found, please, try again"
	fi
done
if [ $choose = n -o $choose = N ];
then
	echo "Do you want to install gspreadcomp main workflow? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_main=$choose
			break
		elif [ $choose = n -o $choose = N ]; then
			echo "Installation of gspreadcomp main workflow denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done
	
	echo "Do you want to install Antimicrobial genes annotation? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_args=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Antimicrobial genes annotation denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done

	echo "Do you want to install Taxonomical assigning? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_taxa=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Taxonomical assigning denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done

	echo "Do you want to install Quality estimation? [Y/N]"
	while :
	do
		read choose
		if [ $choose = y -o $choose = Y ];
		then
			install_qual=$choose
			break
		elif [ $choose = n -o $choose = N ]
		then
			echo "Installation of Quality estimation denied"
			break
		else
			echo "Command not found, please, try again"
		fi
	done

fi

################# INSTALLING CHOSEN MODULES #################


echo -e "\nThe gspreadcomp's installation will begin..\n"


coffe_time


if [ ! -z $install_main ];
then
	echo "-----> installing main gspreadcomp"
	
	################################################################
	## Installing all necessary R and R packs ##
	
	## INSTALL HERE ###
	
	###################
	
	################################################################
	## CREATE ENVIRONMENT AND INSTALLING Plasflow ##
	verify_if_conda_env_exist plasflow_env
	if [ $PRESENT == 'yes' ]
	then :;
	else
	#channels config
	conda config --add channels bioconda
	conda config --add channels conda-forge

	#Env creation
	conda create -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/plasflow_env python=3.5
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/plasflow_env
	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/plasflow_env -c jjhelmus tensorflow=0.10.0rc0
	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/plasflow_env plasflow -c smaegol

	#Deactivate env
	conda deactivate
	fi

	################################################################
	## CREATE ENVIRONMENT AND INSTALLING BLAST ##
	verify_if_conda_env_exist blast_env
	if [ $PRESENT == 'yes' ]
	then :;
	else
	#channels config
	conda config --add channels bioconda

	#Env creation
	conda create -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env
	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/blast_env -c bioconda blast

	#Deactivate env
	conda deactivate
	fi
	
fi

if [ ! -z $install_args ];
then

	echo "-----> installing Antimicrobial genes annotation"
	
	################################################################
	## CREATE ENVIRONMENT AND INSTALLING DeepARG ##
	verify_if_conda_env_exist deeparg_env
	if [ $PRESENT == 'yes' ]
	then :;
	else
	#channels config
	#conda config --add channels bioconda

	#Env creation
	#conda create -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env python=2.7.18
	#conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env
	#mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda diamond==0.9.24
	#pip install deeparg==1.0.2

	conda create -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env python=2.7.18
 	source activate $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env
   	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda diamond==0.9.24
   	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda trimmomatic
   	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda vsearch
   	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda bedtools==2.29.2
   	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda bowtie2==2.3.5.1
   	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda samtools
   	pip install git+https://github.com/gaarangoa/deeparg.git

	#Deactivate env
	conda deactivate
	fi
	
fi

if [ ! -z $install_taxa ];
then

    echo "-----> installing Taxonomical assigning"
	
	################################################################
	## CREATE ENVIRONMENT AND INSTALLING GTDB-TK ##
	verify_if_conda_env_exist gtdbtk_env
	if [ $PRESENT == 'yes' ]
	then :;
	else
	#channels config
	conda config --add channels bioconda

	#Env creation
	conda create -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/gtdbtk_env 
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/gtdbtk_env 
	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/gtdbtk_env -c bioconda gtdbtk
	pip3 install mxnet-mkl==1.6.0 numpy==1.23.1

	#Deactivate env
	conda deactivate
	fi
	
fi

if [ ! -z $install_qual ];
then
    
	echo "-----> installing Quality estimation"
	
	verify_if_conda_env_exist checkm_env
	if [ $PRESENT == 'yes' ]
	then :;
	else
	#channels config
	conda config --add channels bioconda

	#Env creation
	conda create -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/checkm_env python=3.9
	conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/checkm_env

	#Community support all-in-one install
	mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/checkm_env -c bioconda checkm-genome

	#Alternative installation
	#mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/checkm_env numpy matplotlib pysam
	#mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/checkm_env hmmer prodigal pplacer
	#pip3 install checkm-genome

	#Deactivate env
	conda deactivate
	fi
	
fi
