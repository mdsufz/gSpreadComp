#Have mamba installed on base before

#Install databases if want to use the tools
#DeepARG
# deeparg download_data -o /path/to/local/directory/

#GTDB-TK

#CheckM

#VFDB and Victors


################# Installing gspreadcomp #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - R and relevant R packs
## - Plasflow
## - DeepARG
## - GTDB-TK
## - CheckM
## - BLAST
## - Victors and VFDB databases

echo "### gspreadcomp installation started ###\n"
source installation/config           
source installation/installation_utils.sh

## Checking if some tool already have a conda environment created

################################################################
## CREATE ENVIRONMENT AND INSTALLING CheckM ##
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
## CREATE ENVIRONMENT AND INSTALLING DeepARG ##
verify_if_conda_env_exist deeparg_env
if [ $PRESENT == 'yes' ]
then :;
else
#channels config
conda config --add channels bioconda

#Env creation
conda create -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env python=2.7.18
conda activate $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env
mamba install -y --prefix $mSPREAD_DEPENDENCIES_ENVS_PATH/deeparg_env -c bioconda diamond==0.9.24
pip install deeparg==1.0.2

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
