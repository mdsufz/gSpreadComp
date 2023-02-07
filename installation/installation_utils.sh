source installation/config

#### DESCRIPTION #### 
## Call a given installation script from the instalation script path
##
#### PARAMS ####
## $1 the script name
##
#### OUTPUT #### 
## The bash command output to the given script passed as a parameter
call_installation_script() {
    command="$1"
    for i; do
        bash -i "$mSPREAD_INSTALLATION_SCRIPTS_PATH/${command}.sh" "${arg}"
    done
}

#### DESCRIPTION #### 
## Create a environment to mudoger and the paths that will be used
## to store it's dependencies
##
#### OUTPUT #### 
## The conda and mkdir outputs
start_pre_configuration() {

    conda create -y -n mSpreadComp_env r-essentials r-base
    conda activate mSpreadComp_env
    mkdir -p $mSPREAD_CONDA_ENVIRONMENT_PATH/dependencies
    mkdir -p $mSPREAD_CONDA_ENVIRONMENT_PATH/dependencies/conda
    mkdir -p $mSPREAD_CONDA_ENVIRONMENT_PATH/dependencies/conda/envs
    mkdir -p $mSPREAD_CONDA_ENVIRONMENT_PATH/dependencies/cloned_tools
    mkdir -p $mSPREAD_CONDA_ENVIRONMENT_PATH/dependencies/installation_scripts
    
    conda deactivate
}

#### DESCRIPTION #### 
## Verify if a given conda environment exist
##
#### PARAMS ####
## $1 name of environment
##
#### OUTPUT #### 
## if environment's path exist, display a message and pops an error
## otherwise, display the message "$env_path doesn't exist yet"
verify_if_conda_env_exist() {
    env_path=$mSPREAD_DEPENDENCIES_ENVS_PATH/$1
    if [[ -d $env_path ]] && [[ ! -z $1 ]];
    then
        echo "-> Environment $1 already exist in path: $env_path"
        PRESENT='yes'
    else
        PRESENT='no'
        #echo "$env_path doesn't exist yet"
    fi
}


#### DESCRIPTION #### 
## Verify if the main conda environment exists
##
#### PARAMS ####
## $1 name of environment
##
#### OUTPUT #### 
## if environment's path exist, display a message and pops an error
## otherwise, display the message "$env_path doesn't exist yet"
verify_if_main_env_exists() {
    env_path=$mSPREAD_DEPENDENCIES_ENVS_PATH
    if [[ -d $env_path ]] && [[ ! -z $1 ]];
    then
        echo "Environment $1 already exists". # in path: $env_path"
        #exit 1
        main_there='yes'
    else
        echo "$env_path doesn't exist yet"
        main_there='no'
    fi
}

coffe_time(){
echo -e '\n\n'
echo '								'
echo '      (  )   (   )  )			'
echo '       ) (   )  (  (			'	
echo '       ( )  (    ) )			'
echo '       _____________			'
echo '      <_____________> ___		'
echo '      |             |/ _ \		'
echo '      |               | | |		'	
echo '      |               |_| |		'
echo '   ___|             |\___/		'
echo '  /    \___________/    \		'	
echo '  \_____________________/		'

echo -e '\n	This might take a while. Time to grab a coffee...\n'
sleep 3
}

