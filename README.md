# mSpreadComp
a novel automatic approach to compare the spread and plasmid-mediated horizontal transmission of annotated genetic potential to pathogens in microbial communities


# Installation

**1 - Install miniconda**

To bypass conflicting dependencies, the mSpreadComp approach uses miniconda to create automatically orchestrated environments. In addition, [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) is a much faster package manager than conda and is used within the mSpreadComp installation. Consequently, miniconda and mamba are required to be previously installed in your system. Following you have a possible way of installing miniconda and mamba. Please, be aware that mamba works best when installed in your base environment.

```console

#See documentation: https://docs.conda.io/en/latest/miniconda.html

$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ chmod +x Miniconda3-latest-Linux-x86_64.sh

$ ./Miniconda3-latest-Linux-x86_64.sh

$ export PATH=~/miniconda3/bin:$PATH

#Install mamba. See documentation: https://mamba.readthedocs.io/en/latest/installation.html

$ conda install mamba -n base -c conda-forge

```

**2 - Install mSpreadComp**

Once you have miniconda and mamba installed and on your PATH, you can proceed to install mSpreadComp.
The installation script was designed to install and set up all necessary tools and packages.

```console
#clone repository

$ git clone https://github.com/JotaKas/mSpreadComp.git

#Go to the mSpreadComp cloned repository folder
$ cd mSpreadComp

#Make sure you have conda ready and that you are in your base environment.
$ conda activate base
$ echo $CONDA_PREFIX

#You should see something like the following:
/path/to/miniconda3

#Run the installation script as follows
$ bash -i installation/install.sh

#Follow the instructions on the screen:
# Enter "y" if you want to install all modules, otherwise enter "n".
# If you entered "n",enter "y" for each of the modules you would like to install individually.

	The MuDoGeR's installation will begin..





	      (  )   (   )  )			
	       ) (   )  (  (			
	       ( )  (    ) )			
	       _____________			
	      <_____________> ___		
	      |             |/ _ \		
	      |               | | |		
	      |               |_| |		
	   ___|             |\___/		
	  /    \___________/    \		
	  \_____________________/		

	This might take a while. Time to grab a coffee...
```

**3 - Install necessary databases**

**Make sure to run the database setup after MuDoGeR is installed.**

Some bioinformatics tools used within mSpreadComp require specific databases to work. We developed a database download and set up tool to make our lives easier. 

You can choose to install only the databases you intend to use. You can use the flag ```--dbs``` to choose and set up the selected databases (all \[default], install all databases).

Use this script if you want mSpreadComp to take care of everything. 

```console
#Make sure mSpreadComp_env is activated. It should have been created when you ran 'bash -i installation/install.sh'
$ conda activate mSpreadComp_env

#Go to mSpreadComp cloned directory
$ cd mSpreadComp

#Run the database setup script
$ bash -i installation/database-setup.sh --dbs all -o /path/to/save/databases

#You can also check out the database-setup help information
$ bash -i installation/database-setup.sh --help

        mSpreadComp database script v=1.0
        Usage: bash -i database-setup.sh --dbs [module] -o output_folder_for_dbs
		    USE THE SAME DATABASE LOCATION OUTPUT FOLDER FOR ALL DATABASES USED WITH MSPREADCOMP
          --dbs all				download and install the required and optional databases [default]"
          --dbs required              		download and install the required databases (Victors and VFDB) for mSpreadComp
          --dbs optional              		download and install all the optional (ARGs, GTDB-tk, CheckM) databases for mSpreadComp
          --dbs args				download and install the required and the ARGs databases.
          -o path/folder/to/save/dbs		output folder where you want to save the downloaded databases
          --help | -h				show this help message
          --version | -v			show database install script version


```
