## gSpreadComp: Prokaryotic Quantitative Risk Assessment, gene spread, and plasmid-mediated HGT

<p align="center" width="100%">
	<img width="30%" src="/gspreadcomp_logo_noback.png">
</p>

### Description
`gSpreadComp` is a UNIX-based, modular bioinformatics tool designed for prokaryotic quantitative risk assessment, gene spread analysis, and evaluation of plasmid-mediated horizontal transmission to pathogens within microbial communities. It was designed to use sequence data to enhance the traditional quantitative microbial risk assessment (QMRA) models by offering more detailed and integrated analyses of microbial communities in the metagenomic era.

After installation, the user may want to check a detailed tutorial with example input and output data [here](usage_tutorial.md)

### Objectives and Features
- **Six Integrated Modules**: Offers modules for taxonomy assignment, genome quality estimation, ARG annotation, plasmid/chromosome classification, virulence factor annotation, and in-depth downstream analysis, including target-based gene spread analysis and prokaryotic risk assessment.
- **Weighted Average Prevalence (WAP)**: Employs WAP for calculating the spread of target genes at different taxonomical levels or target groups, enabling refined analyses and interpretations of microbial communities.
- **Reference Pathogen Identification**: Compares genomes to the NCBI pathogens database to identify reference pathogens and to determine the pathogenic risk factors.
- **Pathogenic Risk Assessment**: Utilizes the Technique for Order Preference by Similarity to Ideal Solution (TOPSIS) for quantifying pathogenic risk, considering target genes, virulence, and plasmid transmissibility potential.
- **HTML Reporting**: Culminates in a structured HTML report after the complete downstream analysis, providing users with an overview of the results.

### Modular Approach and Flexibility
`gSpreadComp`’s modular nature enables researchers to use the tool's main analysis and report generation steps independently or to integrate only specific pieces of `gSpreadComp` into their pipelines, providing flexibility and accommodating the varying software management needs of investigators.

### Comprehensive Workflow
`gSpreadComp` amalgamates genome annotation, normalization, and sequence comparison in a unified approach, systematically quantifying gene spread and integrating plasmid-mediated gene transfer annotation with the whole microbiome community in a genome-reference independent manner and furnishing a sophisticated QMRA metric.

![ScreenShot](/test_data/01_Kasmanas_gSpread_Fig_1.png)

# Requirements

Before installing and running `gSpreadComp`, ensure that your system meets the following requirements:

## 1. Operating System
- Linux x64 system

## 2. Package Managers
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html): Required for creating environments and managing packages.
- [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html): A faster package manager used within the `gSpreadComp` installation.

## 3. Storage
- Approximately 15 GB for software installation.
- Around 92 GB for the entire database requirements.

# Installation

## Database Management
`gSpreadComp` includes an easy-to-use script for automatic download and configuration of the required databases, with scheduled updates every January and July.

## Modular Nature
Allows for independent usage of its bioinformatic approaches, facilitating integration into existing pipelines and offering flexibility in analysis.

## Compatibility and Requirements
Designed to support Linux x64 systems, requiring approximately 15 GB for software installation and around 92 GB for the entire database requirements.

## 1 - Install miniconda

To bypass conflicting dependencies, the gSpreadComp approach uses miniconda to create automatically orchestrated environments. [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) is a much faster package manager than conda and is used within the gSpreadComp installation. Consequently, miniconda and mamba are required to be previously installed in your system. Below is a possible way of installing miniconda and mamba. Please, be aware that mamba works best when installed in your base environment.

```console
# See documentation: https://docs.conda.io/en/latest/miniconda.html

$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ chmod +x Miniconda3-latest-Linux-x86_64.sh
$ ./Miniconda3-latest-Linux-x86_64.sh
$ export PATH=~/miniconda3/bin:$PATH

# Install mamba. See documentation: https://mamba.readthedocs.io/en/latest/installation.html
$ conda install mamba -n base -c conda-forge
```

## 2 - Install gSpreadComp

Once you have miniconda and mamba installed and on your PATH, you can proceed to install gSpreadComp.
The installation script was designed to install and set up all necessary tools and packages.

```console
# Clone repository
$ git clone https://github.com/JotaKas/gSpreadComp.git

# Go to the gSpreadComp cloned repository folder
$ cd gSpreadComp

# Make sure you have conda ready and that you are in your base environment.
$ conda activate base
$ echo $CONDA_PREFIX

# You should see something like the following:
/path/to/miniconda3

# Run the installation script as follows
$ bash -i installation/install.sh

# Follow the instructions on the screen:
# Enter "y" if you want to install all modules; otherwise, enter "n".
# If you entered "n", enter "y" for each of the modules you would like to install individually.

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

## 3 - Install necessary databases

**Make sure to run the database setup after gSpreadComp is installed.**

Some bioinformatics tools used within gSpreadComp require specific databases to work. We developed a database download and set up tool to make our lives easier. You can choose to install only the databases you intend to use. You can use the flag `--dbs` to choose and set up the selected databases (all [default], install all databases).

Use this script if you want gSpreadComp to take care of everything.

```console
# Make sure gSpreadComp_env is activated. It should have been created when you ran 'bash -i installation/install.sh'
$ conda activate gspreadcomp_env

# Go to gSpreadComp cloned directory
$ cd gSpreadComp

# Run the database setup script
$ bash -i installation/database-setup.sh --dbs all -o /path/to/save/databases

# You can also check out the database-setup help information
$ bash -i installation/database-setup.sh --help

        gSpreadComp database script v=1.0
        Usage: bash -i database-setup.sh --dbs [module] -o output_folder_for_dbs
		    USE THE SAME DATABASE LOCATION OUTPUT FOLDER FOR ALL DATABASES USED WITH gSpreadCOMP
          --dbs all				download and install the required and optional databases [default]"
          --dbs required              		download and install the required databases (Victors and VFDB) for gSpreadComp
          --dbs optional              		download and install all the optional (ARGs, GTDB-tk, CheckM) databases for gSpreadComp
          --dbs args				download and install the required and the ARGs databases.
          -o path/folder/to/save/dbs		output folder where you want to save the downloaded databases
          --help | -h				show this help message
          --version | -v			show database install script version


```

## Usage

### Activating the Conda Environment
Before using `gSpreadComp`, activate the appropriate conda environment using the following command:
```sh
conda activate gSpreadComp_env
```

### Command-Line Usage
`gSpreadComp` provides several modules, each performing a specific task within the pipeline. The quick command-line usage is as follows:
```sh
gspreadcomp --help
```

### Modules and Their Descriptions
`gSpreadComp` comprises several modules, each serving a specific purpose in the genome analysis workflow:

#### 1. Taxonomy Assignment
```sh
gspreadcomp taxonomy [options] --genome_dir genome_folder -o output_dir
```
- Assigns taxonomy to genomes using [GTDBtk v2](https://academic.oup.com/bioinformatics/article/38/23/5315/6758240).
- Options:
  - `--genome_dir STR`: folder with the bins to be classified (in fasta format)
  - `--extension STR`: fasta file extension (e.g. fa or fasta) [default: fa]
  - `-o STR`: output directory
  - `-t INT`: number of threads

#### 2. Genome Quality Estimation
```sh
gspreadcomp quality [options] --genome_dir genome_folder -o output_dir
```
- Estimates genome completeness and contamination using [CheckM](https://genome.cshlp.org/content/25/7/1043).
- Options:
  - `--genome_dir STR`: folder with the genomes to estimate quality (in fasta format)
  - `--extension STR`: fasta file extension (e.g. fa or fasta) [default: fa]
  - `-o STR`: output directory
  - `-t INT`: number of threads [default: 1]
  - `-h --help`: print this message

#### 3. ARG Prediction
```sh
gspreadcomp args [options] --genome_dir genome_folder -o output_dir
```
- Predicts the Antimicrobial Resistance Genes (ARGs) in a genome using [DeepARG](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0401-z).
- Options:
  - `--genome_dir STR`: folder with the genomes to be classified (in fasta format)
  - `--extension STR`: fasta file extension (e.g. fa or fasta) [default: fa]
  - `--min_prob NUM`: Minimum probability cutoff for DeepARG [Default: 0.8]
  - `--arg_alignment_identity NUM`: Identity cutoff for sequence alignment for DeepARG [Default: 35]
  - `--arg_alignment_evalue NUM`: Evalue cutoff for DeepARG [Default: 1e-10]
  - `--arg_alignment_overlap NUM`: Alignment read overlap for DeepARG [Default: 0.8]
  - `--arg_num_alignments_per_entry NUM`: Diamond, minimum number of alignments per entry [Default: 1000]
  - `-o STR`: output directory
  - `-h --help`: print this message

#### 4. Plasmid Prediction
```sh
gspreadcomp plasmid [options] --genome_dir genome_folder -o output_dir
```
- Predicts if a sequence within a fasta file is a chromosome, plasmid, or undetermined using [Plasflow](https://academic.oup.com/nar/article/46/6/e35/4807335).
- Options:
  - `--genome_dir STR`: folder with the genomes to be classified (in fasta format)
  - `--extension STR`: fasta file extension (e.g. fa or fasta) [default: fa]
  - `--threshold NUM`: threshold for probability filtering [default: 0.7]
  - `-o STR`: output directory
  - `-h --help`: print this message

#### 5. Pathogen Alignment
```sh
gspreadcomp pathogens [options] --genome_dir genome_folder -o output_dir
```
- Aligns provided genomes to Virulence Factors databases and formats the output.
- Options:
  - `--genome_dir STR`: folder with the genomes to be aligned against Virulence factors (in fasta format)
  - `--extension STR`: fasta file extension (e.g. fa or fasta) [default: fa]
  - `--evalue NUM`: evalue, expect value, threshold as defined by NCBI-BLAST [default: 1e-50]
  - `-t INT`: number of threads
  - `-o STR`: output directory
  - `-h --help`: print this message

#### 6. Main Analysis
```sh
gspreadcomp gspread [options] -o output_dir
```
- Runs the main `gSpreadComp` to compare spread and plasmid-mediated HGT.
- Options:
  - `--checkm STR`: Path to the formatted Quality estimation dataframe
  - `--gene STR`: Path to the formatted target Gene dataframe to calculate the spread
  - `--gtdbtk STR`: Path to the formatted Taxonomy assignment dataframe
  - `--meta STR`: Path to the formatted Sample's Metadata dataframe
  - `--vf STR`: Path to the formatted Virulence Factors assignment dataframe
  - `--plasmid STR`: Path to the formatted Plasmid prediction dataframe
  - `--nmag INT`: Minimum number of Genomes per Library accepted [default=0]
  - `--spread_taxa STR`: Taxonomic level to check gene spread [default=Phylum]
  - `--target_gene_col STR`: Name of the column from the gene dataset with the Gene_ids to analyse [default=Gene_id]
  - `-t INT`: number of threads
  - `-o STR`: output directory
  - `-h --help`: print this message




