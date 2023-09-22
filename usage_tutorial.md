# gSpreadComp Usage Tutorial

Welcome to the `gSpreadComp` usage tutorial! This guide is designed to help you understand how to use the tool effectively through practical examples and step-by-step instructions.

## Introduction

`gSpreadComp` is designed to work with genomes in fasta format and requires a metadata table in CSV format containing a target variable. For example tables and more details, please refer to [example_metadata_table_link] and [example_genome_fasta_link]. For the recovery of prokaryotic genomes from metagenomic samples, you can use tools like [MuDoGeR](https://github.com/mdsufz/MuDoGeR).

### Prerequisites

Before proceeding with the tutorial, make sure you have:
1. Genomes in fasta format.
2. A metadata table in CSV format with a target variable.

### Modules and Steps

In this tutorial, we will go through the following modules and steps of `gSpreadComp`:
1. **Taxonomy Assignment**
2. **Genome Quality Estimation**
3. **ARG Annotation**
4. **Plasmid Identification**
5. **Virulence Factors (VFs) Annotation**
6. **Downstream Analysis**

We will inspect the main outputs of each module to ensure a comprehensive understanding of the tool's functionalities and results.

Let's get started!

## Taxonomy Assignment Module

The taxonomy module in `gSpreadComp` uses GTDBtk for taxonomy assignment. To run this module, use the `gspreadcomp taxonomy` command. Below are the available options for this module:

```console
gspreadcomp taxonomy --help

Usage: gspreadcomp taxonomy [options] --genome_dir genome_folder -o output_dir
Options:
    --genome_dir STR    folder with the bins to be classified (in fasta format)
    --extension STR     fasta file extension (e.g. fa or fasta) [default: fa]
    -o STR              output directory
    -t INT              number of threads
```

### Running the Taxonomy Module

1. Create a folder for the test run, for example, `test_gspread_run`.
2. Place your genomes in a subfolder within the test run folder, for example, `01_input_genomes`.
3. Create an output folder within the test run folder, for example, `03_gspread_gtdb_taxonomy`.

Assuming you have placed your genomes in `01_input_genomes`, your fasta files have the `fa` extension, and your output folder is `03_gspread_gtdb_taxonomy`, your command will look like:

```console
$ gspreadcomp taxonomy --genome_dir ./01_input_genomes/ --extension fa -o ./03_gspread_gtdb_taxonomy/ -t 25
```

### Output of the Taxonomy Module

After running the Taxonomy Module, you will find the output in the specified output directory, structured as follows:

```
03_gspread_gtdb_taxonomy/
├── align
├── classify
├── gtdb_df_format_gSpread.csv
├── gtdbtk.bac120.summary.tsv -> classify/gtdbtk.bac120.summary.tsv
├── gtdbtk.log
├── gtdbtk_result.tsv
├── gtdbtk.warnings.log
└── identify
```

#### Main Output File: gtdb_df_format_gSpread.csv. Understanding the Output

The `gtdb_df_format_gSpread.csv` file contains taxonomy information for each genome in a format that is ready for integration into subsequent `gSpreadComp` modules. This file is crucial for the downstream analysis and should be retained.
The **gtdbtk_result.tsv**: This file consolidates the results from GTDBtk, providing comprehensive information on taxonomy assignments in the GTDBtk format.
For a detailed description of the other files, the user can go to the [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/38/23/5315/6758240) page.

### Genome Quality Estimation using CheckM

The quality module in `gSpreadComp` uses CheckM to estimate the quality of genomes. To run this module, use the `mspreadcomp quality` command. Below are the available options for this module:

```console
mspreadcomp quality --help

Usage: mspreadcomp quality [options] --genome_dir genome_folder -o output_dir
Options:
    --genome_dir STR    folder with the genomes to estimate quality (in fasta format)
    --extension STR     fasta file extension (e.g. fa or fasta) [default: fa]
    -o STR              output directory
    -t INT              number of threads [default: 1]
```

### Running the Quality Module

1. Ensure you are in the `test_gspread_run` folder created in the previous step.
2. Place your genomes in the `01_input_genomes` subfolder within the test run folder.
3. Create an output folder within the test run folder, for example, `04_gspread_checkm_quality`.

Assuming you have placed your genomes in `01_input_genomes` and your output folder is `04_gspread_checkm_quality`, your command will look like:

```console
$ mspreadcomp quality --genome_dir ./01_input_genomes/ --extension fa -o ./04_gspread_checkm_quality/ -t 25
```

Run this command, and once it's completed, you can proceed to inspect the output in the `04_gspread_checkm_quality` folder.

### Exploring the Output of the Quality Module

After running the Quality Module, you will find the output in the specified output directory, structured as follows:

```
04_gspread_checkm_quality/
├── bins
├── checkm_df_format_gSpread.csv
├── checkm.log
├── lineage.ms
├── outputcheckm.tsv
└── storage
```

#### Main Output File: checkm_df_format_gSpread.csv. Understanding the Output

The `checkm_df_format_gSpread.csv` file contains quality information for each genome in a format that is ready for integration into subsequent `gSpreadComp` modules. This file is crucial for the downstream analysis and should be retained.
The **outputcheckm.tsv** is the main output file from CheckM itself, consolidating the results and providing comprehensive information on genome quality.

For a detailed description of the other files, the user can go to the [**CheckM**](https://genome.cshlp.org/content/25/7/1043) page.



