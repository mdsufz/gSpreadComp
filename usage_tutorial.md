# gSpreadComp Usage Tutorial

Welcome to the `gSpreadComp` usage tutorial! This guide is designed to help you understand how to use the tool effectively through practical examples and step-by-step instructions.

You can use genomes of your own, but if you need some genomes for testing, you can use the one [here](https://e.pcloud.link/publink/show?code=kZoUHsZq72qCkzuJmunUnT0JIdVTkqghAGV).

## Introduction

`gSpreadComp` is designed to work with genomes in fasta format and requires a metadata table in CSV format containing a target variable. For example tables, please refer to [example_metadata_table_link](test_data/02_metadata_gspread_sample.csv) and [example_genome_fasta_link](https://e.pcloud.link/publink/show?code=kZoUHsZq72qCkzuJmunUnT0JIdVTkqghAGV). To recover prokaryotic genomes from metagenomic samples, you can use tools like [MuDoGeR](https://github.com/mdsufz/MuDoGeR).

### Prerequisites

Before proceeding with the tutorial, make sure you have:
1. Genomes in fasta format.
2. A metadata table in CSV format with a target variable. Your metadata table must be formatted correctly, i.e., including a column named "Library" to identify the source sample, a column named "Genome" to identify the genome, and a column named "Target" to identify your target variable. Refer to example [here](test_data/02_metadata_gspread_sample.csv)

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

The taxonomy module in `gSpreadComp` uses [GTDBtk](https://academic.oup.com/bioinformatics/article/38/23/5315/6758240) for taxonomy assignment. To run this module, use the `gspreadcomp taxonomy` command. Below are the available options for this module:

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

The user can find an example of the expected gtdb_df_format_gSpread.csv [here](test_data/gtdb_df_format_gSpread.csv).

The `gtdb_df_format_gSpread.csv` file contains taxonomy information for each genome in a format that is ready for integration into subsequent `gSpreadComp` modules. This file is crucial for the downstream analysis and should be retained.
The **gtdbtk_result.tsv**: This file consolidates the results from GTDBtk, providing comprehensive information on taxonomy assignments in the GTDBtk format.
For a detailed description of the other files, the user can go to the [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/38/23/5315/6758240) page.

### Genome Quality Estimation using CheckM

The quality module in `gSpreadComp` uses [CheckM](https://genome.cshlp.org/content/25/7/1043) to estimate the quality of prokaryotic genomes. To run this module, use the `gspreadcomp quality` command. Below are the available options for this module:

```console
gspreadcomp quality --help

Usage: gspreadcomp quality [options] --genome_dir genome_folder -o output_dir
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

Assuming you have placed your genomes in `01_input_genomes` and your output folder is `04_gspread_checkm_quality`, your command will look like this:

```console
$ gspreadcomp quality --genome_dir ./01_input_genomes/ --extension fa -o ./04_gspread_checkm_quality/ -t 25
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

The user can find an example of the expected checkm_df_format_gSpread.csv [here](test_data/checkm_df_format_gSpread.csv).

The `checkm_df_format_gSpread.csv` file contains quality information for each genome in a format that is ready for integration into subsequent `gSpreadComp` modules. This file is crucial for the downstream analysis and should be retained.
The **outputcheckm.tsv** is the main output file from CheckM itself, consolidating the results and providing comprehensive information on genome quality.

For a detailed description of the other files, the user can go to the [**CheckM**](https://genome.cshlp.org/content/25/7/1043) page.

### ARGs Annotation using DeepARG

The ARGs module in `gSpreadComp` uses [DeepARG](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0401-z) to predict the Antimicrobial Resistance Genes (ARGs) in a genome. To run this module, use the `gspreadcomp args` command. Below are the available options for this module:

```console
gspreadcomp args --help

Usage: gspreadcomp args [options] --genome_dir genome_folder -o output_dir
Options:
    --genome_dir STR    folder with the genomes to be classified (in fasta format)
    --extension STR     fasta file extension (e.g. fa or fasta) [default: fa]
    --min_prob NUM      Minimum probability cutoff for DeepARG [Default: 0.8]
    --arg_alignment_identity NUM   Identity cutoff for sequence alignment for DeepARG [Default: 35]
    --arg_alignment_evalue NUM     Evalue cutoff for DeepARG [Default: 1e-10]
    --arg_alignment_overlap NUM    Alignment read overlap for DeepARG [Default: 0.8]
    --arg_num_alignments_per_entry NUM   Diamond, minimum number of alignments per entry [Default: 1000]
    -o STR              output directory
```

### Running the ARGs Module

1. Ensure you are in the `test_gspread_run` folder created in the previous steps.
2. Place your genomes in the `01_input_genomes` subfolder within the test run folder.
3. Create an output folder within the test run folder, for example, `05_gspread_deeparg_args`.

Assuming you have placed your genomes in `01_input_genomes` and your output folder is `05_gspread_deeparg_args`, your command will look like this:

```console
$ gspreadcomp args --genome_dir ./01_input_genomes/ --extension fa -o ./05_gspread_deeparg_args/
```

Run this command, and once it's completed, you can proceed to inspect the output in the `05_gspread_deeparg_args` folder.

### Inspecting the Output of the ARGs Module

After running the ARGs Module, you will find the output in the specified output directory, structured as follows:

```
05_gspread_deeparg_args/
├── deeparg_df_format_gSpread.csv
├── deeparg_df_combined_raw.csv
├── genome_name_1.fa
│   ├── genome_name_1.fa_deeparg_out.align.daa
│   ├── genome_name_1.fa_deeparg_out.align.daa.tsv
│   ├── genome_name_1.fa_deeparg_out.mapping.ARG
│   └── genome_name_1.fa_deeparg_out.mapping.potential.ARG
├── genome_name_2.fa
│   ├── genome_name_2.fa_deeparg_out.align.daa
│   ├── genome_name_2.fa_deeparg_out.align.daa.tsv
│   ├── genome_name_2.fa_deeparg_out.mapping.ARG
│   └── genome_name_2.fa_deeparg_out.mapping.potential.ARG
└── genomes_with_no_found_deeparg.csv
```

#### Understanding the Output Files and Directories

The user can find an example of the expected deeparg_df_format_gSpread.csv [here](test_data/deeparg_df_format_gSpread.csv).

1. **deeparg_df_format_gSpread.csv**: This is the format-ready main output file, containing formatted ARGs annotation information per genome. It's ready for integration into subsequent `gSpreadComp` modules.
2. **deeparg_df_combined_raw.csv**: This file combines the raw output from DeepARG for all genomes analyzed.
3. **genome_name.fa Directories**: For each genome analyzed, a separate directory is created, named after the genome. To get a detailed description of its content, the user can read the [DeepARG](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0401-z) documentation
4. **genomes_with_no_found_deeparg.csv**: This file lists the genomes for which no ARGs were found by DeepARG.

### Plasmid Identification using PlasFlow

The Plasmid module in `gSpreadComp` uses [PlasFlow](https://academic.oup.com/nar/article/46/6/e35/4807335) to predict if a sequence within a fasta file is a chromosome, plasmid, or undetermined. To run this module, use the `gspreadcomp plasmid` command. Below are the available options for this module:

```console
gspreadcomp plasmid --help

Usage: gspreadcomp plasmid [options] --genome_dir genome_folder -o output_dir
Options:
    --genome_dir STR    folder with the genomes to be classified (in fasta format)
    --extension STR     fasta file extension (e.g. fa or fasta) [default: fa]
    --threshold NUM     threshold for probability filtering [default: 0.7]
    -o STR              output directory
```

### Running the Plasmid Identification Module

1. Ensure you are in the `test_gspread_run` folder created in the previous steps.
2. Place your genomes in the `01_input_genomes` subfolder within the test run folder.
3. Create an output folder within the test run folder, for example, `06_gspread_plasmids`.

Assuming you have placed your genomes in `01_input_genomes` and your output folder is `06_gspread_plasmids`, your command will look like this:

```console
$ gspreadcomp plasmid --genome_dir ./01_input_genomes/ --extension fa -o ./06_gspread_plasmids/
```

Run this command, and once it's completed, you can proceed to inspect the output in the `06_gspread_plasmids` folder. Below is the expected output structure and explanation of each output file.

### Inspecting the Output of the Plasmid Module

After running the Plasmid Module, you will find the output in the specified output directory, structured as follows:

```
06_gspread_plasmids/
├── genome_name_1.fa
│   ├── genome_name_1.fa_plasflow_out.tsv
│   ├── genome_name_1.fa_plasflow_out.tsv_chromosomes.fasta
│   ├── genome_name_1.fa_plasflow_out.tsv_plasmids.fasta
│   └── genome_name_1.fa_plasflow_out.tsv_unclassified.fasta
├── genome_name_2.fa
│   ├── genome_name_2.fa_plasflow_out.tsv
│   ├── genome_name_2.fa_plasflow_out.tsv_chromosomes.fasta
│   ├── genome_name_2.fa_plasflow_out.tsv_plasmids.fasta
│   └── genome_name_2.fa_plasflow_out.tsv_unclassified.fasta
├── genomes_with_no_found_plasflow.csv
└── plasflow_combined_format_gSpread.csv
```

#### Understanding the Output Files and Directories

The user can find an example of the expected plasflow_combined_format_gSpread.csv [here](test_data/plasflow_combined_format_gSpread.csv).

1. **genome_name.fa Directories**: For each genome analyzed, a separate directory is created, named after the genome. It contains the following files:
   - **genome_name.fa_plasflow_out.tsv**: This file contains the PlasFlow results in tab-separated values format.
   - **genome_name.fa_plasflow_out.tsv_chromosomes.fasta**: This file contains sequences predicted to be chromosomes.
   - **genome_name.fa_plasflow_out.tsv_plasmids.fasta**: This file contains sequences predicted to be plasmids.
   - **genome_name.fa_plasflow_out.tsv_unclassified.fasta**: This file contains sequences that could not be classified as either plasmids or chromosomes.
   For a detailed description of the files, the user can read the [Plasflow](https://academic.oup.com/nar/article/46/6/e35/4807335?login=true) documentation
2. **genomes_with_no_found_plasflow.csv**: This file lists the genomes for which no sequences were found by PlasFlow. Hopefully, it will be empty.

3. **plasflow_combined_format_gSpread.csv**: This is the format-ready main output file containing formatted PlasFlow results per genome. It's ready for integration into subsequent `gSpreadComp` modules.

### Pathogens Annotation using Virulence Factors Databases

The Pathogens module in `gSpreadComp` aligns the provided genomes against selected Virulence Factors databases and formats the output. 
The pathogens module essentially uses BLAST to align your genomes with defined Virulence Factors databases. Here, the user can find the [Victors](https://academic.oup.com/nar/article/47/D1/D693/5144967?) database and the [VFDB](https://academic.oup.com/nar/article/50/D1/D912/6446532) database.

To run this module, use the `gspreadcomp pathogens` command. Below are the available options for this module:

```console
gspreadcomp pathogens --help

Usage: gspreadcomp pathogens [options] --genome_dir genome_folder -o output_dir
Options:
    --genome_dir STR    folder with the genomes to be aligned against Virulence factors (in fasta format)
    --extension STR     fasta file extension (e.g. fa or fasta) [default: fa]
    --evalue NUM        evalue, expect value, threshold as defined by NCBI-BLAST [default: 1e-50]
    --vf STR            select the virulence factors database to be used (e.g. victors, vfdb or both) [default: both]
    -t INT              number of threads
    -o STR              output directory
```

### Running the Pathogens Module

1. Ensure you are in the `test_gspread_run` folder created in the previous steps.
2. Place your genomes in the `01_input_genomes` subfolder within the test run folder.
3. Create an output folder within the test run folder, for example, `07_gspread_pathogens`.

Assuming you have placed your genomes in `01_input_genomes` and your output folder is `07_gspread_pathogens`, your command will look like this:

```console
$ gspreadcomp pathogens --genome_dir ./01_input_genomes/ --extension fa -o ./07_gspread_pathogens/ --vf both -t 25
```

Run this command, and once it's completed, you can proceed to inspect the output in the `07_gspread_pathogens` folder.

### Inspecting the Output of the Pathogens Module

After running the Pathogens Module, you will find the output in the specified output directory. Below is the expected output structure and explanation of each output file.

```console
07_gspread_pathogens/
├── genome_name_1.fa
│   ├── vfdb_genome_name_1.out
│   └── victors_genome_name_1.out
├── genome_name_2.fa
│   ├── vfdb_genome_name_2.out
│   └── victors_genome_name_2.out
├── vfdb_format_gSpread.csv
├── vfdb_headers.txt
├── vfdb_merged.out
├── vfdb_per_genome_unique_count.csv
├── victors_format_gSpread.csv
├── victors_db_headers.txt
├── victors_merged.out
└── victors_per_genome_unique_count.csv
```

#### Understanding the Output Files and Directories

The user can find an example of the expected victors_format_gSpread.csv [here](test_data/victors_format_gSpread.csv). An equivalent output is generated if the user uses the VFDB instead of Victors database.

1. **vfdb_format_gSpread.csv & victors_format_gSpread.csv**: These are the format-ready main output files containing formatted virulence factors results per genome. They're ready for integration into subsequent `gSpreadComp` modules.
2. **vfdb_headers.txt & victors_db_headers.txt**: These files contain the headers for the VFDB and Victors databases respectively.
3. **vfdb_merged.out & victors_merged.out**: These files contain the merged results of the alignments against the VFDB and Victors databases, respectively.
4. **vfdb_per_genome_unique_count.csv & victors_per_genome_unique_count.csv**: These files contain the count of unique virulence factors per genome for the VFDB and Victors databases respectively.

#### Main Output Files: vfdb_format_gSpread.csv.csv & victors_format_gSpread.csv.csv

The `vfdb_format_gSpread.csv.csv` and `victors_format_gSpread.csv.csv` files contain virulence factors results for each genome in a format that is ready for integration into subsequent `gSpreadComp` modules. These files are crucial for downstream analysis and should be retained.

## Using custom files not generated with gSpreadComp

If the user wants to generate custom, quality, taxonomy, gene annotation, plasmid identification, or Virulence Factors annotation files outside gSpreadComp, it is important to maintain the table formatting. 

All the examples of input tables used by the gSpread module are [here](test_data).

 **Most important is to keep the column naming as seen in the example tables.**

## gSpread Module: Main Analysis and Downstream Processing

The `gspread` module is the final step in the `gSpreadComp` pipeline, integrating the previous modules' outputs to comprehensively analyze gene spread, plasmid-mediated horizontal gene transfer, and pathogenic risk.

### Usage:

To view the available options for the `gspread` module, use the following command:

```bash
gspreadcomp gspread --help
```

This will display the available parameters and their descriptions.

### Running the Module:

If you've followed our tutorial steps sequentially, you should have the required inputs ready for the `gspread` module. Here's how to execute the module with the processed outputs:

```bash
gspreadcomp gspread --gtdbtk ./03_gspread_gtdb_taxonomy/gtdb_df_format_gSpread.csv  --checkm ./04_gspread_checkm_quality/checkm_df_format_gSpread.csv --gene ./05_gspread_deeparg_args/deeparg_df_combined_gSpreadformat.csv --meta ./02_metadata_gspread_sample.csv --plasmid ./06_gspread_plasmids/plasflow_output_combined.csv --vf ./07_gspread_pathogens/victors_ann_format.csv -t 25 -o ./08_gspread_results/ --target_gene_col Gene_id
```

### Inspecting the Output of the gSpread Module

After running the command, the `gspread` module will produce several output files in the specified output directory (`./08_gspread_results/` in our example):

```
08_gspread_results/
├── common_tax_target.csv
├── gSpread_risk_report.html
├── gene_pairwise_comp_results
├── gene_spread_results
├── genome_quality_norm
├── hgt_events_results
├── mags_complete_annotation.csv
├── mags_summary_results.csv
├── network_vis_files
└── pathogens_results
```

Among these, the `gSpread_risk_report.html` provides a comprehensive report on the risk assessment, while other files and directories contain detailed results from various analyses performed by the module.

Several of the generated files are used in the report. However, a detailed inspection of all the outputs may also help the user.

The main files are the mags_complete_annotation.csv and the mags_summary_results.csv, which contain a result compilation and a risk metric for every given genome.

MORE DETAILS ARE COMING!
