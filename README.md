# Switch seq
Analysis of tite-seq data using mean bin approach described in [A.M.Phillips, et al. 2021. eLife.](https://doi.org/10.7554/eLife.71393)

## Input files
In order to run the following analysis, you will need the following input files:

- designed sequences file: csv file containing design ids and the amino acid sequence of the design
- unique sequence files for each sample: csv files with unique sequences in each sample
- flow statistics file: csv file with sorting statistics for the tite-seq sort (can be exported from sorter after sorting)

## Step 1: Get unique sequences matching library designs in each sample

Run fastq_to_counts.py using the designed sequence file and the directory containing the unique sequence files as inputs.

For now, this script parses directly from Azenta unique amino acid sequences. In the future, will update to first pair reads and match sequences from raw fastq files. 

Outputs count_table.tsv containing design id and counts for that sequence for all sort concentrations and bins.

## Step 2: Format flow statistics table

First format flow statistics file by adding a column for library_id, sample_id, and concentration to relevant sorted samples. See example input files for example. 

Run format_fluor_csv.py using the above csv file as input. 

Outputs flow_table.csv containing sample_id, library_id, concentration of binder, cell count, mean_log10_APC, std_log_10_APC, bin, and integer concentration for each sorted concentration and bin. 

## Step 3: Calculate kd values.

Run kd_inference.py with the counts_table.tsv and flow_table.csv generated above as inputs. 

Outputs Kds.csv with design_id, Kd, A (max--expression), B (min--background), errs, covs, and Kds in Molar. 

For interactive plotting and viewing fits, use the associated [jupyter notebook](plotting/plot_curve_fits.ipynb).
