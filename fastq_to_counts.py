####IMPORT####
import os
import pandas as pd
####

'''Script to parse fastq input data in various forms for tite-seq analysis

Adpated from Hugh Haddox: https://github.com/Haddox/prot_stab_analysis_pipeline

and Angela Phillips: https://github.com/amphilli/CH65-comblib/tree/main/Kd_Inference

Returns a counts table with a row for each design and column for each sample
'''

####USER DEFINED INPUTS####
ngs_input = 'Azenta unique seqs' #'Azenta unique seqs', 'Azenta fastq', 'fastq'
designed_sequences_file = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240304_ngs_data/20240304_all_seqs_r1_order.csv'
summary_df = ''
outdir = '/Users/stephaniecrilly/test'

unique_seqs_dir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240304_ngs_data/30-992575946/UniqueSeq'

fastq_dir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240304_ngs_data/30-992575946/00_fastq'
pear_path = '/Users/stephaniecrilly/pear-src-0.9.11/src/pear' #update with path to python package
five_prime_flanking_seq = 'CATATG'
three_prime_flanking_seq = 'CTCGAG'
####

os.makedirs(outdir, exist_ok=True)

#####
#start from counts table for each sample (Azenta unique seqs)
if ngs_input == 'Azenta unique seqs':

    #add amino acids from display construct to design sequences
    designed_seqs_df = pd.read_csv(designed_sequences_file)
    designed_seqs_df = designed_seqs_df.rename(columns={'Name':'design'})
    designed_seqs_df['sequence'] = 'SASHM' + designed_seqs_df['aa_sequence'].astype(str) + 'LEGGG'
    
    unique_seqs_files = os.listdir(unique_seqs_dir)

    dfs_to_concat = []

    for file in unique_seqs_files:
        if file.endswith('Unique_AA.csv'):
            sample_name = file.split('_U')[0]
            unique_seqs_df = pd.read_csv(f"{unique_seqs_dir}/{file}")
            matching_seqs_df = pd.merge(designed_seqs_df, unique_seqs_df, left_on='sequence', right_on='Unique Amino Acid', how='left')
            matching_seqs_df = matching_seqs_df.drop('Unique Amino Acid', axis=1)
            matching_seqs_df = matching_seqs_df.rename(columns={' Unique Amino Acid Count':'count'})
            matching_seqs_df = matching_seqs_df.fillna(0)
            matching_seqs_df.to_csv(f"{outdir}/{sample_name}_matching_seqs.csv")
            matching_seqs_df['sample_name'] = sample_name
            matching_seqs_df = matching_seqs_df.pivot(index='design', columns='sample_name', values='count')
            dfs_to_concat.append(matching_seqs_df)

    counts_df = pd.concat(dfs_to_concat, axis=1)
    counts_df.to_csv(f"{outdir}/count_table.tsv", sep='\t')

#####
#start from fastq samples already separated into diff samples (Azenta fastq)
#adapted from Hugh Haddox protease analysis scripts
#alternatively look into flashpy: https://github.com/ponnhide/flashpy
if ngs_input == 'Azenta fastq':
     #pair forward and reverse reads
    paired_FASTQ_files_dir = os.path.join(outdir, 'paired_FASTQ_files')
    os.makedirs(paired_FASTQ_files_dir)

    files = os.listdir(fastq_dir)

    r1_files = []
    for file in files:
        if '_R1_' in file:
            r1_files.append(file)

    for r1_file in r1_files:
        paired_filename = r1_files.split('_R1')[0]
        r2_file = r1_file.replace('R1_', 'R2_')

        outfile_prefix = os.path.join(paired_FASTQ_files_dir, f"{paired_filename}")
        logfile = '{0}.log'.format(outfile_prefix)

        cmd = [
            pear_path,
            '-f', r1_file,
            '-r', r2_file,
            '-o', outfile_prefix,
            '-j', '20'
        ]

        log = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        out, err = log.communicate()
        with open(logfile, 'wb') as f:
            f.write(out)

#parse fastq files for matching sequences and assemble counts table

#####
#add option to start from pooled sequencing data and parse by UMI 
#if ngs_input == 'fastq':


#consider reporting stats on sequences matching library

                    
