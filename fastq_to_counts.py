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
ngs_input = 'Azenta unique seqs' #'Azenta unique seqs', 'fastq'
designed_sequences_file = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240304_ngs_data/20240304_all_seqs_r1_order.csv'
summary_df = ''
outdir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240507_r1-hs_tite-seq_rep2_ngs'

unique_seqs_dir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240507_r1-hs_tite-seq_rep2_ngs/30-1023075317/UniqueSeq'

####

os.makedirs(outdir, exist_ok=True)

#####
#start from counts table for each sample (Azenta unique seqs)
if ngs_input == 'Azenta unique seqs':

    #add amino acids from display construct to design sequences
    designed_seqs_df = pd.read_csv(designed_sequences_file)
    designed_seqs_df = designed_seqs_df.rename(columns={'Name':'design'})
    
    unique_seqs_files = os.listdir(unique_seqs_dir)

    dfs_to_concat = []

    for file in unique_seqs_files:
        if file.endswith('Unique_AA.csv'):
            sample_name = file.split('_U')[0]
            unique_seqs_df = pd.read_csv(f"{unique_seqs_dir}/{file}")

            design_counts_dict = {}

            for i, seq in enumerate(designed_seqs_df['aa_sequence']):
                counts = unique_seqs_df[unique_seqs_df['Unique Amino Acid'].str.contains(seq, na=False)][' Unique Amino Acid Count'].sum()
                design_counts_dict[designed_seqs_df.at[i, 'design']] = counts

            matching_seqs_df = pd.DataFrame.from_dict(design_counts_dict, orient='index', columns=['count']).reset_index().rename(columns={'index':'design'})    
            matching_seqs_df = matching_seqs_df.fillna(0)
            matching_seqs_df.to_csv(f"{outdir}/{sample_name}_matching_seqs.csv")
            matching_seqs_df['sample_name'] = sample_name
            matching_seqs_df = matching_seqs_df.pivot(index='design', columns='sample_name', values='count')
            dfs_to_concat.append(matching_seqs_df)

    counts_df = pd.concat(dfs_to_concat, axis=1)
    counts_df.to_csv(f"{outdir}/count_table.tsv", sep='\t')

#####

#####
#add option to start from pooled sequencing data and parse by UMI 
#if ngs_input == 'fastq':


#consider reporting stats on sequences matching library

                    
