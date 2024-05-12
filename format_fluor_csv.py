
'''
returns df:
columns: sample_id, library_id, concentration_float, cell count, mean_log10_APC, std_log10_APC, bin, concentration
'''

#import
import pandas as pd
import numpy as np

####USER DEFINED INPUTS####
fluor_csv = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240507_r1-hs_tite-seq_rep2_ngs/20240412_tite_seq_stats.csv'

#list of gate names corresponding to bins 1-4
gate_names = ['bin1', 'bin2', 'bin3', 'bin4']

#list of corresponding bin number
bin_names = [1, 2, 3, 4]

outdir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240507_r1-hs_tite-seq_rep2_ngs'
####

####FUNCTIONS####

fluor_df = pd.read_csv(fluor_csv)

#replace empty bins with no events in gate with 0
fluor_df = fluor_df.replace('####', 0)

#drop rows with no id info--conditions used for setting sort gates
fluor_df = fluor_df.dropna(subset=['sample_id'], axis=0)

#transform pM concentrations to -log10(M)
fluor_df['concentration_float'] = -1*np.log10(fluor_df['concentration'].astype(float)*1e-12)

def get_vals_per_bin(df, colname, bin_num):

    #subset on relevant columns
    df = df[['sample_id', 'library_id', 'concentration_float', f'{colname} Events', f'{colname} Mean APC-A', f'{colname} SD APC-A']]
    
    #log10 transform fluorescence values
    df[[f'{colname} Mean APC-A', f'{colname} SD APC-A']] = df[[f'{colname} Mean APC-A', f'{colname} SD APC-A']].astype(float).apply(np.log10)

    #rename columns
    df = df.rename(columns={f'{colname} Events':'cell count', f'{colname} Mean APC-A':'mean_log10_APC', f'{colname} SD APC-A':'std_log10_APC'})
    
    df['bin'] = bin_num
    df['sample_id'] = df['sample_id'] + '-' + str(bin_num)
    
    #create column for concentration category
    df['concentration'] = df['sample_id'].str.split('-', expand=True)[2].astype(int)

    return df
####

dfs_to_concat = []

for (gate, bin) in zip(gate_names, bin_names):
    df = get_vals_per_bin(fluor_df, gate, bin)
    dfs_to_concat.append(df)

flow_df = pd.concat(dfs_to_concat)
flow_df.to_csv(f'{outdir}/flow_table.csv')
