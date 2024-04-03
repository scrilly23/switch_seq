'''ad hoc script for cleaning up csv file for creating fluorescence table

standardize input file from sort export for future runs

returns df:
columns: sample_id, construct (library), concentration, events, mean_log10_APC, std_log10_APC, bin
'''

#import
import pandas as pd
import numpy as np

####USER DEFINED INPUTS####
fluor_csv = '/Users/stephaniecrilly/Library/CloudStorage/Box-Box/kortemmelab/home/scrilly/helix_sliding/20240129_r1_hs_tite-seq_test/r1_hs_tite_seq_fluor_sample_ids.csv'

#list of gate names corresponding to bins 1-4
gate_names = ['P5', 'P6', 'P4', 'P3']

#list of corresponding bin number
bin_names = [1, 2, 3, 4]

outdir = '/Users/stephaniecrilly/test'
####

####FUNCTIONS####

fluor_df = pd.read_csv(fluor_csv)

#replace empty bins with no events in gate with 0
fluor_df = fluor_df.replace('####', 0)
print(fluor_df.dtypes)

def get_vals_per_bin(df, colname, bin_num):

    #subset on relevant columns
    df = df[['sample_id', 'construct', 'concentration', f'{colname} Events', f'{colname} Mean APC-A', f'{colname} SD APC-A']]
    
    #log10 transform fluorescence values
    df[[f'{colname} Mean APC-A', f'{colname} SD APC-A']] = df[[f'{colname} Mean APC-A', f'{colname} SD APC-A']].astype(float).apply(np.log10)

    #rename columns
    df = df.rename(columns={f'{colname} Events':'events', f'{colname} Mean APC-A':'mean_log10_APC', f'{colname} SD APC-A':'std_log10_APC'})
    
    df['bin'] = bin_num
    df['sample_id'] = df['sample_id'] + '-' + str(bin_num)

    return df
####

dfs_to_concat = []

for (gate, bin) in zip(gate_names, bin_names):
    df = get_vals_per_bin(fluor_df, gate, bin)
    dfs_to_concat.append(df)

flow_df = pd.concat(dfs_to_concat)
flow_df.to_csv(f'{outdir}/flow_table.csv')
