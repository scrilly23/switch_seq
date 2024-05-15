import os
import pandas as pd
import numpy as np

####USER DEFINED INPUTS####
kds_rep1_path = '/Users/stephaniecrilly/Kortemme_lab/switch_seq/Kds.csv'
kds_rep2_path = '/Users/stephaniecrilly/Kortemme_lab/switch_seq/r1-hs_tite-seq_rep2_fixed_swapped_samples_Kds.csv'
prot_stab_path = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240304_ngs_data/30-993725825/rocklin_analysis/20240304_r1hs_rep1_ Trypsin_protein_w_design_type.csv'
exp_path = '/Users/stephaniecrilly/Kortemme_lab/switch_seq/r1-hs-pilot_rep1_expression.csv'

outdir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240514_r1-hs_tite_seq_reps_filtered'
exp_id = 'r1-hs_tite_seq_rep1_rep2_filtered_bbs'

#cutoffs for filtering designs
upper_log_kd = 0.5
lower_log_kd = -0.5

upper_ss = 1.5
lower_ss = 1.0

upper_exp = np.inf
lower_exp = 0

A_min = 1e4

r2_cutoff = 0.7

####

if not os.path.exists(outdir):
   os.makedirs(outdir)

#read in data
kds_rep1_df = pd.read_csv(kds_rep1_path).rename({'geno':'design'}, axis=1)
kds_rep2_df = pd.read_csv(kds_rep2_path).rename({'geno':'design'}, axis=1)
prot_stab_df = pd.read_csv(prot_stab_path).rename({'Name':'design'}, axis=1)
exp_df = pd.read_csv(exp_path).rename({'Name':'design'}, axis=1)

#normalize Kds to expected kds of controls for each experiment
kds_rep1_df['Kds_norm_to_expected'] = ((kds_rep1_df['Kds'] - kds_rep1_df.loc[kds_rep1_df['design'] == 'bm01_ALFA_t3', 'Kds'].item())
                                    .where(kds_rep1_df['design'].str.contains('min_0'), other=(kds_rep1_df['Kds'] - kds_rep1_df.loc[kds_rep1_df['design'] == 'bm01_ALFA_t6', 'Kds'].item())))

kds_rep2_df['Kds_norm_to_expected'] = ((kds_rep2_df['Kds'] - kds_rep2_df.loc[kds_rep2_df['design'] == 'bm01_ALFA_t3', 'Kds'].item())
                                    .where(kds_rep2_df['design'].str.contains('min_0'), other=(kds_rep2_df['Kds'] - kds_rep2_df.loc[kds_rep2_df['design'] == 'bm01_ALFA_t6', 'Kds'].item())))

#bin normalized Kds
bins = [-10, -1, lower_log_kd, upper_log_kd, 10]
labels = ['lower binding', 'intermediate lower binding', 'expected binding', 'higher binding']

kds_rep1_df['Kds_norm_to_expected_binned'] = pd.cut(kds_rep1_df['Kds_norm_to_expected'], bins=bins, labels=labels)
kds_rep2_df['Kds_norm_to_expected_binned'] = pd.cut(kds_rep2_df['Kds_norm_to_expected'], bins=bins, labels=labels)

#create master df
merge_kds_df = kds_rep1_df.merge(kds_rep2_df, on='design', suffixes=('_rep1', '_rep2'))

prot_stab_kds_df = pd.merge(merge_kds_df, prot_stab_df, on='design', how='left')
all_merge_df = pd.merge(prot_stab_kds_df, exp_df, on='design', how='left')

all_merge_df.to_csv(f"{outdir}/{exp_id}_all_merged.csv")

#filter designs
filtered_df = all_merge_df[(all_merge_df['Kds_norm_to_expected_binned_rep1'] == 'expected binding') &
                           (all_merge_df['Kds_norm_to_expected_binned_rep2'] == 'expected binding') &
                           (all_merge_df['stabilityscore'] < upper_ss) &
                           (all_merge_df['stabilityscore'] > lower_ss) &
                           (all_merge_df['As_rep1'] > A_min) &
                           (all_merge_df['As_rep2'] > A_min) &
                            (all_merge_df['errs_rep1'] > r2_cutoff) &
                           (all_merge_df['errs_rep2'] > r2_cutoff)]

filtered_df.to_csv(f"{outdir}/{exp_id}_filtered_designs.csv")