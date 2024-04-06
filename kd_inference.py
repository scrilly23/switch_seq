'''Modified from Angela Phillips: https://github.com/amphilli/CH65-comblib/blob/main/Kd_Inference/scripts/kd_inference.py

Fits Kds to tite-seq counts data

Inputs: counts_table, flow_table
 '''
# importing things
import pandas as pd
import numpy as np
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
mpl.rcParams['figure.dpi'] = 300

####USER DEFINED INPUTS####
counts_table_path = '/Users/stephaniecrilly/test/count_table.csv'
flow_table_path = '/Users/stephaniecrilly/test/flow_table.csv'
bounds_dict = {}
construct = 'r1-hs'

####

counts_table = pd.read_csv(counts_table_path, sep="\t")
flow_table = pd.read_csv(flow_table_path)


# # reading in counts
# df = pd.read_csv(counts_table, sep="\t")
# df["geno"] = df.geno.apply(lambda x: f"{int(x):0{conf['num_mutations']}d}")
# # reading in fluorescence data
# fluor = pd.read_csv(snakemake.input.fluorescence, sep="\t")
# fluor = fluor[fluor.replicate == replicate]


# plotting with error
# define function

#why is c negative????

def extractKd(concentrations, bins_mean, bins_std):
    """
        Arguments: concentrations and exponential mean bin numbers (bin number: 1, 2, 3, 4)
        Return: -log10(Kd), and the r^2
    """
    popt, pcov = scipy.optimize.curve_fit(sigmoid, concentrations,
                                          bins_mean,
                                          p0=[(-9), 10**(4), 10**(2)], #initial guess for parameters
                                          sigma=bins_std, absolute_sigma=True,
                                          bounds=[(bounds_dict['bounds_log10Kd_min'],
                                                   bounds_dict['bounds_A_min'],
                                                   bounds_dict['bounds_B_min']),
                                                  (bounds_dict['bounds_log10Kd_max'],
                                                   bounds_dict['bounds_A_max'],
                                                   bounds_dict['bounds_B_max'])],
                                            nan_policy='omit' #do calculations without NaNs
                                          ) #got rid if maxfev=400000 for now, number of iterations for fitting
    return(-1*popt[0], popt[1], popt[2], 1 - np.sum((sigmoid(concentrations, *popt) - bins_mean)**2)/np.sum((bins_mean - bins_mean.mean())**2), np.sqrt(np.diag(pcov))[0])


def sigmoid(c, Kd, A, B):
    return np.log10(A * (10**c/((10**c)+(10**Kd))) + B)


def compute_Kds(flow_table, counts_table):
    #sample_info = sample_info[~sample_info.concentration_float.isna()].copy()
    #get parameters needed from input tables to compute Kd
    nb_bins = flow_table.bin.nunique()
    nb_concs = flow_table.concentration.nunique()
    concentrations = flow_table.concentration.unique()
    nb_genos = counts_table.shape[0]

    #initialize empty arrays
    probas = np.zeros((nb_bins, nb_genos, nb_concs))
    counts = np.zeros((nb_bins, nb_genos, nb_concs))
    cells = np.zeros((nb_bins, nb_concs))
    meanfluor, stdfluor = np.zeros((2, nb_bins, nb_concs))
    
    #my attempt to create the above arrays from dfs directly 
    #make counts array
#     id_conc_bin = flow_table[['sample_id', 'concentration', 'bin']]
#     id_conc_bin = id_conc_bin.transpose()
#     id_conc_bin.columns = id_conc_bin.iloc[0]
#     id_conc_bin = id_conc_bin.drop(id_conc_bin.index[0])

#     print(id_conc_bin)

#     dfs_to_concat = [counts_table, id_conc_bin]

#     counts_id_conc_bin = pd.concat(dfs_to_concat).to_numpy()

#     return counts_id_conc_bin

# test = compute_Kds(flow_table, counts_table)
# print(test.shape)
    #convert df to array



#     #populate empty arrays with relevant data from input tables
    for bb, gate in enumerate(range(1, nb_bins+1)):
        for cc, conc in enumerate(concentrations):
            counts[bb, :, cc] = counts_table[f"{construct}-{str(conc)}-{gate}"]
            cells[bb, cc] = sample_info[(sample_info.concentration == conc) #make dict of concentration number with actual concentration #add column to flow table
                                        & (sample_info.bin == gate)]["cell count"].iloc[0] #in sample_onfo tsv conc is an object and bin is an int
            meanfluor[bb, cc] = fluor[(fluor.concentration == conc)
                                      & (fluor.bin == gate)]["mean_log10_PEs"].iloc[0]
            stdfluor[bb, cc] = fluor[(fluor.concentration == conc)
                                     & (fluor.bin == gate)]["std_log10_PEs"].iloc[0]

#     probas = counts / (counts.sum(axis=1)[:, None, :]) * cells[:, None, :] #what is probas? #seems to be norm of counts to cells as well as total number of counts
#     probas = probas / probas.sum(axis=0)[None, :, :]
#     mean_log10_fluor = (probas * meanfluor[:, None, :]).sum(axis=0)
#     std_log10_fluor = np.sqrt((stdfluor[:, None, :]**2 * probas**2
#                                + meanfluor[:, None, :]**2 * probas**2 / (1e-22 + counts)).sum(axis=0))


#     # fit of Kd
#     Kds, A, B, err, cov = np.zeros((5, nb_genos))
#     for s in range(nb_genos):
#         notnanindex = [ii for ii in range(nb_concs)
#                        if not np.isnan(mean_log10_fluor[s, ii] + std_log10_fluor[s, ii])]
#         if len(notnanindex) < 4: #i think this is filtering out an genos which have <4 data points for fluor
#             Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
#         # not enough reads
#         if np.sum(counts.sum(axis=0) > snakemake.config['min_number_counts']) < 4:
#             Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
#         else:
#             Kds[s], A[s], B[s], err[s], cov[s] = extractKd(concentrations[notnanindex],
#                                                            mean_log10_fluor[s,
#                                                                             notnanindex],
#                                                            std_log10_fluor[s, notnanindex])

#     return Kds, A, B, err, cov, mean_log10_fluor, std_log10_fluor, concentrations


# #plot a subset of the data
# #need to rework
# df_subset = df.sample(n=8)
# xs = np.linspace(-14, -6, 100)
# Kds, As, Bs, errs, covs, mlog10, slog10, concentrations = compute_Kds(
#     sample_info, fluor, df_subset)

# xlim = (-max(sample_info.concentration_float)-0.5,
#         -min(sample_info[sample_info.concentration_float != 0].concentration_float)+0.5)
# fig, ax = plt.subplots(4, 2, figsize=(10, 10), sharex=True)
# for ii, ax in enumerate(ax.flatten()):
#     ax.errorbar(x=concentrations,
#                 y=mlog10[ii, :], yerr=slog10[ii, :], label=df_subset.geno.iloc[ii])
#     ax.plot(xs, sigmoid(xs, -Kds[ii], As[ii], Bs[ii]))
#     ax.set_xlim(xlim)
#     ax.set_xlabel("$\log_{10}(\mathrm{Concentration})$")
#     ax.set_ylabel("Est. Mean Fluorescence")
# plt.savefig(snakemake.output.plot_test_curve)


# # Compute Kds for the full dataset
# #set up option to run this if above looks good
# Kds, As, Bs, errs, covs, mean_log10_PE, std_log10_PE, concs = compute_Kds(
#     sample_info, fluor, df)
# df["log10Kd"] = Kds
# df["A"] = As
# df["B"] = Bs
# df["r2"] = errs
# df["sigma"] = covs

# # save the Kds values
# df[["geno", "log10Kd", "A", "B", "r2", "sigma"] +
#    [f"mean_log10PE{cc}" for cc in range(0, mean_log10_PE.shape[1])] +
#    [f"std_log10PE{cc}" for cc in range(0, mean_log10_PE.shape[1])] +
#    (["Mean fluorescence expression", "Std fluorescence expression"]
#     if ('F' in sample_info.concentration.unique()) else [])].to_csv(
#     snakemake.output.tsv,
#     sep="\t", index=False)


# with pd.option_context('mode.use_inf_as_na', True): #could be useful for negative inf where no counts
#     fig, ax = plt.subplots()
#     if conf['do_plots']:
#         # Kd distribution
#         sns.histplot(x="log10Kd", data=df.dropna(subset=["log10Kd"]),
#                      ax=ax, label="All intermediates")
#         # No mutations
#         if df.geno.apply(lambda x: '1' not in x).sum() > 0:
#             ax.axvline(x=df[df.geno.apply(lambda x: '1' not in x)
#                             ].log10Kd.iloc[0], label="Original", c="g")
#         # All mutations
#         if df.geno.apply(lambda x: '0' not in x).sum() > 0:
#             ax.axvline(x=df[df.geno.apply(lambda x: '0' not in x)
#                             ].log10Kd.iloc[0], label="Variant", c="r")
#         ax.legend()
#     plt.savefig(snakemake.output.plot_Kd_distribution)

#     # Mean expression distribution
#     fig, ax = plt.subplots()
#     if "Mean fluorescence expression" in df and conf['do_plots']:
#         sns.histplot(x="Mean fluorescence expression",
#                      data=df.dropna(subset=["Mean fluorescence expression"]),
#                      label="All intermediates", ax=ax)
#         if df.geno.apply(lambda x: '1' not in x).sum() > 0:
#             ax.axvline(x=df[df.geno.apply(lambda x: '1' not in x)]["Mean fluorescence expression"].iloc[0],
#                        label="Original", c="g")
#         if df.geno.apply(lambda x: '0' not in x).sum() > 0:
#             ax.axvline(x=df[df.geno.apply(lambda x: '0' not in x)]["Mean fluorescence expression"].iloc[0],
#                        label="Variant", c="r")
#         ax.legend()
#     plt.savefig(snakemake.output.plot_fluo_distribution)

#     # Correlation btw expression and Kd
#     if "Mean fluorescence expression" in df and conf['do_plots']:
#         sns.jointplot(x="Mean fluorescence expression",
#                       y="log10Kd",
#                       data=df.dropna(
#                           subset=["Mean fluorescence expression", "log10Kd"]),
#                       kind="hex", bins='log')
#     else:
#          plt.subplots()
#     plt.savefig(snakemake.output.plot_corr_fluo_Kd)

#     # Correlation btw Kd error and Kd
#     if conf['do_plots']:
#         df["$\log_{10}(\mathrm{error})$"] = np.log(1 - df["r2"])
#         sns.jointplot(x="$\log_{10}(\mathrm{error})$",
#                       y="log10Kd",
#                       data=df.dropna(
#                           subset=["$\log_{10}(\mathrm{error})$", "log10Kd"]),
#                       kind="hex", bins='log')
#     else:
#         plt.subplots()
#     plt.savefig(snakemake.output.plot_corr_Kd_error)


