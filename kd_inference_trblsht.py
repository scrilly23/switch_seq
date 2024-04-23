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
import sys
import seaborn as sns

np.set_printoptions(threshold=sys.maxsize)
mpl.rcParams['figure.dpi'] = 300

####USER DEFINED INPUTS####
#counts_table_path = '/Users/stephaniecrilly/test/count_table_ctrl.tsv'
counts_table_path = '/Users/stephaniecrilly/test/count_table.csv'
flow_table_path = '/Users/stephaniecrilly/test/flow_table.csv'
bounds_dict = {'bounds_log10Kd_min':-15, 'bounds_A_min':100, 'bounds_B_min':1, 'bounds_log10Kd_max':-3, 'bounds_A_max':1000000, 'bounds_B_max':100000}
construct = 'r1-hs'
min_num_counts = 20
outdir = '/Users/stephaniecrilly/test'

####

#counts_table = pd.read_csv(counts_table_path, sep="\t")
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

#def extractKd(concentrations, bins_mean, bins_std):
def extractKd(concentrations, bins_mean):
    """
        Arguments: concentrations and exponential mean bin numbers (bin number: 1, 2, 3, 4)
        Return: -log10(Kd), and the r^2
    """
    popt, pcov = scipy.optimize.curve_fit(sigmoid, concentrations,
                                          bins_mean,
                                          #p0=[(-12), 10**(4), 10**(2)], #initial guess for parameters
                                          #sigma=bins_std, absolute_sigma=True,
                                          bounds=[(bounds_dict['bounds_log10Kd_min'],
                                                   bounds_dict['bounds_A_min'],
                                                   bounds_dict['bounds_B_min']),
                                                  (bounds_dict['bounds_log10Kd_max'],
                                                   bounds_dict['bounds_A_max'],
                                                   bounds_dict['bounds_B_max'])],
                                                   maxfev=400000
                                           # nan_policy='omit' #do calculations without NaNs
                                          ) #got rid if maxfev=400000 for now, number of iterations for fitting
    return(-1*popt[0], popt[1], popt[2], 1 - np.sum((sigmoid(concentrations, *popt) - bins_mean)**2)/np.sum((bins_mean - bins_mean.mean())**2), np.sqrt(np.diag(pcov))[0])


def sigmoid(c, Kd, A, B):
    return np.log10(A * (10**c/((10**c)+(10**Kd))) + B)


def compute_Kds(flow_table, counts_table):

    #get parameters needed from input tables to compute Kd
    nb_bins = flow_table.bin.nunique()
    nb_concs = flow_table.concentration.nunique()
    flow_table.concentration_float = -1*np.log10(flow_table.concentration_float*1e12)
    concentrations = flow_table.concentration_float.unique()
    nb_genos = counts_table.shape[0]

    #initialize empty arrays
    probas = np.zeros((nb_bins, nb_genos, nb_concs))
    counts = np.zeros((nb_bins, nb_genos, nb_concs))
    cells = np.zeros((nb_bins, nb_concs))
    meanfluor, stdfluor = np.zeros((2, nb_bins, nb_concs))
    
    #populate empty arrays with relevant data from input tables
    for bb, gate in enumerate(range(1, nb_bins+1)):
        for cc, conc in enumerate(flow_table.concentration.unique()):
            #if sample doesn't exist because there were no events in that bin to send for ngs populate with nan
            if f"{construct}-{conc}-{gate}" in counts_table.columns:
                counts[bb, :, cc] = counts_table[f"{construct}-{conc}-{gate}"]
            else:
                counts[bb, :, cc] = np.nan 
            cells[bb, cc] = flow_table[(flow_table.concentration == conc) 
                                        & (flow_table.bin == gate)]["cell count"].iloc[0] 
            meanfluor[bb, cc] = flow_table[(flow_table.concentration == conc)
                                      & (flow_table.bin == gate)]["mean_log10_APC"].iloc[0]
            stdfluor[bb, cc] = flow_table[(flow_table.concentration == conc)
                                     & (flow_table.bin == gate)]["std_log10_APC"].iloc[0]

    print('text')
    print(counts)
    print(cells)
    print(meanfluor)
    print(stdfluor)

    print('probas test')
    probas = counts / (counts.sum(axis=1)[:, None, :])
    print(probas)
    probas = counts / (counts.sum(axis=1)[:, None, :]) * cells[:, None, :] #calc fraction of total reads in bin for seq s, scaled by num sorted cells in that bin
    print(probas.shape)
    print(probas)
    #good up to here 20240421 2:30 pm
    probas = probas / probas.sum(axis=0)[None, :, :] #normalized over four bins for each concentration
    mean_log10_fluor = (probas * meanfluor[:, None, :]).sum(axis=0)
    #print(mean_log10_fluor)
    std_log10_fluor = np.sqrt((stdfluor[:, None, :]**2 * probas**2
                               + meanfluor[:, None, :]**2 * probas**2 / (1e-22 + counts)).sum(axis=0)) #where is 1e-22 coming from
    #print(std_log10_fluor)

    concentrations = -concentrations #why invert values? was throwing error before
    # fit of Kd
    Kds, A, B, err, cov = np.zeros((5, nb_genos))
    for s in range(nb_genos):
        notnanindex = [ii for ii in range(nb_concs)
                       if not np.isnan(mean_log10_fluor[s, ii] + std_log10_fluor[s, ii])]
        if len(notnanindex) < 4: #i think this is filtering out an genos which have <4 data points for fluor
            Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
        # not enough reads
        if np.sum(counts.sum(axis=0) > min_num_counts) < 4: #setting min number of counts to 20 based on snakemake config file
            Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5 #filtering out genos with fewer than min counts for less than four samples
        else:
            #print(concentrations[notnanindex])
            #print(mean_log10_fluor[s, notnanindex])
            # Kds[s], A[s], B[s], err[s], cov[s] = extractKd(concentrations[notnanindex],
            #                                                mean_log10_fluor[s,
            #                                                                 notnanindex],
            #                                                std_log10_fluor[s, notnanindex])
            Kds[s], A[s], B[s], err[s], cov[s] = extractKd(concentrations[notnanindex],
                                                            mean_log10_fluor[s,
                                                                             notnanindex])

    return Kds, A, B, err, cov, mean_log10_fluor, std_log10_fluor, concentrations

test_concentrations = np.array([1.00e+05, 3.16e+04, 1.00e+04, 3.16e+03, 1.00e+03, 3.16e+02, 1.00e+02,
  3.16e+01, 1.00e+01, 3.16e+00, 1.00e+00, 3.16e-01])
# # test_std_log10 = np.array([2.32810988, 3.90291414, 2.73507449, 2.32034969, 2.21755796, 2.71159665,
# #  2.42383509, 2.52671346])
# test_mean_log10 = np.array([104155, 89253, 73857, 71115, 55644, 32912,
#   16179, 5989, 2148, 1135])
# test_mean_log10 = np.log10(test_mean_log10)

# #need to fix nans so can get as many concentrations as possible for each bin not capped at 8 concentrations
#test_Kds, test_A, test_B, test_err, test_cov, test_mean, test_std, test_conc = compute_Kds(flow_table, counts_table)
# print(test_Kds)                             
# #plot a subset of the data
# #need to rework
#df_subset = counts_table.sample(n=5)
df_subset = counts_table[counts_table.design == '08644_ALFA_53_7_min_2_53_bm01_loop']

# xs = np.linspace(-15, -3, 100)


# xlim = (-max(test_concentrations)-0.5,
#         -min(test_concentrations)+0.5)
# fig, ax = plt.subplots(4, 2, figsize=(10, 10), sharex=True)
# for ii, ax in enumerate(ax.flatten()):
#    #ax.errorbar(x=concentrations,
#                 #y=mlog10[ii, :], yerr=slog10[ii, :], label=df_subset.geno.iloc[ii])
#     ax.plot(xs, sigmoid(xs, -test_Kds, test_A, test_B))
#     ax.set_xlim(xlim)
#     ax.set_xlabel("$\log_{10}(\mathrm{Concentration})$")
#     ax.set_ylabel("Est. Mean Fluorescence")
# plt.show()


# Compute Kds for the full dataset
#set up option to run this if above looks good

df = pd.DataFrame()
Kds, As, Bs, errs, covs, mean_log10_PE, std_log10_PE, concs = compute_Kds(
    flow_table, df_subset)
print('output')
print(df_subset)
print(mean_log10_PE)
print(mean_log10_PE[0])
sns.scatterplot(x=test_concentrations, y=mean_log10_PE[0])
plt.xscale('log')
plt.show()
plt.clf()
df["log10Kd"] = Kds
df["A"] = As
df["B"] = Bs
df["r2"] = errs
df["sigma"] = covs
print(df)

# save the Kds values
df[["geno", "log10Kd", "A", "B", "r2", "sigma"] +
   [f"mean_log10PE{cc}" for cc in range(0, mean_log10_PE.shape[1])] +
   [f"std_log10PE{cc}" for cc in range(0, mean_log10_PE.shape[1])]].to_csv(
    f'{outdir}/all_kds.csv',
    sep="\t", index=False)


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


