import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.optimize
from scipy.optimize import curve_fit
import os

#From A Phillips and ZY Weinberg. Fits Kd to tite-seq data using mean bin approach and plots estimated mean log fluorescence and fit curve. 
#See: https://github.com/amphilli/CH65-comblib/tree/main/Kd_Inference

####USER DEFINED VARIBLES####
counts_table_path = '/Users/stephaniecrilly/test/count_table.csv'
counts_table = pd.read_csv(counts_table_path, sep='\t')

flow_table_path = '/Users/stephaniecrilly/test/flow_table.csv'
flow_table = pd.read_csv(flow_table_path)

outdir = '/Users/stephaniecrilly/Kortemme_lab/switch_seq'
####

####FUNCTIONS####

def extractKd(concentrations, bins_mean, bins_std):
    """
        Arguments: concentrations and exponential mean bin numbers (bin number: 1, 2, 3, 4)
        Return: -log10(Kd), and the r^2
    """
    popt, pcov = scipy.optimize.curve_fit(sigmoid, concentrations,
                                          bins_mean,
                                          p0=[(-9), 10**(4), 10**(2)],
                                          sigma=bins_std, absolute_sigma=True,
                                          bounds=[(-12,
                                                   100,
                                                   1),
                                                  (-7,
                                                   1000000,
                                                   100000)],
                                          maxfev=400000)
    return(-1*popt[0], popt[1], popt[2], 1 - np.sum((sigmoid(concentrations, *popt) - bins_mean)**2)/np.sum((bins_mean - bins_mean.mean())**2), np.sqrt(np.diag(pcov))[0])


def sigmoid(c, Kd, A, B):
    return np.log10(A * (10**c/((10**c)+(10**Kd))) + B)


def compute_Kds(sample_info, fluor, df):
    sample_info = sample_info[~sample_info.concentration_float.isna()].copy()
    nb_bins = sample_info.bin.nunique()
    nb_concs = sample_info.concentration.nunique()
    concentrations = sample_info.concentration_float.unique()
    nb_genos = len(df)
    probas = np.zeros((nb_bins, nb_genos, nb_concs))
    counts = np.zeros((nb_bins, nb_genos, nb_concs))
    cells = np.zeros((nb_bins, nb_concs))
    meanfluor, stdfluor = np.zeros((2, nb_bins, nb_concs))
    for gate, conc in zip(sample_info['bin'],
                                      sample_info['concentration']):
        counts[gate-1, :, conc-1] = df[f"r1-hs-{conc}-{gate}"]
        cells[gate-1, conc-1] = sample_info[(sample_info.concentration == conc)
                                    & (sample_info.bin == gate)]["cell count"].iloc[0]
        meanfluor[gate-1, conc-1] = fluor[(fluor.concentration == conc)
                                  & (fluor.bin == gate)]["mean_log10_PEs"].iloc[0]
        stdfluor[gate-1, conc-1] = fluor[(fluor.concentration == conc)
                                 & (fluor.bin == gate)]["std_log10_PEs"].iloc[0]

    probas = counts / (counts.sum(axis=1)[:, None, :]) * cells[:, None, :]
    probas[np.isnan(probas)] = 0.
    probas = probas / probas.sum(axis=0)[None, :, :]
    mean_log10_fluor = (probas * meanfluor[:, None, :]).sum(axis=0)
    std_log10_fluor = np.sqrt((stdfluor[:, None, :]**2 * probas**2
                               + meanfluor[:, None, :]**2 * probas**2 / (1e-22 + counts)).sum(axis=0))

    # replace the "0" concentration by an arbitrary large value (here 20) and invert the values
    concentrations = -concentrations
    concentrations[concentrations ==
                   0] = concentrations[concentrations == 0] - 20
    # fit of Kd
    Kds, A, B, err, cov = np.zeros((5, nb_genos))
    for s in range(nb_genos):
        notnanindex = [ii for ii in range(nb_concs)
                       if not np.isnan(mean_log10_fluor[s, ii] + std_log10_fluor[s, ii])]
        if len(notnanindex) < 4:
            Kds[s], A[s], B[s], err[s], cov[s] = [np.nan]*5
        # not enough reads
        Kds[s], A[s], B[s], err[s], cov[s] = extractKd(concentrations[notnanindex],
                                                           mean_log10_fluor[s,
                                                                            notnanindex],
                                                           std_log10_fluor[s, notnanindex])

    return Kds, A, B, err, cov, mean_log10_fluor, std_log10_fluor, concentrations

####

#make outdir
if not os.path.exists(outdir):
   os.makedirs(outdir)

#Read in counts and flow data
counts = (pd.read_csv(counts_table_path, sep='\t')
          .rename({'design':'geno'}, axis=1)
          .replace([np.inf, -np.inf], np.nan)
          .dropna())
fluor = (pd.read_csv(flow_table_path, sep=',', index_col=0)
         .rename({'sample_id': 'sample_name',
                  'mean_log10_APC': 'mean_log10_PEs',
                  'std_log10_APC': 'std_log10_PEs'}, axis=1)
         .replace([np.inf, -np.inf], np.nan)
         .dropna())
sample_info = fluor[['sample_name', 'concentration', 'concentration_float', 'bin', 'cell count']]
sample_info.concentration_float = -1*np.log10(sample_info.concentration_float*1e-12)
fluor = fluor[['sample_name', 'concentration', 'bin', 'mean_log10_PEs', 'std_log10_PEs']]

#compute Kds
Kds, As, Bs, errs, covs, mean_log10_PE, std_log10_PE, concs = compute_Kds(
        sample_info, fluor, counts)

#write file with calculated Kds for each design
fits = pd.DataFrame({'Kds':Kds, 'As':As, 'Bs':Bs, 'errs':errs, 'covs':covs})
counts_and_fits = pd.concat([counts, fits], axis=1)

counts_and_fits = counts_and_fits[['geno', 'Kds', 'As', 'Bs', 'errs', 'covs']]

#convert Kds from -log10(Kd) to Kd
counts_and_fits['Kds_M'] = 10**(-counts_and_fits['Kds'])

counts_and_fits.to_csv(f'{outdir}/Kds.csv', index=False)
