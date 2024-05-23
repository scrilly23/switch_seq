
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import seaborn as sns

####FUNCTIONS####

#all assume no cooperativity, hill_coeff=1

#this function does not force curve through zero and does not assume top and bottom values
def funcHyp(x, b, a, c):
    return a + ( (b-a) / (1 + (c/x)) )

#this function is for normalized data, assuming top=1 and bottom=0
#assumes no cooperativity, hill_coeff=1
def funcHyp_norm(x, c):
    return ( (x) / (c+x)) 

def sigmoid(b_tot, Kd, p_tot, bg):
    '''
    Equation for a sigmoidal curve
    For curve fit with this function, 
    pass array of y vals as log10 transformed fluorescence
    paass array of x vals as -log10 transformed concentration of binder
    
    Parameters:
    b_tot: -log10 transformed binder concentration values (Molar)
    Kd: dissociation constant
    p_tot: total protein concentration (upper asymptote)
    bg: background signal (lower aymptote)

    '''
    return np.log(p_tot * (10**b_tot/((10**b_tot)+(10**Kd))) + bg) #equivalent to A.M.Phillips sigmoid function


def background_correct(input_df, d, bg_vals):
    '''
    Subtracts background signal from data.
    
    Parameters:
    input_df (df): Data frame of raw cytometry data
    d (str): Name of data column to background subtract
    bg_vals (list of floats/ints): List of values from control sample to subtract as background, should be in same sample order as df samples
     
    Returns:
    df: Returns original df with a new column "BG_values" and a new column "d_bg_corrected" with BG values subtracted from data 
    
    '''

    assert input_df.shape[0] == len(bg_vals), "Number of bg_vals does not match number of data points."

    new_df = input_df.copy()

    new_df['BG_values'] = bg_vals

    new_colname = d+'_bg_corrected'
    new_df[new_colname] = new_df[d] - new_df['BG_values']

    return new_df

def min_max_normalize(input_df, d):
    '''
    Performs min to max normalization on df column of interest.

    Parameters:
    input_df (df): Data frame of cytometry data
    d (str): Name of data column to normalize

    Returns:
    df: Returns original df with a new column
    '''
    new_df = input_df.copy()

    new_colname = d+'_norm'

    new_df[new_colname] = (new_df[d] - new_df[d].min()) / (new_df[d].max() - new_df[d].min())

    return new_df