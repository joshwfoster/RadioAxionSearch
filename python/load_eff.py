import sys, os, h5py, corner
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats, signal, ndimage, interpolate

import astropy
from astropy.io import fits



##################################################
###   Load the Data at Specified Downbinning   ###
##################################################

def downsample_stack(stack, downsample, shift = 0):
    freqs = stack[0][shift:]
    sig = stack[1][shift:]
    ref = stack[2][shift:]
    
    max_index = len(freqs) / downsample * downsample
    num_intervals = len(freqs) / downsample
    
    out = np.zeros((3, num_intervals))
    
    for i, item in enumerate([freqs, sig, ref]):
        item = np.mean(item[:max_index].reshape(num_intervals, downsample), axis = 1)
        out[i] = np.copy(item)
        
    return out


def load(file, downsample = 1, shift = 0):
    stack = np.vstack((np.load(file)['arr_0'].item().values()))
    return downsapmle_stack(stack, downsample, shift = shift)
