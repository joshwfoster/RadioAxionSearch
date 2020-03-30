import os
import numpy as np
from scipy import ndimage, interpolate


###################################################################
###   Get the Expected Flux Density from a Calibration Source   ###
###################################################################

def get_expected(freq, source):

    '''
    Calculate the frequency-dependent expected flux
    density for a calibration source

    :param freq: the radio frequency [GHz]
    :param source: the name of the calibration source ['3C286' or '3C48']

    :returns: the expected flux density [Jy]
    '''

    params_3C286 = np.array([1.2481, -0.4507, -0.1798, 0.0357])
    params_3C48  = np.array([1.3253, -0.7553, -0.1914,  0.0498])

    if source == '3C48':
        a0 = params_3C48[0]
        a1 = params_3C48[1]
        a2 = params_3C48[2]
        a3 = params_3C48[3]
    elif source == '3C286':
        a0 = params_3C286[0]
        a1 = params_3C286[1]
        a2 = params_3C286[2]
        a3 = params_3C286[3]
    else:
	print 'Invalid Source'
        return None

    return 10**(a0 + a1*np.log10(freq) + a2 * np.log10(freq)**2 + a3 * np.log10(freq)**3)

##########################################################################
###   Load the Data From the Data Dictionary and Make the Data Stack   ###
##########################################################################

def get_stack(data_dir, data_tag):
    
    data_stack = np.zeros((4, 0))

    freqs = np.zeros((0))
    
    xx_sig_accepted = np.zeros((0))
    xx_sig_on = np.zeros((0))
    xx_sig_off = np.zeros((0))
    xx_ref_accepted = np.zeros((0))
    xx_ref_on = np.zeros((0))
    xx_ref_off = np.zeros((0))

    yy_sig_accepted = np.zeros((0))
    yy_sig_on = np.zeros((0))
    yy_sig_off = np.zeros((0))
    yy_ref_accepted = np.zeros((0))
    yy_ref_on = np.zeros((0))
    yy_ref_off = np.zeros((0))
                               
        
    for filename in os.listdir(data_dir):
        if data_tag in filename:
            print filename
            
            data = np.load(data_dir + filename)['arr_0'].item() 
            
            XX = data['XX_Out']
            YY  = data['YY_Out']
            
            xx_tcal = data['XX_TCal']
            yy_tcal = data['YY_TCal']

            freqs = np.append(freqs, data['freqs'])
                               
            xx_sig_accepted = np.append(xx_sig_accepted, XX[0])
            xx_sig_on = np.append(xx_sig_on, XX[1])
            xx_sig_off = np.append(xx_sig_off, XX[2])
            xx_ref_accepted = np.append(xx_ref_accepted, XX[3])
            xx_ref_on = np.append(xx_ref_on, XX[4])
            xx_ref_off = np.append(xx_ref_off, XX[5])


            yy_sig_accepted = np.append(yy_sig_accepted, YY[0])
            yy_sig_on = np.append(yy_sig_on, YY[1])
            yy_sig_off = np.append(yy_sig_off, YY[2])
            yy_ref_accepted = np.append(yy_ref_accepted, YY[3])
            yy_ref_on = np.append(yy_ref_on, YY[4])
            yy_ref_off = np.append(yy_ref_off, YY[5])
            
            
    sort_order = np.argsort(freqs)
    freqs = freqs[sort_order]
    
    xx_sig_accepted = xx_sig_accepted[sort_order]  
    xx_ref_accepted = xx_ref_accepted[sort_order]
    yy_sig_accepted = yy_sig_accepted[sort_order]
    yy_ref_accepted = yy_ref_accepted[sort_order]

    xx_sig_on = xx_sig_on[sort_order]
    xx_sig_off = xx_sig_off[sort_order]
    xx_sig = (xx_sig_on + xx_sig_off) /2.
    
    xx_ref_on = xx_ref_on[sort_order]
    xx_ref_off = xx_ref_off[sort_order]
    xx_ref = (xx_ref_on + xx_ref_off) / 2.
    
    yy_sig_on = yy_sig_on[sort_order]
    yy_sig_off = yy_sig_off[sort_order]
    yy_sig = (yy_sig_on + yy_sig_off) /2.
 
    yy_ref_on = yy_ref_on[sort_order]
    yy_ref_off = yy_ref_off[sort_order]
    yy_ref = (yy_ref_on + yy_ref_off) /2.

    xx_ref_tsys = xx_tcal * ( xx_ref_off / (xx_ref_on - xx_ref_off) + .5)
    yy_ref_tsys = yy_tcal * ( yy_ref_off / (yy_ref_on - yy_ref_off) + .5) 
    ref_tsys = (xx_ref_tsys + yy_ref_tsys) / 2

    xx_sig_tsys = xx_tcal * ( xx_sig_off / (xx_sig_on - xx_sig_off) + .5)
    yy_sig_tsys = yy_tcal * ( yy_sig_off / (yy_sig_on - yy_sig_off) + .5)
    sig_tsys = (xx_sig_tsys + yy_sig_tsys) / 2


    xx_ta = (xx_sig - xx_ref) / xx_ref * ndimage.median_filter(xx_ref_tsys, size = 31)
    yy_ta = (yy_sig - yy_ref) / yy_ref * ndimage.median_filter(yy_ref_tsys, size = 31)


    ta = (xx_ta + yy_ta) / 2
    
    total_temp = ref_tsys + sig_tsys + ta
    
    sig = (xx_sig + yy_sig)/2
    ref = (xx_ref + yy_ref)/2
    
    data_stack = np.vstack((freqs, ta, sig, ref, xx_sig_accepted, xx_ref_accepted, yy_sig_accepted, yy_ref_accepted))
    return data_stack


##################################################
###   Load the Data at Specified Downbinning   ###
##################################################

def downsample_stack(stack, downsample, shift = 0):
    freqs = stack[0][shift:]
    ta = stack[1][shift:]
    sig = stack[2][shift:]
    ref = stack[3][shift:]
    xx_sig_accepted = stack[4][shift:]
    xx_ref_accepted = stack[5][shift:]
    yy_sig_accepted = stack[6][shift:]
    yy_ref_accepted = stack[7][shift:]
    
    sig_accepted = (xx_sig_accepted + yy_sig_accepted) / 2
    ref_accepted = (xx_ref_accepted + yy_ref_accepted) / 2
    
    max_index = len(freqs) / downsample * downsample
    num_intervals = len(freqs) / downsample
    
    out = np.zeros((6, num_intervals))
    
    for i, item in enumerate([freqs, ta, sig, ref, sig_accepted, ref_accepted]):
        item = np.mean(item[:max_index].reshape(num_intervals, downsample), axis = 1)
        out[i] = np.copy(item)
        
    return out


def load(data_dir, data_tag, downsample = 1, shift = 0):
    downsample_stack(get_stack(data_dir, data_tag), downsample, shift = shift)
