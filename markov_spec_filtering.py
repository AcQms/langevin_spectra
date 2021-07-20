import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
import sys
import os

def chi_rescaled(w,k,t1):
    return 1/(k-t1*w**2+1j*w)



def Peak(x,f):
    """
    Peak(x,f) returns (amplitude, peak position, FWHM)
    
    """

    results = np.zeros([3])
    f_max = np.amax(f)
    w0 = x[np.argmax(f)]

    f_fwhm = .5*f_max
    f_fwhm_index = np.where(f == f_max)[0][0]
    w_fwhm1 = x[np.argmin(abs(f[:f_fwhm_index]-f_fwhm))]
    w_fwhm2 = x[f_fwhm_index + np.argmin(abs(f[f_fwhm_index:]-f_fwhm))]
    width = w_fwhm2-w_fwhm1
    results[0] = f_max
    results[1] = w0
    #results[2] = w_fwhm1
    #results[3] = w_fwhm2
    results[2] = width

    return results

t1=float(sys.argv[1])
a=float(sys.argv[2])
b=float(sys.argv[3])

window_weight=0

data = np.zeros([2,int(50000000)])


f,spec=np.load('/scratch/voncanaa93/results_markov_versuch/results_a=%.3f_b=%.2f_t1=%.2f/trial_markov.npy'%(a,b,t1),allow_pickle=True)
data=f,f**2*spec

if t1==100.:
    window_weight=10
elif t1==10.:
    window_weight=30
elif t1==1.:
    window_weight=80
elif t1==0.1:
    window_weight=75
elif t1==0.01:
    window_weight=200

dw=data[0][1]-data[0][0]
window_length=Peak(data[0],data[1])[2]/window_weight
sigma=window_length/dw
gauss_filtered = gaussian_filter1d(data[1],sigma)
filtered_shape = Peak(data[0],gauss_filtered)

analytic = -2*data[0]*chi_rescaled(data[0],1,t1).imag
analytic_shape = Peak(data[0],analytic)

w_shift_rel = (analytic_shape[1]-filtered_shape[1])/analytic_shape[1]
fwhm_shift_rel = (analytic_shape[2]-filtered_shape[2])/analytic_shape[2]

np.save('markov_peak_params_sim_t1=%.2f_a=%.3f_b=%.2f.npy'%(t1,a,b),filtered_shape)
np.save('markov_peak_params_ana_t1=%.2f.npy'%(t1),analytic_shape)
np.save('markov_w_shift_rel_t1=%.2f_a=%.3f_b=%.2f.npy'%(t1,a,b),w_shift_rel)
np.save('markov_fwhm_shift_rel_t1=%.2f_a=%.3f_b=%.2f.npy'%(t1,a,b),fwhm_shift_rel)
                                                                                                                                                                                                                                                                                                                                