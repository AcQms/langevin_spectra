import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
import sys


def chi_mem_rescaled(w,k,t1,t2):
    kernel = 1/(1j*w*t2+1)
    return 1/(k-t1*w**2+1j*w*kernel)



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
t2=float(sys.argv[2])
a=float(sys.argv[3])
b=float(sys.argv[4])

data_mem=np.zeros([2,int(50000000)])
window_weight=0

f,spec=np.load('/scratch/voncanaa93/results_mem_versuch/results_a=%.3f_b=%.2f_t1=%.2f_t2=%.2f/trial.npy'%(a,b,t1,t2),allow_pickle=True)
data_mem=f,f**2*spec

if t1 and t2 == 0.01:
    window_weight=100
elif t1 == 0.01:
    window_weight=30
elif t1 == 0.1:
    window_weight=75
elif t1 == 1.0:


dw=data_mem[0][1]-data_mem[0][0] 
window_length=Peak(data_mem[0],data_mem[1])[2]/window_weight 
sigma=window_length/dw 
gauss_filtered = gaussian_filter1d(data_mem[1],sigma) 
filtered_shape = Peak(data_mem[0],gauss_filtered) 
 
analytic_mem = -2*data_mem[0]*chi_mem_rescaled(data_mem[0],1,t1,t2).imag 
analytic_shape = Peak(data_mem[0],analytic_mem) 
 
w_shift_rel = (analytic_shape[1]-filtered_shape[1])/analytic_shape[1] 
fwhm_shift_rel = (analytic_shape[2]-filtered_shape[2])/analytic_shape[2] 
 
np.save('mem_peak_params_sim_t1=%.2f_t2=%.2f_a=%.3f_b=%.2f.npy'%(t1,t2,a,b),filtered_shape) 
np.save('mem_peak_params_ana_t1=%.2f_t2=%.2f.npy'%(t1,t2),analytic_shape) 
np.save('w_shift_rel_t1=%.2f_t2=%.2f_a=%.3f_b=%.2f.npy'%(t1,t2,a,b),w_shift_rel) 
np.save('fwhm_shift_rel_t1=%.2f_t2=%.2f_a=%.3f_b=%.2f.npy'%(t1,t2,a,b),fwhm_shift_rel) 
 
                




