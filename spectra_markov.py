import random
import numpy as np
import math
import numba
from numba import jit
from scipy.ndimage.filters import gaussian_filter1d
import sys

@jit(nopython=True,nogil=True)
def integrateRK_rescaled(nsteps,stride,dt,x0,v0,t1,a,b,c):

    x = np.zeros((nsteps,),dtype=np.float64)
    x[0] = x0
    v = np.zeros((nsteps,),dtype=np.float64)
    v[0] = v0
    f = math.sqrt(2/dt)
    #xi = math.sqrt(1/td)

    for i in range(nsteps):
        for j in range(stride):

            F = f*random.gauss(0.0,1)

            kx1=dt*v0
            kv1=dt*(-v0-a*x0**3-b*x0**2-c*x0+F)*1/t1

            x1=x0+kx1/2
            v1=v0+kv1/2

            kx2=dt*v1
            kv2=dt*(-v1-a*x1**3-b*x1**2-c*x1+F)*1/t1

            x2=x1+kx2/2
            v2=v1+kv2/2

            kx3=dt*v2
            kv3=dt*(-v2-a*x2**3-b*x2**2-c*x2+F)*1/t1

            x3=x2+kx3/2
            v3=v2+kv3/2

            kx4=dt*v3
            kv4=dt*(-v3-a*x3**3-b*x3**2-c*x3+F)*1/t1

            x0+=(kx1+2*kx2+2*kx3+kx4)/6
            v0+=(kv1+2*kv2+2*kv3+kv4)/6

        x[i] = x0
        v[i] = v0

    return (x,v)


def spectrum1(t,a,b=None,subtract_mean=False):

    dt = t[1]-t[0]
    meana = int(subtract_mean)*np.mean(a)
    a2 = a-meana
    f, fra = FT(t,a2)

    if b is None:
        sf = np.conj(fra)*fra

    else:
        meanb = int(subtract_mean)*np.mean(b)
        b2 = b-meanb
        f, frb = FT(t,b2)
        sf = np.conj(fra)*frb

    sff = sf[f>=0]
    ff = f[f>=0]
    spec = sff/(len(t)*dt)
    return ff,spec



def FT(t,x):

    a, b = np.min(t), np.max(t)
    dt = t[1]-t[0]
    #if (abs((t[1:]-t[:-1] - dt)) > 1e-13).any():
        #print(np.max( abs(t[1:]-t[:-1])))
        #raise RuntimeError("Time series not equally spaced!")
    N = len(t)
    # calculate frequency values for FT
    k = np.fft.fftshift(np.fft.fftfreq(N,d=dt)*2*np.pi)
    # calculate FT of data
    xf = np.fft.fftshift(np.fft.fft(x))
    xf2 = xf*(b-a)/N*np.exp(-1j*k*a)
    return k, xf2

def Peak(x,f):

    results = np.zeros([3])
    f_max = np.max(f)
    w0 = x[np.argmax(f)]
    f_fwhm =.5*f_max
    f_fwhm_index = np.where(f == f_max)[0][0]
    w_fwhm1 = x[np.argmin(abs(f[:f_fwhm_index]-f_fwhm))]
    w_fwhm2 = x[f_fwhm_index + np.argmin(abs(f[f_fwhm_index:]-f_fwhm))]
    width = w_fwhm2 - w_fwhm1
    results[0] = f_max
    results[1] = w0
    results[2] = width

    return results



#@jit(nopython=True,nogil=True)
def Mean_Spectrum(N,nsteps,stride,t1,a,b,c):

    #shape = np.zeros([N,3])
    td = 1
    specs = np.zeros([int(nsteps/2)])
    dt = 1e-2*min(t1,td)
    mean_malus = 0

    time = stride*np.linspace(0,dt*nsteps,nsteps)
    for i in range(N):
        x0 = random.gauss(0.0,1)
        v0 = random.gauss(0.0, 1/np.sqrt(t1))
        x = integrateRK_rescaled(nsteps,stride,dt,x0,v0,t1,a,b,c)[0]
        f,spectrum = spectrum1(time,x,None,True)
        if any(np.isnan(spectrum.real))==False:
            pass #shape[i] = Peak(f,spectrum.real)
        else:
            mean_malus += 1
            #pass
        specs += np.nan_to_num(spectrum.real,copy=True)

    mean_spec = 1/(N-mean_malus) * specs

    return (f,mean_spec) #,shape)


N=int(sys.argv[1])
stride=int(sys.argv[2])
nsteps=int(sys.argv[3])
t1=float(sys.argv[4])
a=float(sys.argv[5])
b=float(sys.argv[6])

c=1
#ww=0

results=Mean_Spectrum(N,nsteps,stride,t1,a,b,c)

#if t1==100.:
#    ww=10
#elif t1==10.:
#     ww=30
#elif t1==1.:
#    ww=80
#elif t1==.1:
#    ww=75
#elif t1==.01:
#    ww=200

#dw=results[0][1]-results[0][0]
#window_length=Peak(results[0],results[0]**2*results[1])[2]/ww
#sigma=window_length/dw
#smooth=gaussian_filter1d(results[0]**2*results[1],sigma)




np.save('trial_markov.npy',(results[:2]))#,smooth))
#np.save('shape_raw.npy', results[2])                                                                                                                                                                                                                                           