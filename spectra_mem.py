
import random
import numpy as np
import math
import numba
from numba import jit
from scipy.ndimage.filters import gaussian_filter1d
import sys

@jit(nopython=True,nogil=True)
def integrate_1_exp_rescaled(nsteps,stride,dt,t1,t2,x0,v0,y,a,b,c):
    """ nsteps = number of steps in terms of dt
        dt  = timestep
        m = mass in atomic units
        gamma = friction coefficient in u/ps
        tgamma = memory time in ps
        x0 = initial position
        v0 = initial velocity
        y = initial force
        k = hook s constant
        kbT = k_B T default 2.494 in kJ/mol
    """
    tm = t1 ; tg = t2
    x=np.zeros((nsteps,),dtype=np.float64)
    v=np.zeros((nsteps,),dtype=np.float64)
    #xi_sigma1=math.sqrt(2*tau_d)
    xi_factor=math.sqrt(2/dt)
    x[0]=x0
    xx=x0
    v[0]=v0
    vv=v0


    for i in range(nsteps):
        for j in range(stride):
            xi=xi_factor*random.gauss(0.0, 1)

            kx1=dt*vv
            kv1=dt*((y-xx)/tg-(c*xx+b*xx**2+a*xx**3))/tm
            k1=-dt*((y-xx)/tg-xi)

            x1=xx+kx1/2
            v1=vv+kv1/2
            y1=y+k1/2

            kx2=dt*v1
            kv2=dt*((y1-x1)/tg-(c*x1+b*x1**2+a*x1**3))/tm
            k2=-dt*((y1-x1)/tg-xi)


            x2=xx+kx2/2
            v2=vv+kv2/2
            y2=y+k2/2


            kx3=dt*v2
            kv3=dt*((y2-x2)/tg-(c*x2+b*x2**2+a*x2**3))/tm
            k3=-dt*((y2-x2)/tg-xi)

            x3=xx+kx3
            v3=vv+kv3
            y3=y+k3

            kx4=dt*v3
            kv4=dt*((y3-x3)/tg-(c*x3+b*x3**2+a*x3**3))/tm
            k4=-dt*((y3-x3)/tg-xi)

            xx+=(kx1+2*kx2+2*kx3+kx4)/6
            vv+=(kv1+2*kv2+2*kv3+kv4)/6
            y+=(k1+2*k2+2*k3+k4)/6


        x[i]=xx
        v[i]=vv

    return x,v


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


#@jit(nopython=True,nogil=True)
def Mean_Spectrum(N,nsteps,stride,t1,t2,a,b,c):

    #shape = np.zeros([N,3])
    specs = np.zeros([int(nsteps/2)])
    td = 1
    dt = 1e-2*min(t1,t2,td)
    mean_malus = 0

    time = stride*np.linspace(0,dt*nsteps,nsteps)
    for i in range(N):
        x0 = random.gauss(0.0,1)
        v0 = random.gauss(0.0, 1/np.sqrt(t1))
        y = np.random.normal(0.0,1/np.sqrt(t2))
        x = integrate_1_exp_rescaled(nsteps,stride,dt,t1,t2,x0,v0,y,a,b,c)[0]
        #x = integrateRK_rescaled(nsteps,stride,dt,x0,v0,t1,a,b,c)[0]
        f, spectrum = spectrum1(time,x,None,True)
        if any(np.isnan(spectrum.real))==False:
            pass #shape[i] = Peak(f,spectrum.real)
        else:
            mean_malus += 1
            #pass
        specs += np.nan_to_num(spectrum.real,copy=True)
    #specs_no_nans=specs[~np.isnan(specs).any(axis=1)]
    mean_spec = 1/(N-mean_malus) * specs

    return (f,mean_spec) #,shape)

N=int(sys.argv[1])
stride=int(sys.argv[2])
nsteps=int(sys.argv[3])
t1=float(sys.argv[4])
t2=float(sys.argv[5])
a=float(sys.argv[6])
b=float(sys.argv[7])

c=1
#ww=0
results=Mean_Spectrum(N,nsteps,stride,t1,t2,a,b,c)

#if t1 and t2 ==.01:
#    ww=250
#elif t1 == .01:
#    ww=30
#elif t1 == .1:
#    ww=60
#elif t1 == 1:
#    ww=20
#elif t1 == 10:
#    ww=20    
#elif t1 == 100:
#    ww=10

#df = results[0][1]-results[0][0]
#window_length=Peak(results[0],results[0]**2*results[1])[2]/ww
#sigma=window_length/df
#smooth=gaussian_filter1d(f**2*spec[0],sigma)
#smooth=gaussian_filter1d(results[0]**2*results[1],sigma)

np.save('trial.npy',(results[:2])) #,smooth))
#np.save('shape_raw',results[2])
#np.save('smooth',(results[0],smooth))

