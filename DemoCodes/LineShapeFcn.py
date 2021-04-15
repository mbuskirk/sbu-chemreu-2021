#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 12:49:45 2021

@author: tka3
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand, randn
from scipy.integrate import cumtrapz
from scipy.signal import butter, filtfilt

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams.update({'figure.figsize': (10,6)})

#%%
# Define Fourier Transform Pair ========================================

def dft(t,f):
	"""
	tomlab.dft
	
	This function performs the discrete fourier transform of the 
	function f(t) evalulated at times t. f amd t are numpy arrays.
	t is assumed to be an evenly spaced grid of times that is a power 
	of 2 in length. It uses numpy's fft code, but also takes into 
	account for the spectral phase artifacts of FFT due to the fact 
	that the first time point is likely not at zero. This is done by 
	multiplying the FFT data by exp(-j*w*t[1])
	
	Returns a tuple with the frequency axis as well:
	(w,F), where w is the (angular) frequency axis and F(w) is the
	Fourier transform.
	"""
	dt = t[1]-t[0]
	
	F = np.fft.fft(f) #fft f(t)
	F = np.fft.fftshift(F) #fft shift to put negative frequency components first in the array
	w = 2*np.pi*np.fft.fftfreq(len(F),dt) #get frequency axis
	w = np.fft.fftshift(w) #fft shift returned frequencies as well
	F = np.exp(-(0+1j)*t[0]*w)*F #multiply by linear spectral phase to compensate time offset, since the first element is not 0
	
	return (w,F)  
	
def idft(w,F,t0 = None):
	"""
	tomlab.idft
	
	This function performs the inverse operations of dft. See documentation
	for tomlab.dft(t,f) for more information. 
	
	Returns a tuple with the time axis:
	(t,f), where t is the time and f(t) is the function of time.
	
	There is ambiguity in the time grid with inverse FFT in general.
	The optional argument t0 lets you specify the first element of the
	time array. If it is left blank, then it is an evenly spaced grid
	around t =0
	"""
	dw = w[1]-w[0]
	N = len(w);
	dt = 2*np.pi/((w[2]-w[1])*N)
	F = np.fft.ifftshift(F) #inversefftshift both f and w
	w = np.fft.ifftshift(w)
	
	if t0 is None:
		t = np.linspace(-(N-1)*dt/2, (N-1)*dt/2,N)		
	else:
		t = t0 + np.linspace(0,(N-1)*dt,N)
			
	f = np.fft.ifft(np.exp((0+1j)*w*t[0])*F)
	
	return (t,f)

# %%
# Define frequency fluctuation problem
  

w21 = 10 # transition frequency [ar]
BW1 = 1  # Bandwidth of fluctuations (not precisely defined)
DW_rawamp = 100 # Amplitude of raw (unfiltered) fluctuations
Nmult = 100     # Number of runs to do between displayed molecules.
nplot = 5       # Number of molecular trajectories to plot.
N = Nmult*nplot # Total number of molecules to do
BWN = 0*0.5     # Bandwidth of inhomogeneous broadening (not precisely defined)
delta = BWN*randn(N) # Array Random inhomogenous detunings sampled from from a normal distribution.

# Initialize t vector and empty frequency fluctuation and coherence
trun = np.linspace(0,100,2**16)
dw = np.zeros( (N,len(trun)) )
rho21 = np.zeros( (N,len(trun)) ) + 1j*np.zeros( (N,len(trun)) ) 

#Loop through molecules
for n in range(N):
    dwraw = DW_rawamp*(rand(len(trun))-0.5)  # array of unfiltiered frequency fluctuations - completely Markovian  
    
    # Apply Gaussian filter to raw fluctuations
    w,F = dft(trun,dwraw)   
    F_filt = F*np.exp(-w**2/BW1**2)
    t,dw1 = idft(w,F_filt,t0 = trun[0])
    
    #Add inhomogeneous broadening
    dw[n,:] = np.real(dw1) + delta[n]

    #Calculate coherence
    rho21[n,:] = np.exp(-1j*cumtrapz(w21+dw[n,:],t, initial = 0))
    
    # Make it so all coherences start in phase.
    if n != 0:
        rho21[n,:] = (rho21[0,0]/rho21[n,0])*rho21[n,:] # make all initially in phase

# Sum frequency fluctuations - i.e. perform ensemble average        
DW = np.mean(dw, axis = 0) 

# Calculate histogram, and FFCF via autocorrelation.    
dw_hist,be = np.histogram(dw, bins = 50)
bc = (be[1:] + be[0:-1])/2
ACF = np.correlate(DW,DW, mode = 'full')
ACF = ACF[round(len(ACF)/2)-1:]

# Calculate absorption spectrum as Fourier transform of coherence of the ensemble average.
w,S = dft(t,np.mean(np.conj(rho21),axis = 0))

# %%
# Make Plots ===========================================================
fig, (ax1,ax2) = plt.subplots(nrows = 1,ncols = 2, sharex = False)
for n in range(nplot+1):
    if n != nplot:
        ax1.plot(t,dw[Nmult*n,:])
    else :
        ax1.plot(t,dw[Nmult*n-1,:])

ax1.plot(t,DW, color = 'black',linewidth = 2)
ax1.grid()
ax1.set_xlabel('$t$')
ax1.set_ylabel('$\delta \omega(t)$')

ax2.plot(dw_hist,bc, color = 'black')
ax2.set_ylim(ax1.get_ylim())
ax2.grid()

# ax3.plot(t,ACF)
# ax3.set_xlim( (0,10) )
# ax3.grid()

fig, (ax4,ax5, ax6) = plt.subplots(3, sharex = False)
 
for n in range(nplot):
    if n != nplot:
        ax4.plot(t,np.real(rho21[Nmult*n,:]))
    else :
        ax4.plot(t,np.real(rho21[Nmult*n-1,:]))
        
ax4.set_xlim( (0,10/BW1) )
ax4.set_xlabel('t')

ax5.plot(t,np.real(np.mean(rho21, axis = 0)))
ax5.set_ylabel('$\\rho_{21}(t)$')
ax5.set_xlim( (0,10/BW1) )
ax5.set_xlabel('t')

ax6.plot(w,np.real(S))
ax6.set_xlim((w21-4,w21+4))
ax6.grid()
ax6.set_xlabel('$\omega$')
ax6.set_ylabel('Absorption Spectrum')

plt.show()

