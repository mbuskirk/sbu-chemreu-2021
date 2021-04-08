#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpticalBlochEquations.py

This script integrates the Optical Bloch equations and plots the results.

Also makes plots of the lineshape for various paramters!

Origin: CHE 525, S21 Light/Matter Interaction Lectures
Author: Tom Allison
"""

#%%
# Preliminaries ===============================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import matplotlib

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams.update({'figure.figsize': (10,6)})

# %%
# Define global variables

Omega = 1   #Everything is referenced to the Rabi frequency Omega, so might as well make it one!
T1inv = (1/1000)*Omega #spontaneous emission rate
T1 = 1/T1inv  
T2inv = (1/2)*T1inv + 0*Omega # total dephasing rate, including spontaneous emission.
T2 = 1/T2inv
Delta = 0*Omega #detuning
weq = -1

# Initial conditions
sigma_init = 0
w_init = weq

if T1*T2*Omega**2 > 1:
    trun = np.linspace(0,5*2*np.pi/Omega,200)
    tlabel = '$\Omega t$'
    tplot = trun*Omega
else :
    trun = np.linspace(0,5*T1,200)
    tlabel = '$t/T_1$'
    tplot = trun/T1


Ysoln = np.zeros( (len(trun), 2) ) + 1j*np.zeros( (len(trun), 2) )
Ysoln[0,:] = [sigma_init,w_init]

# %%
# Define RHS function =========================================================

def RHS(t,Y):
    # Unpack Y
    # print('In RHS, Y = ')
    # print(Y)
    sigma = Y[0]
    w = Y[1]
    
    Ydot = np.zeros( (2, ) ) + 1j*np.zeros( (2, ) )
    
    Ydot[0] = (1j*Delta-(1/T2))*sigma - (1j/2)*Omega*w
    Ydot[1] = -(1/T1)*(w-weq) + 2*np.imag(np.conj(Omega)*sigma)
    
    return Ydot


r = ode(RHS) 
r.set_integrator('zvode', method='bdf')
r.set_initial_value(y = [0,-1], t=0.0)


for i in range(1,len(trun)):
    Ysoln[i,:] = r.integrate(trun[i])
    
sigma_soln = Ysoln[:,0]
w_soln = Ysoln[:,1]

rho11 = (1-w_soln)/2
rho22 = (1+w_soln)/2

fig, (ax1,ax2) = plt.subplots(2, sharex = True)
ax1.set_title('$\Omega^2 T_1 T_2 =$ ' + str(round(T1*T2*Omega**2,2)) \
              +', $\Delta/\Omega =$ ' + str(round(Delta/Omega,2)))
ax1.plot(tplot,rho11,label = '$\\rho_{11}(t)$')
ax1.plot(tplot,rho22,label = '$\\rho_{22}(t)$')

ax1.set_xlabel(tlabel)
ax1.set_ylabel('Populations')

ax1.set_ylim(-0.1,1.1)
ax1.legend()
ax1.grid()

ax2.plot(tplot,np.real(sigma_soln),label = 'Re$[\sigma(t)]$')
ax2.plot(tplot,np.imag(sigma_soln),label = 'Im$[\sigma(t)]$')
ax2.set_ylim(-1.1,1.1)
ax2.legend()
ax2.grid()
ax2.set_xlabel(tlabel)
ax2.set_ylabel('Coherence')

# %%
# Plot real and imaginary parts of the susceptability ==================

# Amp = 1 #Factor that would appear out in front in real units

# FWHM = 2/T2*np.sqrt((1+T1*T2*Omega**2)) #FWHM including power broadening

# DeltaPlot = np.linspace(-5*FWHM,5*FWHM,90000)

# Chi = T2*(1j-DeltaPlot*T2)*1/(1 + (DeltaPlot*T2)**2 + T1*T2*Omega**2)
# fig, (ax3,ax4) = plt.subplots(2, sharex = True)
# ax3.set_title( '$T_2 = $ ' + str(round(T2,2)) + ', $\Omega^2 T_1 T_2 =$ ' + str(round(T1*T2*Omega**2,2)) )
# ax3.plot(DeltaPlot,np.imag(Chi))
# ax3.grid()
# ax3.set_xlabel('$\Delta = \omega - \omega_{21}$')
# ax3.set_ylabel("$\chi''$ [arb. units]")

# ax4.plot(DeltaPlot,np.real(Chi))
# ax4.grid()
# ax4.set_xlabel('$\Delta = \omega - \omega_{21}$')
# ax4.set_ylabel("$\chi'$ [arb. units]")

# # %%
# # Plot real and imaginary parts of the susceptability for a range of T2 ====

# T2inv_plot = np.linspace(T1inv,1/T2,3)
# T2plot = 1/T2inv_plot

# fig, (ax5,ax6) = plt.subplots(2, sharex = True)
# for T2i in T2plot:
#     Chi = T2i*(1j-DeltaPlot*T2i)*1/(1 + (DeltaPlot*T2i)**2) # + T1*T2i*Omega**2)   
#     ax5.plot(DeltaPlot,np.imag(Chi), label = '$1/T_2$ = ' + str(round(1/T2i,2)))
#     ax5.grid()
#     ax5.set_xlabel('$\Delta = \omega - \omega_{21}$')
#     ax5.set_ylabel("$\chi''$ [arb. units]")
#     ax5.legend(loc = 'upper right')
    
#     ax6.plot(DeltaPlot,np.real(Chi),label = '$1/T_2 = $' + str(1/T2i))
#     ax6.grid()
#     ax6.set_xlabel('$\Delta = \omega - \omega_{21}$')
#     ax6.set_ylabel("$\chi'$ [arb. units]")
    
plt.show()




