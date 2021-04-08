#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SincDeltaFcn.py

This script shows how a sinc^2 becomes a delta-function in the long-time limit.

Origin: CHE 525, S21, Light-matter interaction lectures
Author: Tom Allison
"""

# %%
# Preliminaries ===============================================================

import numpy as np
import matplotlib.pyplot as plt

# %%
# Define functions ============================================================

def cn(delta,t): 
    return 1*(1-np.exp(1j*delta*t))/(delta)

# %%
# Make Plots ==================================================================

delta = np.linspace(-5,5,300)

times = [5,30]

fig = plt.figure()
ax = plt.axes()
for t in times:
    ax.plot(delta,cn(delta,t)*np.conj(cn(delta,t)), label = 't = ' + str(t))

#ax.set_ylim(-1,600)    
ax.set_xlabel('$\Delta = \omega_{n0} - \omega$ [arb. units]')
ax.legend()
ax.grid()
  


