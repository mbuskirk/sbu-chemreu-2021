#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 07:30:29 2021

@author: tka3
"""

#%%
# Import modules =======================================================

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#%%
# Read in data =========================================================
df = pd.read_csv('OHpotentials.csv')
Nroots = len(df.keys()) - 3
Labels = df.keys()[2:]

# %%
# Unpack data into numpy arrays and subtract offset ====================
R = np.asarray(df['R [Angstrom]'])
Eg = np.asarray(df['Eg'])

# subtract dissociation limit for ground state
Eg0au = Eg[-1]
Eg = 27.211399*(Eg-Eg0au)


fig = plt.figure()
ax = plt.axes()
ax.plot(R,Eg, marker = 'o', label = 'Ground State')
for j in range(Nroots):
    Ejvec = np.asarray(df[df.keys()[3+j]]) #unpack excited states
    Ejvec = 27.211399*(Ejvec-Eg0au)
    ax.plot(R,Ejvec, marker = 'o', label = Labels[j+1])
ax.set_xlabel('R [Angstrom]')
ax.set_ylabel('U(R) [eV]')
ax.grid()
ax.legend()
plt.show()
