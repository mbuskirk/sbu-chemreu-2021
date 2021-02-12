"""
Lab2_makeplots.py

This script reads in data from a .csv file and plots the results.

Origin: CHE525 problem development, Spring 2021
Author: Tom Allison
"""

#%%
# Import modules =======================================================

#import psi4
import numpy as np  
import matplotlib.pyplot as plt
import pandas as pd

#%%
# Read in data from file and unpack int numpy arrays ====================
df = pd.read_csv('Lab2starter_2_RigidRotation.csv')   # read in data from the spreadsheet file into a DataFrame object.
angles = np.asarray(df['H-O-O-H Dihedral Angle']) # extract the column labeled 'H-O-O-H dihedral angle', convert to numpy array, and store in variable angles.
Eth_rigid = np.asarray(df['E_rigid'])         # extract the column labeled 'E_rigid', convert to numpy array, and store in variable E_rigid

df2 = pd.read_csv('Lab2starter_3_RelaxedRotation.csv') # now read data for relaxed rotation
Eth_relaxed = np.asarray(df2['E_relaxed'])     # extract the column labeled 'E_relaxed', convert to numpy array, and store in variable E_relaxed

#%%
# Subtract off ofsets and ocnvert to eV  ===================================
Eth_rigid = (Eth_rigid - Eth_rigid[0])*27.2114 # subtract off 180 deg. value and convert to eV
Eth_relaxed = (Eth_relaxed - Eth_relaxed[0])*27.2114 # subtract off 180 deg. value and convert to eV

#%%
# Plot data using Matplotlib ============================================

fig = plt.figure() # initialize new figure.
ax = plt.axes()    # initialize new axes

ax.plot(angles,Eth_rigid, marker = 'o', label = 'rigid')  
ax.plot(angles,Eth_relaxed, marker = 'o', label = 'relaxed')
ax.set_ylabel('Energy relative to planar geometry [eV]')
ax.set_xlabel('Dihedral angle [deg.]')
ax.set_title('HOOH dihedral angle')
ax.grid()   # turn grid on
ax.legend() # show legend
plt.show()  # show plot. Need this line if not running from IPython (e.g. Spyder) or Jupyter. 
