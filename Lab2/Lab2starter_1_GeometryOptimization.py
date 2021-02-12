"""
Lab2starter_1_GeometryOptimization.py

This code uses Psi4 to perform a geometry optimization for 
hydrogen peroxide at the Hartree Fock level.

This is a starter code for Lab2.

Origin CHE 525, S21 Problem Development
Author: Tom Allison, Ben Levine
"""
#%%
# Import modules =======================================================
import psi4
import numpy as np  
import matplotlib.pyplot as plt
import pandas as pd

#%%
# Set up Psi4 ==========================================================

psi4.core.clean()
psi4.core.clean_options()
psi4.set_memory('4000 MB')  # Can make this much larger on Seawulf, each compute node has more than 100 GB RAM
psi4.set_num_threads(4)    # Can make this much larger on Seawulf, each compute node can support 28 threads.
                            # But it doesn't help much for small molecules...
psi4.core.set_output_file('Lab2starter_1_GeometryOptimization.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

#%%
# Define QC method to use and angles to scan ===========================

method = 'scf/6-31G'

#%%
# Set up Z-matrix string and perform initial optimization ==============
HOOH_string= """
    o
 o   1 oo2     
 h    1 ho3         2 hoo3      
 h    2 ho4         1 hoo4          3 dih4   
 
oo2=        1.4
ho3=        1.0
hoo3=       107.0
ho4=        1.0
hoo4=       104.0
dih4=       180.0"""

# Define molecule  
HOOH = psi4.geometry(HOOH_string) 

# Run the geometry optimization 
E0, wfn0 = psi4.optimize(method, molecule = HOOH, return_wfn = True) # perform initial optimization

# The object wfn0 now contains the optimized geometry and associated wave function
# You can output this information to a molden file as follows
psi4.molden(wfn0,"Lab2starter_1_GeometryOptimization.molden")


    
