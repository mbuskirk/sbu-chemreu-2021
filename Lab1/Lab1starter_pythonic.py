"""
Lab1starter_pythonic.py

This code executes all the functions needed for Lab1 in a slick, automated
way using Python data structures like dictionaries and lists. 

If you are confused by this, look at the other two starter codes.

Origin CHE 525, S21 Problem Development
Author: Tom Allison
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
psi4.core.set_output_file('Lab1starter_pythonics.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

#%%
# Define methods to use and .xyz files to read =========================
methods = ['scf/3-21G', 'scf/6-311G*'] # List of methods to use. Last element will be used for optimization!
mol_names = ['OCS_linear' , 'OCS_bent']

#%%
# Initialize Python dictionaries to store results ======================
    
molecules = {} # initialize empty dictionary for storing molecule objects
E0 = {}        # initialize empty dictionary for storing un-optimized energies
Eopt = {}      # initialize empty dictionary for storing optimized energies
wfn_opt = {}   # initialize empty dictionary for storing optimized wave functions

#%%
# Load molecules from using either .xyz files or cut/paste Z-matrices ==
# ======================================================================
# psi4.geometry() creates a molecule object, in both cases stored as an
# element in the dictionary defined above
# Here we will not user a for loop to illustrated doing it different ways,
# but obviously this could be further automated, which would be sensible
# if running many molecules...

# example of loading geometry from Z-matrix cut/paste
molecules['OCS_bent'] = psi4.geometry("""
 o
 c    1 co2     
 s    2 sc3         1 sco3      
 
co2=        1.380000
sc3=        1.780000
sco3=       120.000
""")
molecules['OCS_bent'].reset_point_group('c1') # turn symmetry off!
E0['OCS_bent'] = {} # initialize an empty dictionary to store the E0 results calculated with the different methods for each molecule.
Eopt['OCS_bent'] = 0 # initialize optimized energy to scalar 0.

#example of loading geometry from .xyz file
f = open('OCS_linear.xyz') #use Python's 'open' to open xyz file.
molecules['OCS_linear'] = psi4.geometry(f.read())
molecules['OCS_linear'].reset_point_group('c1') # turn symmetry off!
E0['OCS_linear'] = {} # initialize an empty dictionary to store the E0 results
Eopt['OCS_linear'] = 0 # initialize optimized energy to scalar 0

#%%
# Run calculations. Illustrate using a for loop ========================

for mol_name in mol_names:
    print('Working on ' + mol_name)    
    for meth_j in methods:                 #loop over methods
        E0[mol_name][meth_j] = psi4.energy(meth_j, molecule = molecules[mol_name])
    
    #use if statements to control the special cases of water_170 and methanol_eclipsed
    if mol_name == 'OCS_bent':
        psi4.set_module_options('optking',{'frozen_bend': '1 2 3'})
        
    Eopt[mol_name], wfn_opt[mol_name] = psi4.optimize(methods[-1], molecule = molecules[mol_name], return_wfn = True)
    psi4.driver.molden(wfn_opt[mol_name], mol_name + '.molden')
    
#%%
# Print energy results to screeen ======================================

print('Energies at initial positions with different methods:')
print()
print(E0)
print()
print(methods[-1] + ' energies after optimization:')
print(Eopt)

#%%
# Save optimized results to a spreadsheet using Pandas
df = pd.DataFrame(data = Eopt, index = [0])
df.to_csv('Lab1starter_Eopt.csv')  # save dataframe to .csv spreadsheet file.

#%%
# Save initial results to a spreadsheet using Pandas
df = pd.DataFrame(E0)
df.to_csv('Lab1starter_E0.csv')  # save dataframe to .csv spreadsheet file.
 
# Or of course you can extract the results from the raw .dat file
# generated by Psi4.
