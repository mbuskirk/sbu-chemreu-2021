"""
Lab3starter_OptimizationAndFrequencies.py

This code uses Psi4 to perform a geometry optimization for H2O, and  
at the Hartree Fock level, and then compute the vibrational frequencies.

This is a starter code for Lab3.

Origin CHE 525, S21 Problem Development
Author: Tom Allison, Ben Levine
"""
#%%
# Import modules =======================================================
import psi4
import numpy as np
import CHE525_Vib as vib
import pandas as pd

#%%
# Set up Psi4 ==========================================================

psi4.core.clean()
psi4.core.clean_options()
psi4.set_memory('4000 MB')  # Can make this much larger on Seawulf, each compute node has more than 100 GB RAM
psi4.set_num_threads(4)    # Can make this much larger on Seawulf, each compute node can support 28 threads.
                            # But it doesn't help much for small molecules...
psi4.core.set_output_file('MyLab3.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

#%%
# Define QC method to use and angles to scan ===========================

methods = ['scf/6-31G*','mp2/6-31G*','b3lyp/6-31G*','ccsd(t)/6-31G*']
#%%
# Set up Z-matrix string and perform initial optimization ==============
CO_string= open('CO.xyz')

# Define molecule  
CO = psi4.geometry(CO_string.read()) 


for method in methods:
    print('\n\n\n\n\n\n')
    print('CO '+method)
    # Run the geometry optimization 
    E0, wfn0 = psi4.optimize(method, molecule = CO, return_wfn = True) # perform initial optimization
    
    # This option will tell Psi4 that, upon calculating the frequencies, it 
    # should output the a molden file that will include information so you can 
    # visualize the normal modes.
    psi4.set_options({"normal_modes_write": True})
    
    # This options provides Psi4 with the name for the file.  (Psi4 will append
    # some additional characters to the file name.) 
    psi4.set_options({"writer_file_label": "CO"})
    
    # Compute the virbational frequencies. Note that Psi4 will remember the 
    # last calculation that you did, and by default will start using that 
    # (optimized) geometry.   
    psi4.frequencies(method)
    
    # Compute IR intensities using using CHE525_vib package
    # It may be helpful to adjust the step size eps to make sure the results are
    # stable with respect to eps, especially if you are getting strange numbers!
    vib_results = vib.dipder(method,CO, eps = 0.01) 
    
    # Display results at terminal
    print(vib_results['om'])
    print(vib_results['IR_Intensity'])
    
    # Save data fo file with Pandas. Very easy when you already have a Python dictionary!
    df = pd.DataFrame(data = vib_results) # inialize DataFrame object using vib_results dictionary
    df.to_csv('CO'+method.split('/')[0]+'vib_results.csv') # Save .csv file with the data.


# Set up Z-matrix string and perform initial optimization ==============
CH3Cl_string= open('CH3Cl.xyz')

# Define molecule  
CH3Cl = psi4.geometry(CH3Cl_string.read()) 
for method in methods:
    print('\n\n\n\n\n\n')
    print('CH3Cl '+method)
    # Run the geometry optimization 
    E0, wfn0 = psi4.optimize(method, molecule = CH3Cl, return_wfn = True) # perform initial optimization
    
    # This option will tell Psi4 that, upon calculating the frequencies, it 
    # should output the a molden file that will include information so you can 
    # visualize the normal modes.
    psi4.set_options({"normal_modes_write": True})
    
    # This options provides Psi4 with the name for the file.  (Psi4 will append
    # some additional characters to the file name.) 
    psi4.set_options({"writer_file_label": "CH3Cl"})
    
    # Compute the virbational frequencies. Note that Psi4 will remember the 
    # last calculation that you did, and by default will start using that 
    # (optimized) geometry.   
    psi4.frequencies(method)
    
    # Compute IR intensities using using CHE525_vib package
    # It may be helpful to adjust the step size eps to make sure the results are
    # stable with respect to eps, especially if you are getting strange numbers!
    vib_results = vib.dipder(method,CH3Cl, eps = 0.01) 
    
    # Display results at terminal
    print(vib_results['om'])
    print(vib_results['IR_Intensity'])
    print(E0)
    print(wfn0)
    
    # Save data fo file with Pandas. Very easy when you already have a Python dictionary!
    df = pd.DataFrame(data = vib_results) # inialize DataFrame object using vib_results dictionary
    df.to_csv('CH3Cl'+method.split('/')[0]+'vib_results.csv') # Save .csv file with the data.

