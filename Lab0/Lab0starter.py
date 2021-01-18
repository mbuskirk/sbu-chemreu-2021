"""
Lab0starter.py

This is starter code for lab0 of CHE 525. Lab0 is mostly just about getting
started with running code on Seawulf, so most work here is done for the students already.

Origin: CHE 525, S21 problem development
Author: Tom Allison
"""

#%%
# Import modules ===============================================================
import psi4
import numpy as np  
import matplotlib.pyplot as plt
import h5py

#%%
# Set up Psi4 =================================================================

psi4.core.clean()
psi4.core.clean_options()
psi4.set_memory('4000 MB')  # Can make this much larger on Seawulf, each compute node has more than 100 GB RAM
psi4.set_num_threads(8)    # Can make this much larger on Seawulf, each compute node can support 24 threads. 

#%%
# Set up molecule and file i/o ================================================

molecule_name = 'ethene'
method = 'scf/3-21G'
outfilename = molecule_name + '-scf'
psi4.core.set_output_file(outfilename + '.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

C2H4 = psi4.geometry(""" 
 C     0.000000     0.000000     0.000000
 H     0.000000     0.000000     1.070000
 C     1.156144     0.000000    -0.667500
 H    -0.943102     0.000000    -0.544500
 H     1.156144     0.000000    -1.756500
 H     2.099246     0.000000    -0.123000
   """)  # Z-matrix from molden in XYZ format

#%%
# Calculate energy and optimize geometry and output results ============
E0 = psi4.energy(method)                                                # calcualte energy at not optimized position
E_opt, wfn_opt = psi4.optimize(method, molecule=C2H4,return_wfn = True) # calculate energy and `wave function' at optimized position

psi4.driver.molden(wfn_opt, outfilename + '.molden') #save wave function to file in molden format

print()
print('============ End Psi4 Output ====================')
print('=================================================')
print('=================================================')
print()
print('Hello World')
print()
print(molecule_name + ' energy before optimization is ' + str(E0))
print(molecule_name + ' energy after optimization is ' + str(E_opt))
print('Calculation complete look at the output files to see more results')
print()




