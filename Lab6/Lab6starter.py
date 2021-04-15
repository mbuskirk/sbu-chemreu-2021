"""
Lab6starter.py

This code uses Psi4 to compute the change in energy
for the reaction 2 O_3 -> 3 O_2 at various levels of 
theory. 

Origin CHE 525, S21 Problem Development
Author: Ben Levine
"""
#%%
# Import modules =======================================================
import psi4
import numpy as np  

#%%
# Set up Psi4 ==========================================================

psi4.core.clean()
psi4.core.clean_options()
psi4.set_memory('4000 MB')  # Can make this much larger on Seawulf, each compute node has more than 100 GB RAM
psi4.set_num_threads(4)    # Can make this much larger on Seawulf, each compute node can support 28 threads.
                            # But it doesn't help much for small molecules...
psi4.core.set_output_file('Lab6starter.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

#%%
# Define QC method to use ==============================================

# we will optimize with this method
method1 = 'b3lyp/cc-pvdz'

# we will then calculate single-point energies at the optimized structure with this set of functionals
functionals = {'ccsd(t)','pbe','b3lyp'}

# and this basis
basis = 'cc-pvdz'

#%%
# Set up Z-matrix string and perform initial optimization ==============
O3_string= """
    o
 o   1 oo2
 o   2 oo3 1 ooo3     

oo2=        1.2
oo3=        1.3
ooo3=       120.0"""

# Define the molecules.  
O3 = psi4.geometry(O3_string) 
    
# Run a geometry optimization
EO3, wfnO3 = psi4.optimize(method1, molecule = O3, return_wfn = True)

# Here we create an empty Python dictionary to store the computed energies 
EnergiesO3 = {}

# loop over all functionals and compute the energy
for funcl in functionals:
    # concatinate the functional with the basis before we pass it to psi4
    method2 = funcl + '/' + basis

    # compute energy
    E_funcl, wfn_funcl = psi4.energy(method2, molecule = O3, return_wfn = True)
 
    # print molden output
    psi4.molden(wfn_funcl,"Lab6starter-O3."+funcl+".molden")

    # store energy in python dictionary
    EnergiesO3[funcl] = E_funcl

#%%
# Set up Z-matrix string and perform initial optimization ==============
# The first line here sets the charge to zero and the spin multiplicity to 3 (triplet)
# This is necessary in this case, because the ground state of O2 is a spin triplet, not
# a spin singlet, which is the default.
O2_string= """
0 3
    o
 o   1 oo2

oo2=        1.2"""

# Because O2 is open shell, we will use an unrestricted reference for all calculations
psi4.set_options({'reference':'uhf'})

# Define the molecules.  Here the format function replaces
# placeholder, {0}, with the value of theta
O2 = psi4.geometry(O2_string)

# Run a calculation of the energy without optimizing the geometry
EO2, wfnO2 = psi4.optimize(method1, molecule = O2, return_wfn = True)

# Here we create an empty Python dictionary to store the computed energies
EnergiesO2 = {}

# loop over all functionals and compute the energy
for funcl in functionals:
    # concatinate the functional with the basis before we pass it to psi4
    method2 = funcl + '/' + basis

    # compute energy
    E_funcl, wfn_funcl = psi4.energy(method2, molecule = O2, return_wfn = True)

    # print molden output
    psi4.molden(wfn_funcl,"Lab6starter-O2."+funcl+".molden")

    # store energy in python dictionary
    EnergiesO2[funcl] = E_funcl

# Now we will compute the change in energy (Delta U) for the reaction

# First initial another empty Python dictionary
DeltaU = {}

# loop over all functional and computer energy of reaction:
# 2 O3 -> 3 O2

for funcl in functionals:
    # compute energy
    DeltaU[funcl] = 3.0 * EnergiesO2[funcl] - 2.0 * EnergiesO3[funcl]

# Finally, we will print out a CSV file with all Delta U values for each functional
with open('Lab6starter.csv', 'w') as f:
    for funcl in functionals:
        f.write("%s,%e\n"%(funcl,DeltaU[funcl]))




