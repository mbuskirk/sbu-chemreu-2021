"""
Lab4starter-CAS.py

This code uses Psi4 to perform a PES scan of the
bond dissociation of H2 at the CASSCF level with a
two electron/two orbital active space. 

Origin CHE 525, S21 Problem Development
Author: Tom Allison, Ben Levine
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
psi4.set_num_threads(1)    # Can make this much larger on Seawulf, each compute node can support 28 threads.
                            # But it doesn't help much for small molecules...
psi4.core.set_output_file('Lab4starter-CAS.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

#%%
# Define QC method to use and angles to scan ===========================

method = 'casscf/cc-pvdz'

#%%
# Set up Z-matrix string and perform initial optimization ==============
# This geometry was taken from the previous geometry optimization.
# Specifically, we opened the previous molden file in molden, opened
# molden's Z-matrix editor, and read the values.  The dihedral angle
# is replaced by the placeholder {0}.  Below, you will see how we 
# can replace this placeholder with a series of bond distances to scan the 
# PES.
H2_string= """
    h
 h   1 hh2     

hh2=        {0}"""

# Here is a list of bond distances that we want to compute, in Angstron
distances = [0.7,0.8,0.9,1.0,1.1,1.2,1.5,2.0,3.0,5.0,10.0,20.0]

# This creates an empty 2-D array (matrix) that will store the energies 
# that come from our calculations.  The argument (0,2) tells empty that 
# we would request a 2 x 0 array, meaning that each row will have 2 
# elements, but initially there are zero rows.  We will add rows below.
Er = np.empty((0,2))

# This for loop will repeat for every element of the array angles, 
# defined on the previous line.  For each iteration, the variable theta
# will be set to a different angle (first 180, then 165, etc...)
for r in distances:

    # Define the molecules.  Here the format function replaces 
    # placeholder, {0}, with the value of theta
    H2 = psi4.geometry(H2_string.format(r)) 

    H2.reset_point_group('c1') # turn symmetry off!

    psi4.set_options({'fail_on_maxiter':'false'})
    
    # Here we are setting up the active space.  restrict_docc is the number of 
    # orbitals that will be restricted to be doubly occupied.  active is the
    # number of active orbitals.  In this case, we want 2 active orbitals,
    # specifically the bonding and antibonding combinations of the 1s 
    # atomic orbitals.  There are no core orbitals below, so irestrict_docc 
    # is 0.  
    psi4.set_options({'restricted_docc':[0]})
    psi4.set_options({'active':[2]})

    # Run a CASSCF calculation.  
    Ei, wfni = psi4.energy(method, molecule = H2, return_wfn = True)

    # This adds the bond distance (r) and energy (Ei) to the array (Er)
    # that we created above.  Each row will contain 2 values: r and Ei
    Er = np.append(Er, [[r,Ei]], axis=0)

    # The object wfn0 now contains the optimized geometry and 
    # associated wave function.  You can output this information 
    # to a molden file as follows.  Here str converts r from a 
    # number to a string, and the + sign concatenates several
    # strings into one long sting.  For example, if r is 1.1, the 
    # filename will be something like "Lab4starter-CAS.1.1.molden").  This 
    # allows us to save a new molden file for each iteration of the loop.
    psi4.molden(wfni,"Lab4starter-CAS."+str(r)+".molden")

# Here we print the results array (Er) to a CSV file that can be
# read later by a python plotting program or the spreadsheet of your
# choice.  CSV files are also human readable.
np.savetxt("Lab4starter-CAS.csv", Er, delimiter=",",header='H-H Distance,E',comments='')
