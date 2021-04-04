"""
Lab5starter.py

This script is the starter code for Lab5. It performs excited-state 
calcuations on the OH radical.

Origin: CHE525 Lab Development, Spring 2021
Author: Tom Allison
"""

#%%
# Import modules =======================================================
import psi4
import numpy as np  
import pandas as pd

#%%
# Set up Psi4 ==========================================================

psi4.core.clean()
psi4.core.clean_options()
psi4.set_memory('4000 MB')  # Can make this much larger on Seawulf, each compute node has more than 100 GB RAM
psi4.set_num_threads(4)    # Can make this much larger on Seawulf, each compute node can support 28 threads.
                            # But it doesn't help much for small molecules...
psi4.core.set_output_file('PotentialScan.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

#%%
# Define QC method to use and set up calculation =======================

method = 'eom-ccsd/aug-cc-pvdz'

#Define Z-matrix string with adjustable bond length
OHstring = """
O
H 1 {0}
"""
# Do initial molecule initialization in order to extract the irreps.
OH = psi4.geometry(OHstring.format(2)) 
irreps = OH.irrep_labels()

# Set up vector with number of roots to find per irreducible representation
# ------------[A1, A2, B1  B2] 
RootsVector = [ 2,  2,  2,  0]
Nroots = np.sum(RootsVector)

# Set up open-shell calculation options. If having trouble converging, you
# can try a uhf referncer or change the method for guessing the starting point
# for the reference scf calculation d 
psi4.set_options({'reference': 'rohf', 'freeze_core':'true', \
    'roots_per_irrep':RootsVector, 'guess':'read'}) #Only need A1

# Set up internuclear distances to scan.
Rvals = [0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.7,1.9,2.1,2.3,2.5,3.0]

# %%
# Do calculation! ======================================================

# Set up empty 2D array to store the excited state energies and empty 1D
# list to store the ground state energy and the state labels.
Ee = np.zeros( (len(Rvals),Nroots) )
Eg0 = []
Labels = []

i = 0 #loop counter
# loop over internuclear distance
for R in Rvals:
     n = 0
     OH = psi4.geometry(OHstring.format(R)) #set internuclear distance
     print('Working on R = ' + str(R))
     psi4.energy(method, molecule = OH) # perform calculation.
     Eg0.append(psi4.core.variable('CCSD TOTAL ENERGY'))  #append ground state energy to list.   
     print('Ground State Energy = ' + str(Eg0[-1]))
     #loop through excited states
     for n in range(Nroots):
        Labels.append('CC Root ' + str(n+1)) # Append state label, could get more descriptive here.
        Ee[i,n] = psi4.core.variable('CC ROOT %d TOTAL ENERGY' %(n+1)) #store excited state enegies to 2D array            
    
     i += 1

# %%
# Save Data for potential scan =========================================

#Build dictionary to save with Pandas
D = {'R [Angstrom]': Rvals,'Eg': Eg0}
for n in range(Nroots):
    D[Labels[n]] = Ee[:,n]
        
# Save data        
df = pd.DataFrame(data = D) 
df.to_csv('OHpotentials.csv')

# %%
# Calculate energies and dipoles  near equilibrium geometry of the =====
# lowest excited A1 (\Sigma^+) state. ==================================
print('Potential Scan Done, Calculating States and Transition Dipoles at Equilibrium...')
OH = psi4.geometry(OHstring.format(1.0)) # set geometry near A1 state equilibrium.

psi4.core.set_output_file('EquilibriumSpectrum.dat') # change output file name.
psi4.energy(method,molecule = OH)         # calculate energy (put results in new output file).

psi4.core.set_output_file('OscStrength.dat', False) # change output file name to store osc. strength data.
psi4.properties(method, molecule = OH, properties= ['oscillator_strength']) # calculate oscillator strengths.






