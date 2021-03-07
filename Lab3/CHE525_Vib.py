#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CHE525_Vib.py

This package contains functions for easily extracting a molecule's dipole moment
and calculating the derivative of the dipole moment, and resulting "IR intensities"
for a molecule.

Origin: Stony Brook University CHE525 Problem development, S21
Author: Tom Allison
"""

import psi4
import numpy as np
import qcelemental as qcel
#from qcelemental import Datum

def get_dipole(wfn_in):
    """This function returns a 3-component numpy array for dipole moment extracted 
       from  a psi4 wave function object. It doesn't calculate anything. It just 
       pulls the data out of the wfn in a more convenient format."""
    psi4.oeprop(wfn_in,'DIPOLE')
    D = np.zeros( (3,) ) # initialize empty numpy array to store dipole vector components
    D[0] = wfn_in.variable(' DIPOLE X')
    D[1] = wfn_in.variable(' DIPOLE Y')
    D[2] = wfn_in.variable(' DIPOLE Z')
    D = (0.3934303)*D # convert to atomic units e a_0
    return D

def dipder(method, mol, eps = 0.01):
    """Dipole Derivatives
       
       This function performs vibrational analysis for a psi4 molecule object. 
       It fills a gap in standard Psi4, providing dipole derivatives and IR
       intensities. The derivatives are calculated using finite difference via
       displacing the molecule geometry a small amount in the direction of 
       the normal mode coordinates. 
       
       Required input arguments:
       method = QC method to use for the calculation, e.g. b3lyp/6-31G*
       mol = psi4 molecule object                                  
       
       Optional argument eps lets you set the (fractional) step size to take for
       the finite difference. Default value of 0.01 takes a step size of 1% of the normal mode.
       Making this too large reduces the accuracy of the finite difference as an approximation
       of the derivative. Making this too small can run into problems with numerical accuracy
       in the dipole calculation. 1% is a compromise that works reasonably well
       for many cases.
       
       This code returns a dictionary with the following keys:
       
       'dipder': Dipole derivatives for small change in normal mode coordinate. In units of dipole [au] * length [au]. Similar to what psi4 usually calculates for RHF
       'om': Vibrational frequencies in cm^{1}. Excludes rotations and translations.
       'IR_Intensity': Infrared transition intensities in km/mol
       
       """
    uconv_kmmol = (qcel.constants.get("Avogadro constant") * np.pi * 1.e-3 * qcel.constants.get("electron mass in u") *\
                   qcel.constants.get("fine-structure constant")**2 * qcel.constants.get("atomic unit of length") / 3)   # code taken from https://github.com/psi4/psi4/blob/master/psi4/driver/qcdb/vib.py
    
    XYZ0 = mol.geometry() #psi4 matrix object, initial geometry
    XYZ0_str = mol.save_string_xyz()
    mol = psi4.geometry(XYZ0_str) # redefine molecule in terms of XYZ matrix since having problems when molecule defined in terms of Z-matrix ?!
    XYZ1 = mol.geometry() # initialize displaced geometry
    XYZ0_np = XYZ0.to_array() #extract XYZ matrix  in numpy format

    #initialize dictionary with empty lists
    results = {}
    results['dipder'] = [] #dipole derivatives in a.u.
    results['om'] = []     #frequencies in cm^{-1}
    results['IR_Intensity'] = [] # IR intensities in km/mol
    
    E, wfn0 = psi4.frequency(method, molecule = mol, return_wfn = True)
    mol.fix_orientation(True) # turn reorientation off! Without this you will get nonsense if psi4 reorients the molecule between dipole calculations!
    D0 = get_dipole(wfn0)
    print('Initial dipole:')
    print(D0)
    TRV = wfn0.frequency_analysis['TRV'].data
    w = wfn0.frequency_analysis['w'].data
    om = wfn0.frequency_analysis['omega'].data
    
    for j in np.arange(len(om)):
        if TRV[j] == 'V':
            wn = w[:,j] #extract normal mode coordinates for normal mode n.
            #print(wn)
            results['om'].append(om[j])
            print('Working on mode ' + str(j) + ' with frequency ' + str(om[j]))
            disp_mat = np.zeros(np.shape(XYZ0))
            for i in np.arange(np.shape(XYZ0)[0]):                 
                disp_mat[i,:] = eps*wn[3*i:(3*i+3)] 
            XYZ1_np = XYZ0_np + disp_mat
            XYZ1 = XYZ1.from_array(XYZ1_np) 
            mol.set_geometry(XYZ1) #update geometry
            E, wfn1 = psi4.energy(method, molecule = mol, return_wfn = True)
            D1 = get_dipole(wfn1)
            print('D1 = ')
            print(D1)
            dDdw = (D1-D0)/(eps) # calculate dmu/deps. Note that this is not dD/dX or dD/dW, but instead (nearly) equivalent to the dot product of psi4's normal dipder with w.
            results['dipder'].append(dDdw)
            results['IR_Intensity'].append(uconv_kmmol*np.dot(dDdw,dDdw)) #Under construction! Need to get prefactor right!
    
    return results