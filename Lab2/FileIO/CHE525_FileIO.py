"""
CHE525_FileIO.py

This script is for exploring differen fileIO options for CHE52. Uses
an example of saving/retrieving three arrays.

Origin: CHE525_S21 Problem development
Author: Tom Allison
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py

print()
print('h5py version is ' + h5py.__version__)
print()

#Make some arrays
X = np.asarray([1, 2, 3, 4, 5])
Y = np.asarray([0, 12, 50, 80, 100])
Z = np.asarray([0, 11, 13, 8,  5])

#Option 1 put arrays in giant matrix and dump to csv with np.savetxt 
A = np.asarray([X, Y, Z]) 
np.savetxt('np.csv', A.T, fmt='%.2f', delimiter=',', header=" X, Y, Z")
B = np.loadtxt('np.csv', delimiter = ',')
print(B)
fig = plt.figure()
ax = plt.axes()
ax.plot(B[:,0],B[:,1], marker = 'o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('File IO w/ numpy')

#Option 2, use Pandas
d = {'X': X, 'Y': Y, 'Z': Z}
df_in = pd.DataFrame(data = d)
df_in.to_csv('pd.csv')
df_out = pd.read_csv('pd.csv')

print(df_out)
print(np.asarray(df_out['X']))
fig = plt.figure()
ax = plt.axes()
ax.plot(df_out['X'],df_out['Y'], marker = 'o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('File IO w/ pandas')

#Option 3, use h5py

#create file
hf_in = h5py.File('h5test.h5','w')
#dt = h5py.string_dtype(encoding='utf-8') #Doesn't work with h5py version < 3 it seems.


#create main data sets
hf_in.create_dataset(name = 'X', data = X)
hf_in.create_dataset(name = 'Y', data = Y) 
hf_in.create_dataset(name = 'Z',  data = Z)
hf_in.create_dataset(name = 'matrix', data = A.T)

# create metadata group
metagrp = hf_in.create_group('run info')
metagrp.create_dataset(name = 'course number', data = 'che525') #Can write string to see with HDFviewer, but having trouble reading with h5py.
metagrp.create_dataset(name = 'matrix shape', data = np.shape(A.T))

hf_in.close()

hf_out = h5py.File('h5test.h5', 'r')
X5 = hf_out['X'][:]
Y5 = hf_out['Y'][:]
Z5 = hf_out['Z'][:]

fig = plt.figure()
ax = plt.axes()
ax.plot(X5,Y5, marker = 'o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('File IO forw/ h5py')

hf_out.close()

plt.show()



