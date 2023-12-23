# -*- coding: utf-8 -*-
#https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


mols = [read("QN.xyz")]
mols[0].positions -= mols[0].positions.min(axis=0)

#axis = mols[0].positions.mean(axis=0)
axis = [1,1,1] # Rotate around Z axis

X, Y = [], []
for theta in np.linspace(0., 6.3, 100):
    X.append(theta)
    Y.append(np.cos(theta / 2.0))
    mols.append(mols[0].copy())
    for i in range(mols[-1].positions.shape[0]):
        pos = mols[-1].positions[i]
        
        #print(np.dot(rotation_matrix(axis, theta), pos)) # check point gets rotated individually
        
        mols[-1].positions[i] = np.dot(rotation_matrix(axis, theta), pos)
        
for i, mol in enumerate(mols):
    mol.write("Out.xyz", append = (i!=0))


plt.plot(X, Y)

print("Theta value for 1 full rotation:", X[np.argmin(Y)])