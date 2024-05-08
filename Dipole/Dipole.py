# -*- coding: utf-8 -*-
"""
Created on Wed May  8 11:27:11 2024

@author: Alex
"""

from ase.build import molecule
import numpy as np
import matplotlib.pyplot as plt

H2O = molecule("H2O")

xaxis = 1
yaxis = 2

plt.scatter(H2O.positions[:,xaxis], H2O.positions[:,yaxis])
print(H2O.get_chemical_symbols())
charges = np.array([-1, 0.5, 0.5])

dipole = charges*H2O.positions
print(dipole)

dipole = dipole.sum(axis=0)
print(dipole)

# Find the direction to draw the arrow

if np.sum(np.abs(H2O.positions[0] - dipole)) > np.sum(np.abs(H2O.positions[0])):
    plt.arrow(dipole[xaxis], dipole[yaxis], -dipole[xaxis], -dipole[yaxis], width= 0.02, label="dipole")
else:
    plt.arrow(0, 0, dipole[xaxis], dipole[yaxis], width= 0.02, label="dipole")

plt.title("The convention in chemistry is that the arrow representing the dipole moment goes from positive to negative.")

plt.show()