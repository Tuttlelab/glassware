# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 23:21:57 2023

@author: Alex
"""
import pandas
import numpy as np

def readin(fpath):
    f = open(fpath, 'r', encoding="utf8", errors='ignore')
    content = f.read()
    return content

class topol_parser:
    def write_topol(self, filename, bonds, bond_data, angles, LJ):
        with open(filename, 'w') as otop:
            otop.write("; Topology written by topol_parse.write_topol\n")
            otop.write("; Element MW Z\n")
            for atom, mass, charge in zip(self.Atoms, self.masses, self.charges):
                otop.write(f"{atom}\t{mass[0]}\t{charge[0]}\n")
            otop.write(self.spacer+"\n")
            otop.write("; Bonds: i, j, rmin (nm), K\n")
            for bond in bonds:
                bond = list(bond)
                bond[0] = int(bond[0])
                bond[1] = int(bond[1])
                otop.write("\t".join([str(x) for x in bond]))
                #otop.write("\t")
                
                otop.write("\n")
            otop.write(self.spacer+"\n")
            otop.write("; Angles: i, j, k,  theta_min (°), K\n")
            for angle in angles:
                angle = list(angle)
                angle[0] = int(angle[0])
                angle[1] = int(angle[1])
                angle[2] = int(angle[2])
                angle[3] = np.rad2deg(angle[3]) 
                otop.write("\t".join([str(x) for x in angle]))
                otop.write("\n")
            otop.write(self.spacer+"\n")
            otop.write("; Dihedrals: i, j, k,j,  theta_min (°), K\n")
            otop.write(self.spacer+"\n")
            otop.write("; LJ terms: i, j, sigma, eps ; The LJ potential has its minimum at a distance of r = rmin = 2^(1/6)sigma, where the PE has the value V = -eps\n")
            otop.write("; This one we'll do NAMD-style\n")
            otop.write(";i, j, rmin, Emin \n")
            for index in LJ.index:
                #otop.write(f"{index}\t{col}\t{self.LJ_Rmin.at[index, col]}\t{self.LJ_Emin.at[index, col]} \n")
                otop.write("\t".join([str(x) for x in LJ.loc[index]]))
                otop.write("\n")
                
    def parse_topol(self):
        topol = readin(self.topol_file)
        topol = topol.split(self.spacer)
        #print(topol)
        if len(topol) > 0:
            elements = topol[0].split("\n")
            for line in elements:
                if len(line.replace(" ", "")) < 1:
                    continue
                if line.replace(" ", "")[0] == ";":
                    continue
                line = line.split()
                atom = line[0]
                self.Atoms.append(atom)
                self.masses.append(float(line[1]))
                self.charges.append(float(line[2]))
        self.masses = np.array(self.masses).reshape(-1, 1)
        self.charges = np.array(self.charges).reshape(-1, 1)
        
        if len(topol) > 1:
            bonds = topol[1].split("\n")
            for line in bonds:
                if len(line.replace(" ", "")) < 1:
                    continue
                if line.replace(" ", "")[0] == ";":
                    continue
                line = line.split()
                #self.bonds.append([int(line[0]), int(line[1]), float(line[2]), float(line[3])])
                self.bonds.append([int(line[0]), int(line[1])])
                i, j, rmin, K = int(line[0]), int(line[1]), float(line[2]), float(line[3])
                index = [self.Atoms[i], self.Atoms[j]]
                index.sort()
                index = f"{index[0]}-{index[1]}"
                self.bond_data.at[index, "i"] = self.Atoms[i]
                self.bond_data.at[index, "j"] = self.Atoms[j]
                self.bond_data.at[index, "Rmin"] = rmin
                self.bond_data.at[index, "K"] = K
            self.bonds = np.array(self.bonds)
        
        if len(topol) > 2:
            angles = topol[2].split("\n")
            for line in angles:
                if len(line.replace(" ", "")) < 1:
                    continue
                if line.replace(" ", "")[0] == ";":
                    continue
                line = line.split()
                self.angles.append([int(line[0]), int(line[1]), int(line[2]), np.deg2rad(float(line[3])), float(line[4])])
        self.angles = np.array(self.angles)
        
        if len(topol) > 3:
            dihedrals = topol[3].split("\n")
            for line in dihedrals:
                if len(line.replace(" ", "")) < 1:
                    continue
                if line.replace(" ", "")[0] == ";":
                    continue
                line = line.split()
                self.dihedrals.append([int(line[0]), int(line[1]), int(line[2]), int(line[3]), np.deg2rad(float(line[4])), float(line[5])])
        self.dihedrals = np.array(self.dihedrals)
        
        if len(topol) > 4:
            LJs = topol[4].split("\n")
            for line in LJs:
                if len(line.replace(" ", "")) < 1:
                    continue
                if line.replace(" ", "")[0] == ";":
                    continue
                line = line.split()
                i = line[0]
                j = line[1]
                self.LJ.at[f"{i}-{j}", "i"] = line[0]
                self.LJ.at[f"{i}-{j}", "j"] = line[1]
                self.LJ.at[f"{i}-{j}", "Rmin"] = float(line[2])
                self.LJ.at[f"{i}-{j}", "Emin"] = float(line[3])
                
    def __init__(self, topol):
        self.topol_file = topol
        self.Atoms = []
        self.masses = []
        self.charges = []
        self.bonds = []
        self.bond_data = pandas.DataFrame()
        self.angles = []
        self.dihedrals = []
        self.LJ = pandas.DataFrame()
        self.spacer = "-"*50
        