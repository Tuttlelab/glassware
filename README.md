# Tuttlelab/glassware
## A repo for frequently used code snippets to help settle new undergrads and MSc students in our lab.

They can all be installed at once by running:
```bash
git clone https://github.com/Tuttlelab/glassware
cd glassware
pip install .
```

### hpctools

Parsing the archie-west HPC queue.

```python
import hpctools

hpc = hpctools.hpctools()
hpc.count_jobs()
> {"jobs":number of queued jobs, "cores": requesting n total cores}

hpc.IsRunning("job name", exact = False)
> True/False if any job with "job name" is in the queue
hpc.IsRunning("job name", exact = True)
> True/False if any job in the queue has the name "job name" 
```

### rmsd_tools

Parsing the archie-west HPC queue.

```python
import rmsd_tools

Test = np.random.random((10000, 18, 3))
RMSD = rmsd_tools()
RMSD.minimize_rotation_and_translation(Test[0], Test[i])
```

# peptideutils

Some tools for working with peptides in python

```python
import peptideutils

peptideutils.GenerateDatasetIndex(3)
> ["AAA", "AAC", ... "YYY"]

peptideutils.charge("AKRDF")
> 1
```



# Gromacs tools

### Parse SASA from an output .xvg file and return an AP score based on the first and last values
```python
sasa = np.genfromtxt("SASA.xvg", comments="@", skip_header = 14)
timeseries, sasa = sasa[:,0], sasa[:,1]
AP = sasa[0] / sasa[1]
```

## Parse and ITP file into several pandas DataFrames, does not currently consider dihedrals
```python
def parse_itp(fname):
    text = readin(fname)
    
    section = text.split("[ atoms ]")[1].split("[")[0]
    ATOMS = pandas.DataFrame(columns=["id", "type", "resnr", "residue", "atom", "cgnr", "charge"])
    i = 0
    for line in section.split("\n"):
        if len(line.split(";")[0].strip()) == 0:
            continue
        line = line.split(";")[0] #Only look at the uncommented part
        line=line.split()
        if len(line) < 6:
            continue
        ATOMS.loc[i] = line
        i+=1
    ATOMS["id"] = ATOMS["id"].astype(np.int64)
    ATOMS["resnr"] = ATOMS["resnr"].astype(np.int64)
    
    BONDS = pandas.DataFrame(columns=["i", "j", "func", "length", "fc", "comment"])
    if "[ bonds ]" in text:
        section = text.split("[ bonds ]")[1].split("[ ")[0]
        i = 0
        for line in section.split("\n"):
            if len(line.split(";")[0].strip()) == 0:
                continue
            line = line.split(";")[0] #Only look at the uncommented part
            line=line.split()
            if len(line) < 5:
                continue
            BONDS.loc[i] = line + [";"]
            i+=1
        BONDS["i"] = BONDS["i"].astype(np.int64)
        BONDS["j"] = BONDS["j"].astype(np.int64)
    
    CONSTRAINTS = pandas.DataFrame(columns=["i", "j", "func", "length"])
    if "[ constraints ]" in text:
        section = text.split("[ constraints ]")[1].split("[ ")[0]
        
        i = 0
        for line in section.split("\n"):
            if len(line.split(";")[0].strip()) == 0:
                continue
            line = line.split(";")[0] #Only look at the uncommented part
            line=line.split()
            if len(line) < 4:
                continue
            CONSTRAINTS.loc[i] = line 
            i+=1
        CONSTRAINTS["i"] = CONSTRAINTS["i"].astype(np.int64)
        CONSTRAINTS["j"] = CONSTRAINTS["j"].astype(np.int64)
    
    ANGLES = pandas.DataFrame(columns=["i", "j", "k", "func", "angle", "fc"])
    
    if "[ angles ]" in text:
        section = text.split("[ angles ]")[1].split("[ ")[0]
        i = 0
        for line in section.split("\n"):
            if len(line.split(";")[0].strip()) == 0:
                continue
            line = line.split(";")[0] #Only look at the uncommented part
            line=line.split()
            if len(line) < 5:
                continue
            ANGLES.loc[i] = line
            i+=1
        ANGLES["i"] = ANGLES["i"].astype(np.int64)
        ANGLES["j"] = ANGLES["j"].astype(np.int64)
        ANGLES["k"] = ANGLES["k"].astype(np.int64)
        
    exclusions = pandas.DataFrame(columns=[0,1,2,3])
    if "[ exclusions ]" in text:
        section = text.split("[ exclusions ]")[1].split("[ ")[0]
        i = 0
        for line in section.split("\n"):
            if len(line.split(";")[0].strip()) == 0:
                continue
            line = line.split(";")[0] #Only look at the uncommented part
            line=line.split()
            if len(line) < 2:
                continue
            for col in range(len(line)):
                exclusions.at[i,col] = int(line[col])
            i+=1

    
    return {"Atoms":ATOMS, "Bonds":BONDS, "Constraints": CONSTRAINTS, "Angles": ANGLES, "Exclusions": exclusions}
```



## Combine two itp files into 1 molecule
```python
def combine_itp(itp0, itp1):
    L0 = itp0["Atoms"].shape[0]
    itp1["Atoms"]["id"] += L0
    itp1["Atoms"].index += L0
    itp1["Atoms"]["resnr"] += itp0["Atoms"]["resnr"].max()
    itp1["Bonds"].index += itp0["Bonds"].shape[0]
    itp1["Bonds"]["i"] += L0
    itp1["Bonds"]["j"] += L0
    itp1["Constraints"].index += itp0["Constraints"].shape[0]
    itp1["Constraints"]["i"] += L0
    itp1["Constraints"]["j"] += L0
    itp1["Angles"].index += itp0["Angles"].shape[0]
    itp1["Angles"]["i"] += L0
    itp1["Angles"]["j"] += L0
    itp1["Angles"]["k"] += L0
    itp1["Exclusions"].index += itp0["Exclusions"].shape[0]
    itp1["Exclusions"][0] += L0
    itp1["Exclusions"][1] += L0
    itp1["Exclusions"][2] += L0
    itp1["Exclusions"][3] += L0
    return {"Atoms": pandas.concat((itp0["Atoms"], itp1["Atoms"])), 
            "Bonds": pandas.concat((itp0["Bonds"], itp1["Bonds"])), 
            "Constraints": pandas.concat((itp0["Constraints"], itp1["Constraints"])), 
            "Angles": pandas.concat((itp0["Angles"], itp1["Angles"])),
            "Exclusions": pandas.concat((itp0["Exclusions"], itp1["Exclusions"]))}

combined = combine_itp(itp0, itp1)
```


## Add a bond to an ITP
```python
def add_bond(itp, i, j, b0, fc):
    func=1
    itp["Bonds"].loc[itp["Bonds"].shape[0]] = [i, j, func, b0, fc, "; Additional Bond"]
    return itp
x = x[x["atom"] == "SC2"]
LYS = np.random.choice([i for i in x["id"] if i not in bonded_atoms])
x = combined["Atoms"]
x = x[x["type"] == "SP3"]
HDI = np.random.choice([i for i in x["id"] if i not in bonded_atoms])
add_bond(combined, LYS, HDI, 0.35, 5000)
```

## Add an angle to an ITP
```python
def add_angle(itp, i, j, k, b0, fc):
    func = 2
    print("ADD_ANGLE:", itp["Angles"])
    itp["Angles"].loc[itp["Angles"].shape[0]] = [int(i), int(j), int(k), int(func), b0, fc]
    for make_int in ["i", "j", "k", "func"]:
        itp["Angles"][make_int] = itp["Angles"][make_int].astype(np.int64)
    return itp
```

## Modify the bead type and/or charge of a bead in the itp
```python
def mod_type(itp, i, newtype=None, charge=None):
    index = itp["Atoms"][itp["Atoms"]["id"] == i].index[0]
    if newtype is not None:
        itp["Atoms"].at[index, "type"] = newtype
    if charge is not None:
        itp["Atoms"].at[index, "charge"] = charge
    return itp
```
## Write your dict of DataFrames to a new itp file:
```python
def write_itp(itp, fname):
    oitp = open(fname, 'w')
    oitp.write("""; written with write_itp
[ moleculetype ]
; molname        nrexcl
Combined         1

""")    
    for key in list(itp.keys()):
        oitp.write(f"[ {key.lower()} ]\n")
        
        for index in itp[key].index:
            if key.lower() != "exclusions":
                row = itp[key].loc[index]
                row = [str(x) for x in row]
            else:
                row = itp[key].loc[index].values
                row = [str(x) for x in row]
                row = [x for x in row if x != "nan"]
                
            oitp.write(" \t".join(row))
            oitp.write("\n")            
        oitp.write("\n")    
    oitp.close()
```
# Other functions

#### Calculate angle from 3 cartesian coordinates
```python
import numpy as np
def angle(pointA, pointB, pointC):
    x1x2s = np.power((pointA[0] - pointB[0]),2)
    x1x3s = np.power((pointA[0] - pointC[0]),2)
    x2x3s = np.power((pointB[0] - pointC[0]),2)
    
    y1y2s = np.power((pointA[1] - pointB[1]),2)
    y1y3s = np.power((pointA[1] - pointC[1]),2)
    y2y3s = np.power((pointB[1] - pointC[1]),2)

    cosine_angle = np.arccos((x1x2s + y1y2s + x2x3s + y2y3s - x1x3s - y1y3s)/(2*np.sqrt(x1x2s + y1y2s)*np.sqrt(x2x3s + y2y3s)))

    return np.degrees(cosine_angle)

#check angle function
x = np.array([[  -5.63164 ,      -1.44837   ,     0.00000],
    [  -4.38630    ,    2.49963   ,    -0.00000],
    [     -2.60073  ,     -1.04151   ,     0.00000]])
assert round(angle(*x), 1) == 44.3, "Angle function broken"
```


#### Calculate dihedral from 4 cartesian coordinates
### 
```python
import numpy as np
def new_dihedral(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product
    https://localcoder.org/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
```

#### Centering a PBC system in MDAnalysis by an atom_group
```python
def center(self):
    for i in range(500):
        self.All.wrap()
        #BoxCOG = All.positions.mean(axis=0)
        BoxCOG = self.U.dimensions[:3]/2
        self.proteinCOG = self.HYF_central_layer.center_of_geometry()
        drift = self.proteinCOG - BoxCOG
        self.All.positions = self.All.positions - drift
        if np.abs(drift).sum() < 0.1:
            #print(W, i, "Converged")
            break
        if i >= 499:
            print("Centering did not converge!")
            sys.exit()
    self.All.wrap()
    self.HYF_allatoms.wrap()
```


#### Other

```python
def readin(fname):
    f = open(fname, 'r')
    content = f.read()
    return content
```

```python
import re
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)
```

##### Replace all instances of 'original' with 'new' in file: 'fname'
```python
def replace_in_file2(fname, original, new):
    content = readin(fname)
    content = content.replace(original, new)
    with open(fname, 'w') as outfile:
        outfile.write(content)
```
# TCL
This section to be flushed out.
```tcl
package require pbctools
package require hbonds 

pbc readxst cphmd_out.xst -molid top -all
animate write dcd wrapped.dcd beg 0 end -1 skip 10 waitfor all

mol addfile AAAADEMC_80.pdb type pdb
mol delrep  0 top


set Bilayer [atomselect top "resname POPC POPS"]

hbonds -sel1 $prot -sel2 $wat -writefile yes -plot no -outdir ./ -outfile hbonds_PW.dat -polar yes 
```

