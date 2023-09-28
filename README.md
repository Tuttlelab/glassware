# Tuttlelab/glassware
## A repo for frequently used code snippets

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

# peptideutils

Some tools for working with peptides in python

```python
import peptideutils

peptideutils.GenerateDatasetIndex(3)
> ["AAA", "AAC", ... "YYY"]

peptideutils.charge("AKRDF")
> 1
```

# to add
parsepsf
MDAnalysis centering


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
if round(angle(*x), 1) != 44.3:
    print("Angle function broken")
    sys.exit()
    
```


#### Calculate dihedral from 3 cartesian coordinates
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

