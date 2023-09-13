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
