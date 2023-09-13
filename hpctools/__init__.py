#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 22:37:42 2020

@author: rkb19187
"""

import subprocess
import numpy as np

class hpctools:
    def readin(self, fpath):
        f = open(fpath)
        content = f.read()
        f.close()
        return content
    
    def read_sqme(self):
         self.sqme = str(subprocess.check_output(['squeue', '-u', 'rkb19187', '-o', "%j"]))
    def nGPU(self):
         total = str(subprocess.check_output(['squeue', '-o', "%P"]))
         mine  = str(subprocess.check_output(['squeue', '-u', 'rkb19187', '-o', "%P"]))
         total = total.split("\\n")
         mine = mine.split("\\n")
         return {"Total": total.count("gpu"), "Mine": mine.count("gpu")}
    
    @property
    def ID(self,):
        self.jobs = subprocess.check_output(['squeue', '-u', 'rkb19187', '-o', '"%j -_- %i"'])
        self.jobs = self.jobs.decode()
        IDs = {}
        for line in self.jobs.split("\n")[1:]:
            if len(line.split(" -_- ")) != 2:
                continue
            key = line.split(" -_- ")[0]
            key = key.replace('"', '')
            val = line.split(" -_- ")[1]
            val = int(val.replace('"', ''))
            IDs[key] = val
        return IDs
     
    def IsRunning(self, JobID, exact=False):
        self.read_sqme()
        if exact:
            return JobID in self.sqme.split("\\n")
        else:
            for line in self.sqme.split("\\n"):
                if JobID in line:
                    return True
            return False
    def count_jobs(self):
        self.cores = subprocess.check_output(['squeue', '-u', os.environ["USERNAME"], '-o', "%C"])
        self.cores = self.cores.decode("utf-8")
        self.jobs = len([int(x) for x in self.cores.split("\n")[1:-1]])
        self.cores = sum([int(x) for x in self.cores.split("\n")[1:-1]]) # first line is 'NODES', last is empty
        return {"jobs":self.jobs, "cores":self.cores}
        
# =============================================================================
#         self.read_sqme()
#         sqme = self.sqme.split()
#         sqme =  np.array(sqme)
#         self.nJobs = np.where(sqme == "standard")[0].shape[0]
#         self.waiting = np.where(sqme == "0:00")[0].shape[0]
#         self.running = self.nJobs-self.waiting
# 
#         return (self.running, self.waiting)
# =============================================================================
        
    def __init__(self):
        #self.tmpfname = "/users/rkb19187/hpctools.temp.out"
        pass

if __name__ == "__main__":
    h = hpctools()
    #print(h.IsRunning("CC_C_CC_C__C_N"))
    print(h.count_jobs())
    print(h.ID["Prot_MD_001862_3_s"])
    #print(h.ID)

