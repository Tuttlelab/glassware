# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 00:30:25 2023

@author: rkb19187
"""
#import cupy as cp
import numpy as np

class rmsd_tools:
    def rotation_matrix_from_points(self, m0, m1):
        """Returns a rigid transformation/rotation matrix that minimizes the
        RMSD between two set of points.
    
        m0 and m1 should be (3, npoints) numpy arrays with
        coordinates as columns::
    
            (x1  x2   x3   ... xN
             y1  y2   y3   ... yN
             z1  z2   z3   ... zN)
    
        The centeroids should be set to origin prior to
        computing the rotation matrix.
    
        The rotation matrix is computed using quaternion
        algebra as detailed in::
    
            Melander et al. J. Chem. Theory Comput., 2015, 11,1055
        """
    
        v0 = np.copy(m0)
        v1 = np.copy(m1)
    
        # compute the rotation quaternion
    
        R11, R22, R33 = np.sum(v0 * v1, axis=1)
        R12, R23, R31 = np.sum(v0 * np.roll(v1, -1, axis=0), axis=1)
        R13, R21, R32 = np.sum(v0 * np.roll(v1, -2, axis=0), axis=1)
    
        f = [[R11 + R22 + R33, R23 - R32, R31 - R13, R12 - R21],
             [R23 - R32, R11 - R22 - R33, R12 + R21, R13 + R31],
             [R31 - R13, R12 + R21, -R11 + R22 - R33, R23 + R32],
             [R12 - R21, R13 + R31, R23 + R32, -R11 - R22 + R33]]
    
        F = np.array(f)
    
        w, V = np.linalg.eigh(F)
        # eigenvector corresponding to the most
        # positive eigenvalue
        q = V[:, np.argmax(w)]
    
        # Rotation matrix from the quaternion q
    
        R = self.quaternion_to_matrix(q)
    
        return R
    
    
    def quaternion_to_matrix(self, q):
        """Returns a rotation matrix.
    
        Computed from a unit quaternion Input as (4,) numpy array.
        """
    
        q0, q1, q2, q3 = q
        R_q = [[q0**2 + q1**2 - q2**2 - q3**2,
                2 * (q1 * q2 - q0 * q3),
                2 * (q1 * q3 + q0 * q2)],
               [2 * (q1 * q2 + q0 * q3),
                q0**2 - q1**2 + q2**2 - q3**2,
                2 * (q2 * q3 - q0 * q1)],
               [2 * (q1 * q3 - q0 * q2),
                2 * (q2 * q3 + q0 * q1),
                q0**2 - q1**2 - q2**2 + q3**2]]
        return np.array(R_q)
    
    def minimize_rotation_and_translation(self, p0, p):
        """Minimize RMSD between atoms and target.
    
        Rotate and translate atoms to best match target.  For more details, see::
    
            Melander et al. J. Chem. Theory Comput., 2015, 11,1055
        """
    
        # centeroids to origin
        c = np.mean(p, axis=0)
        p -= c
        c0 = np.mean(p0, axis=0)
        p0 -= c0
    
        # Compute rotation matrix
        R = self.rotation_matrix_from_points(p.T, p0.T)
    
        return np.dot(p, R.T) + c0
    
    def RMSD(self, p0, p):
        return np.sqrt(np.mean(np.linalg.norm(p0 - p, axis=1)**2))
           
    def __init__(self, gpu=False):
        self.gpu = gpu
    
if __name__ == "__main__":
    import time
    RMSD = rmsd_tools()
    
    st = time.time()
    Test = np.random.random((10000, 18, 3))
    
    for i in range(1, Test.shape[0]):
        RMSD.minimize_rotation_and_translation(Test[0], Test[i])
    
    rmsds = np.ones((Test.shape[0], Test.shape[0]))
    
    for i in range(0, Test.shape[0]):
        rmsds[i] = np.sqrt(np.mean(np.linalg.norm(Test[i] - Test, axis=2)**2, axis=1))
    print(rmsds.round(3))
    
    print("Took:", round(time.time()-st, 2), "s")