#!/usr/bin/env/python
import MDAnalysis as mda
#from MDAnalysis.tests.datafiles import PSF, DCD, CRD
from MDAnalysis.analysis import rms,contacts
import MDAnalysis.lib.util
import numpy as np
import pickle
from numpy.linalg import norm
u = mda.Universe('../../../commonfiles/ABCD_unbound_neutral.psf', 'seg.dcd')
u_parent = mda.Universe('../../../commonfiles/ABCD_unbound_neutral.psf', 'parent.dcd')


def com_dist_dimers(u):
    dimer1=u.select_atoms('segid P1')
    dimer2=u.select_atoms('segid P2')
    return norm(dimer1.center_of_mass()-dimer2.center_of_mass())
def calcComDistUniv(u, trajBeg=0, trajEnd=None):
    for ts in u.trajectory[trajBeg:trajEnd]:
             print(com_dist_dimers(u))

calcComDistUniv(u_parent, trajBeg=-1)
calcComDistUniv(u)
