#!/usr/bin/env/python
#com_init.py
#Calculate center of mass distance between two dimers
import MDAnalysis as mda
from MDAnalysis.analysis import rms, contacts
import MDAnalysis.lib.util
import numpy as np
import pickle
from numpy.linalg import norm


def com_dist_dimers(u):
    dimer1=u.select_atoms('segid P1')
    dimer2=u.select_atoms('segid P2')
    return norm(dimer1.center_of_mass()-dimer2.center_of_mass())
ref=mda.Universe('../cg_ABCD-neutral.psf','../cg_ABCD_neutral.pdb')
initial=mda.Universe('ABCD_unbound_neutral.psf','ABCD_unbound_neutral.pdb')
com_ref=com_dist_dimers(ref)
com_intial=com_dist_dimers(initial)
print(com_ref)
print(com_intial)
