#!/usr/bin/env/python

import MDAnalysis as mda
#from MDAnalysis.tests.datafiles import PSF, DCD, CRD
from MDAnalysis.analysis import rms,contacts
import MDAnalysis.lib.util
import numpy as np
import pickle
import os, sys

this=os.getcwd()

u = mda.Universe('../../../commonfiles/ABCD_unbound_neutral.pdb', 'seg.dcd')
u_parent = mda.Universe('../../../commonfiles/ABCD_unbound_neutral.pdb', 'parent.dcd')
ds = pickle.load(open( "../../../commonfiles/contactIndices.npy", "rb" ) )
contref = ds['contref']
def number_of_native_contacts(univ,cont, trajBeg=0, trajEnd=None):
    for ts in univ.trajectory[trajBeg:trajEnd]:
        group1=univ.select_atoms('segid P1')
        group2=univ.select_atoms('segid P2')
        dist= contacts.distance_array(group1.positions,group2.positions)
    #totalcontacts=len(AB_array)
        distances_in_contacts=dist[cont[0],cont[1]]
        total_number_of_contacts=len(cont[0])
        threshold=10
        count = np.count_nonzero(distances_in_contacts < threshold)
        print(count*1.0/total_number_of_contacts)

number_of_native_contacts(u_parent,contref,trajBeg=-1)
number_of_native_contacts(u,contref)


