#!/usr/bin/env/python
#nativefraction_init.py
#Calculate fraction of native contacts between two dimers
import MDAnalysis as mda
from MDAnalysis.analysis import rms, contacts
import MDAnalysis.lib.util
import numpy as np
import pickle
from numpy.linalg import norm


f='./'
ref=mda.Universe(f+'cg_ABCD_neutral.pdb')
u=mda.Universe(f+'ABCD_unbound_neutral.pdb')
beads_AB=u.select_atoms('segid P1')
beads_CD=u.select_atoms('segid P2')
beads_AB_ref=ref.select_atoms('segid P1')
beads_CD_ref=ref.select_atoms('segid P2')

def contacts_within_cutoff(univ, cutoff_dist):
    group1=univ.select_atoms('segid P1')
    group2=univ.select_atoms('segid P2')
    dist= contacts.distance_array(group1.positions,group2.positions)
    print(dist[6,86])
    cont= np.where(dist<cutoff_dist)
    dist_print=dist[cont[0],cont[1]]
    print(cont[0])
    print(cont[1])
    print(dist_print)
    meanref = dist[cont[0],cont[1]].mean()
    return cont
cutoffdist=10
contref=contacts_within_cutoff(ref,cutoffdist)
beads_AB_ref.atoms[contref[0]].write('contactsB.pdb')
beads_CD_ref.atoms[contref[1]].write('contactsD.pdb')
with open("./commonfiles/contactIndices.npy", "wb" ) as file:
  pickle.dump(dict(contref=contref), file)
print(f)
with open("./commonfiles/contactIndices.npy", "rb" ) as file:
  ds = pickle.load(file)
contref = ds['contref']
def fraction_of_native_contacts(univ,cont):
    group1=univ.select_atoms('segid P1')
    group2=univ.select_atoms('segid P2')
    dist= contacts.distance_array(group1.positions,group2.positions)
    #totalcontacts=len(AB_array)
    distances_in_contacts=dist[cont[0],cont[1]]
    total_number_of_contacts=len(cont[0])
    threshold=10
    count = np.count_nonzero(distances_in_contacts < threshold)
    return count*1.0/total_number_of_contacts
print(fraction_of_native_contacts(ref,contref))
print(fraction_of_native_contacts(u,contref))
