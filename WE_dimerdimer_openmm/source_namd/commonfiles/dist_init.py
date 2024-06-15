#!/usr/bin/env/python
#dist_list.py
#Get a list of distances between two dimers that is less than 20 A
import MDAnalysis as mda
from MDAnalysis.analysis import rms, contacts
import MDAnalysis.lib.util
import numpy as np
import pickle
#Okay we nee dto have a reference structure for which we are calculating the distance (can be the basis state or whatever) and an original crystal structure of the two dimers.
# First step is to calculate the distance between two dimers and note down the dimers that have smaller distance between them then a cutoff in the crystal structure.
#Then we store these contacts. Calculate the distance matrix for the given structure and then calculate the mean distance for the same contacts.

f='./'
ref=mda.Universe(f+'cg_ABCD_bound.pdb')
u=mda.Universe(f+'cg_ABCD_unbound.pdb')
beads_AB=u.select_atoms('segid P1')
beads_CD=u.select_atoms('segid P2')
beads_AB_ref=ref.select_atoms('segid P1')
beads_CD_ref=ref.select_atoms('segid P2')
#print(beads_AB_ref)
#print(beads_CD_ref)

def contacts_within_cutoff(univ, group1, group2, cutoff_dist):
	dist= contacts.distance_array(group1.positions,group2.positions)
	cont= np.where(dist<cutoff_dist)
	dist_print=dist[cont[0],cont[1]]
#TODO: output the inetrfacial beads in a new pdb
#TODO: add mean of distances here
	meanref = dist[cont[0],cont[1]].mean()
	#print ((len(cont[0]))
	#print (len(cont[1]))
	#print (len(dist_print))
	return cont, dist_print, meanref

#beads, distance = contacts_within_cutoff(ref, beads_AB, beads_CD, 20)
#length_dis=len(distance)
#print(length_dis)
#print ("cg_AB,cg_CD,dist \n")
#for i in range(length_dis):
#	print (beads[0][i],beads[1][i],distance[i]) 
#	print ("\n")
#rad = 9

rad=12.5
contref, distref, meanref = contacts_within_cutoff(ref, beads_AB_ref, beads_CD_ref, rad)
#print('meanBD', meanBD) # 7.62
#print('contBD', contBD)
#print('distBD', distBD.shape, distBD)

# contInd - (indices_group1, indices_group2) tuple of indices of c-alpha atoms from which to compute RMSD
# dist - distance matrix
def distBetweenContacts(contInd, dist):
  return dist[contInd[0], contInd[1]].mean()

# distance in the reference protein, i.e. cristallographic structure with formed complex
beads_AB_ref.atoms[contref[0]].write('contactsB.pdb')
beads_CD_ref.atoms[contref[1]].write('contactsD.pdb')
with open("./commonfiles/contactIndices.npy", "wb" ) as file:
  pickle.dump(dict(contref=contref), file)
print(f)
with open("./commonfiles/contactIndices.npy", "rb" ) as file:
  ds = pickle.load(file)
contref = ds['contref']

# distances in the current frame
distBD = contacts.distance_array(beads_AB.positions, beads_CD.positions)
#print('distBD', distBD.shape, distBD)

# distances at pre-selected contacts
meanBD = distBetweenContacts(contref, distBD)

#print(meanAB, meanBC,  meanCD)
#print(meanAB + meanBC + meanCD)
print(meanBD)
#print(meanref)
