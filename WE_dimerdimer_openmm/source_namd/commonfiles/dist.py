#!/usr/bin/env/python
#
# dist.py
#
# Calculate distance between key selected atoms at the interface of the protein. 
# We need all timepoints from the current segment, and the final timepoint
# from the parent segment (which is where the current segment starts)
# REQUIREMENT: this script should print, for each frame, either one number (e.g. 11.801), or 3 numbers (for 3D progress coordinate separated by space) 
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, CRD
from MDAnalysis.analysis import rms,contacts
import MDAnalysis.lib.util
import numpy as np
import pickle
import os, sys

this=os.getcwd()
#print(this)
#
# Load the trajectory.

u = mda.Universe('../../../commonfiles/cg_ABCD_unbound.pdb', 'seg.dcd')
u_parent = mda.Universe('../../../commonfiles/cg_ABCD_unbound.pdb', 'parent.dcd')
#print(u.atoms)
#print(u.atoms.segments)

# contInd - (indices_group1, indices_group2) tuple of indices of c-alpha atoms from which to compute RMSD
# dist - distance matrix
def distBetweenContacts(contInd, dist):
  #print(dist[contInd[0], contInd[1]])
  return dist[contInd[0], contInd[1]].mean()
  
# load a-priori chosen atoms at protein interfaces, for which we'll compute the progress coordinates. 
# contAB = (indices_A, indices_B) are the indices of ~10-20 atoms at the interface of A and B
ds = pickle.load(open( "../../../commonfiles/contactIndices.npy", "rb" ) )
#contAB = ds['contAB']
#contBC = ds['contBC']
#contCD = ds['contCD']
contref = ds['contref']

def calcDistUniv(u, trajBeg=0, trajEnd=None):
	for ts in u.trajectory[trajBeg:trajEnd]:
		# MDAnalysis automatically sets the universe to the current traj. 
		# For memory efficiency, MDA loads each frame as requested ... not everything at once

		#print(u)
		#print(u.atoms.segments)
		#ca_D = u.select_atoms('name CA and segid D') # to check the mapping, see 2g33_init.pdb
		#ca_C = u.select_atoms('name CA and segid C')
		#ca_B = u.select_atoms('name CA and segid B')
		#ca_A = u.select_atoms('name CA and segid A')
		beads_AB=u.select_atoms('segid P1')
		beads_CD=u.select_atoms('segid P2')
		# distances in the current frame
		#distAB = contacts.distance_array(ca_A.positions, ca_B.positions)
		#distBC = contacts.distance_array(ca_B.positions, ca_C.positions)
		#distCD = contacts.distance_array(ca_C.positions, ca_D.positions)
		distBD = contacts.distance_array(beads_AB.positions, beads_CD.positions)

		#print('contAB', contAB)
		#print('distAB', distAB)
		# distances at pre-selected contacts
		#meanAB = distBetweenContacts(contAB, distAB) # 6A avg dist
		#meanBC = distBetweenContacts(contBC, distBC) # 8A
		#meanCD = distBetweenContacts(contCD, distCD) # 
		meanBD = distBetweenContacts(contref, distBD) # 8A

		#print(meanAB + meanBC + meanCD)
		#print(meanAB)
		#print(meanBC)
		#print(meanCD)
		print(meanBD)
# first do parent universe.  
# for parent, only do the very last trajectory
calcDistUniv(u_parent, trajBeg=-1)
calcDistUniv(u)
