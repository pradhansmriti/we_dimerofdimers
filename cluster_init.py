#!/usr/bin/env/python
#dist_list.py
#Calculate distance for each individual cluster 
#Then calculate sum of Gaussian distances
import MDAnalysis as mda
from MDAnalysis.analysis import rms, contacts
import MDAnalysis.lib.util
import numpy as np
import pickle
from numpy.linalg import norm
#Okay we nee dto have a reference structure for which we are calculating the distance (can be the basis state or whatever) and an original crystal structure of the two dimers.
# First step is to calculate the distance between two dimers and note down the dimers that have smaller distance between them then a cutoff in the crystal structure.
#Then we store these contacts. Calculate the distance matrix for the given structure and then calculate the mean distance for the same contacts.

f='./'
#ref=mda.Universe(f+'cg_ABCD_bound_final.pdb')
#u=mda.Universe(f+'cg_ABCD_unbound.pdb')
#beads_AB=u.select_atoms('segid P1')
#beads_CD=u.select_atoms('segid P2')
#beads_AB_ref=ref.select_atoms('segid P1')
#beads_CD_ref=ref.select_atoms('segid P2')
#print(beads_AB_ref)
#print(beads_CD_ref)
meanarr=[10.801121726003823, 8.207154846176936, 23.319908330398302, 18.38726461897426]
#print(meanarr[0])
#print(meanarr[3])
variancearr=[0.3794918332730437, 0.5325791974516593, 1.3844753811965562, 1.006144187649773]

def gaussfn(x,meang,varianceg):
    return np.exp(-(x-meang)**2/(2*varianceg))

def quadfn(x,meang,varianceg):
    return (x-meang)**2/(2*varianceg)

def linearfn(x,meang,varianceg):
    return abs(x-meang)/(np.sqrt(varianceg))
def cluster_gaussian1(u):
    AB_cluster=u.select_atoms("bynum 51 136 163")
    CD_cluster=u.select_atoms("bynum 356 470")
    comAB=AB_cluster.center_of_mass()
    comCD=CD_cluster.center_of_mass()
    distcl=norm(comCD-comAB)
    print(distcl)
    print(meanarr[0])
    #print(linearfn(distcl,meanarr[0],variancearr[0]))
    return gaussfn(distcl,meanarr[0],variancearr[0])

def cluster_gaussian2(u):
    AB_cluster=u.select_atoms("bynum 51 136 163")
    CD_cluster=u.select_atoms("bynum 327")
    comAB=AB_cluster.center_of_mass()
    comCD=CD_cluster.center_of_mass()
    distcl=norm(comCD-comAB)
    print(distcl)
    print(meanarr[1])
    #print(linearfn(distcl,meanarr[1],variancearr[1]))
    return gaussfn(distcl,meanarr[1],variancearr[1])


def cluster_gaussian3(u):
    AB_cluster=u.select_atoms("bynum 122 97 55")
    CD_cluster=u.select_atoms("bynum 327")
    comAB=AB_cluster.center_of_mass()
    comCD=CD_cluster.center_of_mass()
    distcl=norm(comCD-comAB)
   # print(gaussfn(distcl,meanarr[2],variancearr[2]))
    return gaussfn(distcl,meanarr[2],variancearr[2])

def cluster_gaussian4(u):
    AB_cluster=u.select_atoms("bynum 124 103")
    CD_cluster=u.select_atoms("bynum 327")
    comAB=AB_cluster.center_of_mass()
    comCD=CD_cluster.center_of_mass()
    distcl=norm(comCD-comAB)
    #print(gaussfn(distcl,meanarr[3],variancearr[3]))
    return gaussfn(distcl,meanarr[3],variancearr[3])

ref=mda.Universe('cg_ABCD-neutral.psf','cg_ABCD_neutral.pdb')
u=mda.Universe('ABCD_unbound_neutral.psf','ABCD_unbound_neutral.pdb')
#dist_U=cluster_gaussian1(u)+cluster_gaussian2(u)+cluster_gaussian3(u)+cluster_gaussian4(u)
dist_ref=cluster_gaussian1(ref)+cluster_gaussian2(ref)+cluster_gaussian3(ref)+cluster_gaussian4(ref)
print(dist_ref)
#print(dist_U)
# distances in the current frame
#distBD = contacts.distance_array(beads_AB.positions, beads_CD.positions)
#print('distBD', distBD.shape, distBD)

# distances at pre-selected contacts
#meanBD = distBetweenContacts(contref, distBD)

#print(meanAB, meanBC,  meanCD)
#print(meanAB + meanBC + meanCD)
#print(meanBD)
#print(meanref)
