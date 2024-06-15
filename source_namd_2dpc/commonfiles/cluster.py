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
from numpy.linalg import norm
#import os, sys

#this=os.getcwd()
#print(this)
#
# Load the trajectory.

u = mda.Universe('../../../commonfiles/ABCD_unbound_neutral.psf', 'seg.dcd')
u_parent = mda.Universe('../../../commonfiles/ABCD_unbound_neutral.psf', 'parent.dcd')
meanarr=[10.801121726003823, 8.207154846176936, 23.319908330398302, 18.38726461897426]
#print(meanarr[0])
#print(meanarr[3])
variancearr=[0.3794918332730437, 0.5325791974516593, 1.3844753811965562, 1.006144187649773]

def gaussfn(x,meang,varianceg):
    return np.exp(-(x-meang)**2/(2*varianceg))/(np.sqrt(2*np.pi*varianceg))

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
    #print(distcl)
    #print(meanarr[0])
    #print(linearfn(distcl,meanarr[0],variancearr[0]))
    return linearfn(distcl,meanarr[0],variancearr[0])

def cluster_gaussian2(u):
    AB_cluster=u.select_atoms("bynum 51 136 163")
    CD_cluster=u.select_atoms("bynum 327")
    comAB=AB_cluster.center_of_mass()
    comCD=CD_cluster.center_of_mass()
    distcl=norm(comCD-comAB)
    #print(distcl)
    #print(meanarr[1])
    #print(linearfn(distcl,meanarr[1],variancearr[1]))
    return linearfn(distcl,meanarr[1],variancearr[1])


def cluster_gaussian3(u):
    AB_cluster=u.select_atoms("bynum 122 97 55")
    CD_cluster=u.select_atoms("bynum 327")
    comAB=AB_cluster.center_of_mass()
    comCD=CD_cluster.center_of_mass()
    distcl=norm(comCD-comAB)
   # print(gaussfn(distcl,meanarr[2],variancearr[2]))
    return linearfn(distcl,meanarr[2],variancearr[2])

def cluster_gaussian4(u):
    AB_cluster=u.select_atoms("bynum 124 103")
    CD_cluster=u.select_atoms("bynum 327")
    comAB=AB_cluster.center_of_mass()
    comCD=CD_cluster.center_of_mass()
    distcl=norm(comCD-comAB)
    #print(gaussfn(distcl,meanarr[3],variancearr[3]))
    return linearfn(distcl,meanarr[3],variancearr[3])

def calcDistUniv(u, trajBeg=0, trajEnd=None):
	for ts in u.trajectory[trajBeg:trajEnd]:
             print(cluster_gaussian1(u)+cluster_gaussian2(u)+cluster_gaussian3(u)+cluster_gaussian4(u)) 
# first do parent universe.  
# for parent, only do the very last trajectory
calcDistUniv(u_parent, trajBeg=-1)
calcDistUniv(u)
