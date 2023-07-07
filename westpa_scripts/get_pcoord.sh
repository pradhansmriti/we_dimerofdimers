#!/bin/bash
#
# get_pcoord.sh
#
# This script is run when calculating initial progress coordinates for new 
# initial states (istates).  This script is NOT run for calculating the progress
# coordinates of most trajectory segments; that is instead the job of runseg.sh.

# If we are debugging, output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

# Make sure we are in the correct directory
cd $WEST_SIM_ROOT
source env.sh
echo $PWD
# Make a temporary file in which to store output from the python script
DIST=$(mktemp)
echo $DIST
# Use a custom python script to calculate the distance between the Na+ and Cl-
# ions. This script looks for files named 'nacl.psf' and 'seg.dcd'.
python $WEST_SIM_ROOT/commonfiles/cluster_init.py > $DIST
echo $DIST

# Pipe the relevant part of the output file (the distance) to $WEST_PCOORD_RETURN
# WEST_PCOORD_RETURN=$(tail -n 1 $DIST)
cat $DIST | tail -n 1 > $WEST_PCOORD_RETURN
echo $WEST_PCOORD_RETURN
# Remove the temporary file to clean up
#rm $DIST

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
