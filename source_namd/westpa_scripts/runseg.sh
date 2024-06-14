#!/bin/bash
#runseg.sh
# Make sure environment is set
#this script is run for each trajectory segment. westpa supplies environment variables unique to each segment.
#west_current_seg_ref: path to where the current trajectory segment's data is stored. Become west_parent_date_ref  if any child segment spawns from this segment.
#west_parent_data_ref: path to file or directory containing data for parent segment.

#source env.sh

# Clean up
#rm -f west.log

# Run w_run
#w_run --work-manager processes "$@" &> west.log!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi
####################### Setup for running teh dynamics ####################
#setup directoryt where data for this segment is stored
cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF
echo $WEST_CURRENT_SEG_DATA_REF

#sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/run_md.py > run_md.py

# Make symbolic links to topology and parameter files. These are not unique to each segment
# QUESTION: Is the pdb file for bounded target file or initial non-bounded basis file
ln -sv $WEST_SIM_ROOT/commonfiles/cg_ABCD_unbound.psf structure.psf
ln -sv $WEST_SIM_ROOT/commonfiles/cg_ABCD_unbound.pdb structure.pdb
ln -sv $WEST_SIM_ROOT/commonfiles/cg_AB.par structure1.par
ln -sv $WEST_SIM_ROOT/commonfiles/cg_CD.par structure2.par
#Also copy the parameter files here

## uncomment the following code when HDF5 framework is off ##
# don't use symbolic links, as the folders are turned into .h5 files
#ln -sv $WEST_SIM_ROOT/common_files/bstate.psf .
#ln -sv $WEST_SIM_ROOT/common_files/bstate.pdb .
#cp $WEST_PARENT_DATA_REF/seg.xml ./parent.xml 
#cp $WEST_PARENT_DATA_REF/seg.dcd ./parent.dcd

# Run the dynamics with OpenMM
# python run_md.py SEED GPU_DEVICE
#echo 'Will run dynamics'
#echo python $WEST_SIM_ROOT/common_files/run_md.py $WEST_RAND16 $WM_PROCESS_INDEX
#python $WEST_SIM_ROOT/common_files/run_md.py $WEST_RAND16 $WM_PROCESS_INDEX


#Run the dynamics with NAMD
#WEighted ensemble  requires dynamics are stochatstic. 
sed "s/RAND/$WEST_RAND16/g" \
  $WEST_SIM_ROOT/commonfiles/md.conf > md.conf

# Trajectory segment will start off where it's parent segment left off. 
ln -sv $WEST_PARENT_DATA_REF/seg.coor ./parent.coor
ln -sv $WEST_PARENT_DATA_REF/seg.dcd  ./parent.dcd
ln -sv $WEST_PARENT_DATA_REF/seg.vel  ./parent.vel
ln -sv $WEST_PARENT_DATA_REF/seg.xsc  ./parent.xsc

# Run the dynamics: Propagate the segment using namd2
$NAMD md.conf > seg.log
echo $PWD
#Calculate and return data
python $WEST_SIM_ROOT/commonfiles/dist.py > $WEST_PCOORD_RETURN





#Calculate pcoord with MDAnalysis
#python $WEST_SIM_ROOT/commonfiles/dist.py > $WEST_PCOORD_RETURN


#cp bstate.pdb $WEST_TRAJECTORY_RETURN
#cp seg.dcd $WEST_TRAJECTORY_RETURN

#cp bstate.pdb $WEST_RESTART_RETURN
#cp seg.xml $WEST_RESTART_RETURN/parent.xml

#cp seg.log $WEST_LOG_RETURN

# Clean up
#rm -f dist.dat run_md.py
rm -f md.conf seg.pdb \
  seg.restart.coord seg.restart.coor.old seg.restart.vel seg.restart.vel.old\
  seg.restart.xsc seg.restart.xsc.old structure.psf
