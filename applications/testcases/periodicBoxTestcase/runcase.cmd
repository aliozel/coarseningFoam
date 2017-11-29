#!/bin/bash

# parallel job using 32 processors. and runs for 1 hours (max) 
#SBATCH -N 8 
#SBATCH --ntasks-per-node=16 
#SBATCH -t 144:00:00
#SBATCH -J alpp0.10LongRun

# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
# #SBATCH --mail-type=begin
# #SBATCH --mail-type=end
# #SBATCH --mail-user=aozel@princeton.edu

# Source Openfoam
source ~/.bashrc

# Source modules
module load openmpi
module load intel/14.0

# Case folder
caseFolder=$PWD

# Casename
caseName=$(basename $PWD)

# Create folder at node 
#export ScratchDir="/scratch/aozel"
#rm -rf $ScratchDir 
#mkdir -p $ScratchDir

# Copy files to node
#rsync -auvz $caseFolder $ScratchDir
#cp -r $caseFolder $ScratchDir

# Run
#cd $ScratchDir/$caseName	
#echo "Running job folder:" $PWD
source $HOME/OpenFOAM/OpenFOAM-2.2.2/etc/bashrc
decomposePar 
./parCFDDEMrun.sh

# Copy results files to  
#rsync -auvz $ScratchDir/$caseName $caseFolder/../
#cd $ScratchDir
#rm -rf $caseName
