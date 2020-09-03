#!/bin/tcsh
#PBS -N phi4_bc
#PBS -l walltime=72:00:00
#PBS -m abef
#PBS -j oe
#PBS -M snthomas01@email.wm.edu
#PBS -l nodes=1:vortex:ppn=12

cd $PBS_O_WORKDIR


set N=8
mpiexec -n $N python binder_cumulant.py >& log.out
echo "Finished"
