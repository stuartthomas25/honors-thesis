#!/bin/tcsh
#PBS -N phi4_bc
#PBS -l walltime=72:00:00
#PBS -m abef
#PBS -j oe
#PBS -M snthomas01@email.wm.edu
#PBS -l nodes=1:vortex:ppn=12

#if ( ! (-e $PBS_O_WORKDIR) ) then
    #cd $PBS_O_WORKDIR


set N=12
mpiexec -n $N python gradient_flow2.py 48 >& log.out
echo "Finished"
