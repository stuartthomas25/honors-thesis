#!/bin/tcsh
#PBS -N gradientflow
#PBS -l walltime=48:00:00
#PBS -m abef
#PBS -j oe
#PBS -M snthomas01@email.wm.edu
#PBS -l nodes=1:vortex:ppn=12



set N=1
setenv OPENMPI_HOME /usr/local/intel64/broadwell/intel-2017/openmpi-3.1.4-ib
setenv LD_LIBRARY_PATH /usr/local/intel64/broadwell/intel-2017/openmpi-3.1.4-ib/lib:/usr/local/intel-2017/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin:/usr/local/intel-2017/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64:/usr/local/intel-2017/compilers_and_libraries_2017.0.098/linux/mkl/lib/intel64:/usr/local/intel-2017/compilers_and_libraries_2017.0.098/linux/tbb/lib/intel64/gcc4.4:/usr/local/intel-2017/debugger_2017/iga/lib:/usr/local/intel-2017/debugger_2017/libipt/intel64/lib:/usr/local/intel-2017/compilers_and_libraries_2017.0.098/linux/daal/lib/intel64_lin:/usr/local/intel-2017/compilers_and_libraries_2017.0.098/linux/daal/../tbb/lib/intel64_lin/gcc4.4:/usr/local/gcc-5.2.0/lib:/usr/local/gcc-5.2.0/lib64:/lib64:/usr/lib64:/usr/local/torque-6.1.1.1/lib:/usr/local/anaconda3-2020.02/lib

if ($?PBS_O_WORKDIR) then
    cd $PBS_O_WORKDIR
else
    echo 'Running on head node'
endif

foreach input (inputs/*)
    set output="outputs/`basename $input yml`csv"
    #nohup mpiexec -n $N ./cpp_code/bin/sweep -i $input -o $output -q >! log.out &
    ./cpp_code/bin/sweep -i $input -o $output -q >& log.out &
end
#mpiexec -n $N python gradient_flow2.py 48 >& log.out

wait
echo "Finished"
