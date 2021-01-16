#!/bin/bash
# declare a name for this job to be sample_job
#PBS -N FEC_ABC_W
# request the parallel queue for this job
#PBS -q gpuq
# request a total of 4 processors for this job (2 nodes and 2 processors per node)
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1
# Request to run for specified hrs:min:secs
##PBS -l walltime=2000:00:00
# combine PBS standard output and error files
#PBS -j oe
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m bea
# specify your email address
#PBS -M rupak@ipr.res.in
#change to the directory where you submitted the job
cd $PBS_O_WORKDIR
#include the full path to the name of your MPI program
gcc -I/soft/fftw-3.3.3/include -L/home/rupak/bin/lib dcst_3d.c -c -std=c99 -lfftw3 -lm
gfortran -I/soft/fftw-3.3.3/include -L/home/rupak/bin/lib fec.f95 dcst_3d.o -lfftw3 -lm
##export OMP_NUM_THREADS=28
ulimit -s unlimited
./a.out




###PBS -l select=1:ncpus=8:mem=100gb:mpiprocs=1:ompthreads=8
###qstat -f 21137.hn01
