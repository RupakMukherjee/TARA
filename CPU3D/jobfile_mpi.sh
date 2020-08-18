#!/bin/bash
# declare a name for this job to be sample_job
#PBS -N Test
# request the parallel queue for this job
#PBS -q longq
# request a total of 4 processors for this job (2 nodes and 2 processors per node)
##PBS -l nodes=2:ppn=1 
#PBS -l select=110:ncpus=1:mpiprocs=1
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
mpif90 -I/home/rupak/bin/include mhd_3d_mpi.f95 -L/home/rupak/bin/lib -lfftw3_mpi -lfftw3 -lm 
#mpif90 -I/home/rupak/bin/include mytest_MPI_transposed.f90 -L/home/rupak/bin/lib -lfftw3_mpi -lfftw3 -lm 
# gfortran -fopenmp -I/soft/fftw-3.3.3/include -L/home/rupak/bin/lib mhd_3d.f95 -lfftw3_omp -lfftw3 -lm
#mpif90 -fopenmp -I/home/rupak/bin/include -I/soft/fftw-3.3.3/include mhd_3d_mpi.f95 -L/home/rupak/bin/lib -lfftw3_omp -lfftw3_mpi -lfftw3 -lm 
# export OMP_NUM_THREADS=26
ulimit -s unlimited
mpirun -np 110 ./a.out
# ./a.out




###PBS -l select=1:ncpus=8:mem=100gb:mpiprocs=1:ompthreads=8
###qstat -f 21137.hn01
