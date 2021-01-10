#!/bin/bash
#PBS -N 512_GPU
#PBS -q serialq 
##PBS -l select=1:ncpus=2:ngpus=1
#PBS -l select=1:ncpus=2:mpiprocs=2:ngpus=2
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR

module load pgi/18.10 
module load cuda91/toolkit/9.1.85
#nvcc  -gencode arch=compute_30,code=sm_30 -c cufft3D.cu
nvcc  -gencode arch=compute_60,code=sm_60 -c cufft2D.cu

#/soft/pgicdk-1510/linux86-64/15.10/bin/pgf95 -acc -c -Minfo=accel -ta=tesla:cc30 -I/soft/pgicdk-1510/linux86-64/15.10/include/ PSMHD3-cufftacc.f95 -lm
pgf95 -fast -c -Minfo=accel -ta=tesla:cc60,managed PSHD2-cufftacc.f95 -lm

#/soft/pgicdk-1510/linux86-64/15.10/bin/pgf95 -acc -Minfo=accel -ta=tesla:cc30 -I/soft/pgicdk-1510/linux86-64/15.10/include/  -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm
pgf95 -fast -Minfo=accel -ta=tesla:cc60,managed -Mcudalib=cufft PSHD2-cufftacc.o cufft2D.o -lm

./a.out

