#!/bin/bash
#Give your job name
#BSUB -J PSMHD3D

#Give required number of processor
#BSUB -n 32

#Give your email id
##BSUB -u <user-mail-id@ipr.res.in>

#Give your current working directory
##BSUB -cwd ~/<your home area>

#Redirect error to some file and %J is Job ID
#BSUB -e err.%J

#Redirect output to some file and %J is Job ID
#BSUB -o out.%J

#Define number of processor per host user wants to execute a job
#BSUB -R "span[hosts=1]"

#cd /home/rupak/udbhav/My_PhD/FFTW/3D_MHD/Subroutine
make -f MHD3D.mak
export OMP_NUM_THREADS=32
./PSMHD3
