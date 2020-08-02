Multi GPU Card Script
===============================================

module load PrgEnv/PGI+OpenMPI/2018-03-29
module load pgi
module load cuda
module load openmpi


export ACCFFT_DIR=/home/nvydyanathan/Work/accfft/accfft-aug20-2018/accfft-master/build/
export FFTW_DIR=/home/nvydyanathan/Work/fftw3.3.8
mpic++ -c accfft-module.C -I $ACCFFT_DIR/include/ -I $FFTW_DIR/include/
mpif90 -acc -Minfo=accel -ta=tesla:cc70,managed -c PSMHD3-accfft_512.f95 
mpif90 -acc -Minfo=accel -ta=tesla:cc70,managed PSMHD3-accfft_512.o accfft-module.o -o PSMHD3mGPU -laccfft -laccfft_utils -laccfft_gpu -laccfft_utils_gpu -L $ACCFFT_DIR/lib -lfftw3_threads -lfftw3_mpi -lfftw3 -L $FFTW_DIR/lib -lgomp -lmpi_cxx -lm -lstdc++ -Mcudalib=cufft


NOTE: For 512^3 grid size, you still run out of GPU memory even on 4 GPUs to allocate all of the arrays. So, you need to use managed memory. But the code will run faster than on a single GPU. On a server with 4 Tesla V100, connected via PCI-e, the code gives a 3.2x speedup over a single GPU run. If you change the code for in-place ffts, time can further be reduced.  

srun --exclusive -N 1 -n 4 -t 08:00:00 -p dgx-1v_32g --pty "bash" or replace dgx-1v_32g with dgx-2v_32g

srun --exclusive -N 1 -n 4 -t 08:00:00 -p hsw_v100_32g --pty "bash" – this gets you node on the V100 cluster – each node has 4 V100 GPUs via PCI-e

mpirun -np <number of mpi processes> ./<exename>
