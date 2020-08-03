To compile for TeslaP100

nvcc  -gencode arch=compute_70,code=sm_70 -c cufft3D.cu
pgf95 -acc -c -Minfo=accel -ta=tesla:cc70 -lm PSMHD3-cufftacc.f95 
pgf95 -acc -Minfo=accel -ta=tesla:cc70 -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm

To compile:
module load pgi/2018/pgi/18.4
module load cuda/9.0
module unload gcc/7.1.0

nvcc  -gencode arch=compute_70,code=sm_70 -c cufft3D.cu
pgf95 -acc -Minfo=accel -ta=tesla:cc70 -c PSMHD3-cufftacc.f95
#pgf95 -acc -Minfo=accel -ta=tesla:cc70,managed -c PSMHD3-cufftacc.f95
pgf95 -acc -Minfo=accel -ta=tesla:cc70 -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm -o a.out
#pgf95 -acc -Minfo=accel -ta=tesla:cc70,managed -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm -o a.out
srun --exclusive --nodes=1 --time=01:00:00 --partition=gpu --reservation=iprplasma_15 --pty "bash"
To run:
./a.out
