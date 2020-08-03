GMHD3D (GPU MHD 3D)
===============================================

An GPU Parallel (NVIDIA CUDA) Benchmarked Three Dimensional Compressible Viscous Resistive MHD code, using Pseudo-Spectral Method with Multiple Time Solvers. Here Adams-Bashforth Subroutine is Parallalised. Other Time Solvers are NOT Parallalised.


Instructions
------------
### Prerequisites
- [Nvidia CUDA Compiler (NVCC)](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html)
- [PGI Community Edition](https://www.pgroup.com/products/community.htm)
- General Purpose Compilers (GCC/GFORTRAN)

### Compilation on [Nvidia Tesla P100](https://www.nvidia.com/en-us/data-center/tesla-p100/)
```shell
nvcc  -gencode arch=compute_70,code=sm_70 -c cufft3D.cu
pgf95 -acc -c -Minfo=accel -ta=tesla:cc70 -lm PSMHD3-cufftacc.f95 
pgf95 -acc -Minfo=accel -ta=tesla:cc70 -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm
```
### Compilation on [Nvidia Quadro K6000](https://www.pny.com/nvidia-quadro-k6000)
```shell
nvcc  -gencode arch=compute_70,code=sm_70 -c cufft3D.cu
pgf95 -acc -c -Minfo=accel -ta=tesla:cc35 -lm PSMHD3-cufftacc.f95 
pgf95 -acc -Minfo=accel -ta=tesla:cc35 -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm
```
### Compilation on HPC [Nvidia](#)
```shell
module load pgi/2018/pgi/18.4
module load cuda/9.0
module unload gcc/7.1.0

nvcc  -gencode arch=compute_70,code=sm_70 -c cufft3D.cu
pgf95 -acc -Minfo=accel -ta=tesla:cc70 -c PSMHD3-cufftacc.f95
#pgf95 -acc -Minfo=accel -ta=tesla:cc70,managed -c PSMHD3-cufftacc.f95
pgf95 -acc -Minfo=accel -ta=tesla:cc70 -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm -o a.out
#pgf95 -acc -Minfo=accel -ta=tesla:cc70,managed -Mcudalib=cufft PSMHD3-cufftacc.o cufft3D.o -lm -o a.out
srun --exclusive --nodes=1 --time=01:00:00 --partition=gpu --reservation=iprplasma_15 --pty "bash"
```
### Run
```shell
./a.out
```
### GPU Status
```shell
nvidia-smi
```
