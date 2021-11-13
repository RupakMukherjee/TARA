# TARA

[![build](https://github.com/RupakMukherjee/TARA/actions/workflows/main.yml/badge.svg)](https://github.com/RupakMukherjee/TARA/actions/workflows/main.yml)
[![DOI](https://zenodo.org/badge/136785629.svg)](https://zenodo.org/badge/latestdoi/136785629)

Multi dimensional Navier-Stokes equation solver using pseudo-spectral and WENO algorithms for spatial discretization and explicit time iterators.

A detailed documentation can be found at https://rupakmukherjee.github.io/TARA/

[Example Movie-1](https://youtu.be/YwlT4K7pnMs) [Example Movie-2](https://youtu.be/E1fPny0DuRo) [Example Movie-3](https://youtu.be/E1lOt9nsKMk) [Example Movie-4](https://youtu.be/k6Lg6q7t_Kc) 

Installation
------------
#### Prerequisites
1. [GNU Make](https://www.gnu.org/software/make/)
2. [FFTW Library](http://www.fftw.org/) (Source already included!)
3. [python3 or higher](https://www.python.org/download/releases/3.0/)
4. [git](https://git-scm.com/)

#### Procedure
First make a clone of the master branch using the following command
```shell
git clone https://github.com/RupakMukherjee/TARA.git
```
Then enter inside the *TARA* directory 
```shell
cd TARA
```
Now complile and build the *TARA* code
```shell
make subsystems
``` 

```shell
make all
``` 
Usage
-----
Upon successful compilation, run the code using following command
```shell
make run
```
## Parameter Setup
Edit the _input.ini_ and run the code again. The basic structure of _input.ini_ is provided below,

```ini
[dimension]
#Options: 1, 2, 3
dim = 2

[architecture]
#Options: cpu, single-gpu, multi-gpu
arch = cpu

[parallel]
#Options: serial, openmp, mpi
mode = serial

[type]
#Options: hydro, mhd
domain = hydro

[problem]
#Options: bare, screen
hydrotype = screen

[particle]
Np = 10

[grid]
Nx = 64
Ny = 64
Nz = 64

[length]
#All lengths will be multiplied by 2\pi.
Lx = 2.0d0
Ly = 2.0d0
Lz = 2.0d0

[resistivity]
nu = 0.010d0
eta = 0.0d0

[time]
time_min = 0.0d0
time_max = 2.0d0
dt = 0.00010d0

[initialProfile]
W0 = 25.0d0
m0 = 2.0d0
d = 3.0d0/128.0d0

[initialB]
B0 = 0.0d0
```

## Contributing
We welcome contributions to this project.

1. Fork it.
2. Create your feature branch (```git checkout -b my-new-feature```).
3. Commit your changes (```git commit -am 'Add some feature'```).
4. Push to the branch (```git push origin my-new-feature```).
5. Create new Pull Request.

## License
Released under the [MIT license](LICENSE).
