!
! Rupak, IPR, August 10, 2017
!
! PSMHD3 is an OPENMP Parallel Benchmarked Three Dimensional Compressible Viscous Resistive MHD code, 
! using Pseudo-Spectral Method with Multiple Time Solvers.
! Here Adams-Bashforth Subroutine is Parallalised. Other Time Solvers are NOT Parallalised.
!
! Wave-form /sim exp(i(kx-wt)). Hence derivative w.r.t. x gives ik only (and NOT -ik).
!
!___________________________________________________________________________________________________________________________________________
!
! To Run this code in udbhav.ipr.res.in: 
!
!#!/bin/bash
!#Give your job name
!#BSUB -J MHD3D
!#Give required number of processor
!#BSUB -n 15
!#Give your email id
!##BSUB -u <rupak@ipr.res.in>
!#Give your current working directory
!##BSUB -cwd ~/<your home area>
!#Redirect error to some file and %J is Job ID
!#BSUB -e err.%J
!#Redirect output to some file and %J is Job ID
!#BSUB -o out.%J
!#Define number of processor per host user wants to execute a job
!#BSUB -R "span[hosts=1]"
!gfortran -fopenmp -I/home/rupak/bin/include 3dmhd.f95 -lfftw3 -lm
!./a.out
!
! Terminal Command: bsub < runjob.sh
!___________________________________________________________________________________________________________________________________________
!
! To Run this code in uday.ipr.res.in:
!
!#!/bin/bash
!# declare a name for this job to be sample_job
!#PBS -N MHD3D
!# request the parallel queue for this job
!##PBS -q parallel
!# request a total of 25 processors for this job (1 node and 25 processors per node)
!#PBS -l nodes=1:ppn=25
!# Request to run for specified hrs:min:secs
!##PBS -l walltime=500:00:00
!# combine PBS standard output and error files
!#PBS -j oe
!# mail is sent to you when the job starts and when it terminates or aborts
!##PBS -m bea
!# specify your email address
!#PBS -M rupak@gmail.com
!#change to the directory where you submitted the job
!cd $PBS_O_WORKDIR
!#include the full path to the name of your OMP program
!gfortran -fopenmp -I/soft/fftw-3.3.3/include -L/soft/fftw-3.3.3/lib 3dmhd.f95 -lfftw3 -lm
!./a.out
!
! Terminal Command: qsub jobfile.sh
!___________________________________________________________________________________________________________________________________________
!
! To Run this code in Desktop: 
!
!gfortran -fopenmp -I/usr/local/include 2dmhd.f95 -lfftw3 -lm
!./a.out
!___________________________________________________________________________________________________________________________________________
!
!\documentclass[12pt]{article} 
!\usepackage{a4wide,color}
!\begin{document}
!
!MHD Equations in Conservative form:\\
!
!\begin{eqnarray*}
!&& \frac{\partial \rho}{\partial t} + \vec{\nabla} \cdot \left(\rho \vec{u}\right) = 0\\
!&& \frac{\partial (\rho \vec{u})}{\partial t} + \vec{\nabla} \cdot \left[ \rho \vec{u} \vec{u} + \left(P + \frac{B^2}{2}\right){\bf{I}} - \vec{B}\vec{B} \right] = \mu \nabla^2 \vec{u}\\
!&& \frac{\partial E}{\partial t} + \vec{\nabla} \cdot \left[ \left( E + P + \frac{B^2}{2} \right)\vec{u} - \vec{u}\cdot\left( \vec{B} \vec{B} \right)  - \eta \vec{B} \times \left(\vec{\nabla} \times \vec{B} \right) \right] = \mu \left(\vec{\nabla} \cdot \vec{u} \right)^2\\
!&& \frac{\partial \vec{B}}{\partial t} + \vec{\nabla} \cdot \left( \vec{u} \vec{B} - \vec{B} \vec{u}\right) = \eta \nabla^2 \vec{B}\\
!\end{eqnarray*}
!
!In Two Dimensions, the MHD Equations in Conservative form becomes:\\
!
!\begin{eqnarray*}
!&& \frac{\partial \rho}{\partial t} + \frac{\partial}{\partial x} (\rho u_x) + \frac{\partial}{\partial y} (\rho u_y) = 0\\
!&& ~\\
!&& ~\\
!&& \frac{\partial \left( \rho u_x \right)}{\partial t} + \frac{\partial}{\partial x} \left[ \rho u_x u_x + P + \frac{B_x^2+B_y^2}{2} - B_x B_x \right] + \frac{\partial}{\partial y} \left[ \rho u_x u_y - B_x B_y \right] = \mu \left( \frac{\partial^2 u_x}{\partial x^2} + \frac{\partial^2 u_x}{\partial y^2} \right)^2\\
!&& \frac{\partial \left( \rho u_y \right)}{\partial t} + \frac{\partial}{\partial x} \left[ \rho u_x u_y - B_x B_y \right] + \frac{\partial}{\partial y} \left[ \rho u_y u_y + P + \frac{B_x^2+B_y^2}{2} - B_y B_y \right] = \mu \left( \frac{\partial^2 u_y}{\partial x^2} + \frac{\partial^2 u_y}{\partial y^2} \right)^2\\
!&& ~\\
!&& ~\\
!&& \frac{\partial E}{\partial t} + \frac{\partial}{\partial x} \left[ \left( E + P + \frac{B_x^2+B_y^2}{2} \right)u_x - u_x B_x B_x - u_y B_x B_y - \eta B_y \left( \frac{\partial B_y}{\partial x} - \frac{\partial B_x}{\partial y} \right) \right] \\
!&& ~~~~~ + \frac{\partial}{\partial y} \left[ \left( E + P + \frac{B_x^2+B_y^2}{2} \right)u_y - u_x B_x B_y - u_y B_y B_y + \eta B_x \left( \frac{\partial B_y}{\partial x} - \frac{\partial B_x}{\partial y} \right) \right] = \mu \left( \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} \right)^2\\
!&& ~\\
!&& ~\\
!&& \frac{\partial B_x}{\partial t} + \frac{\partial}{\partial y} \left( u_y B_x - B_y u_x \right) = \eta \left( \frac{\partial^2 B_x}{\partial x^2} + \frac{\partial^2 B_x}{\partial y^2} \right)\\
!&& \frac{\partial B_y}{\partial t} + \frac{\partial}{\partial x} \left( u_x B_y - B_x u_y \right) = \eta \left( \frac{\partial^2 B_y}{\partial x^2} + \frac{\partial^2 B_y}{\partial y^2} \right)\\
!\end{eqnarray*}
!
!\newpage
!
!Pseudocode for 2DMHD Code.
!
!\begin{eqnarray*}
!&& [At~time~ (t)] \\
!&& \\
!\rho, \rho u_x, \rho u_y, E, B_x, B_y ~~~~~ \Leftarrow ~~~ && \hat{\rho}, \widehat{\rho u_x}, \widehat{\rho u_y}, \hat{E}, \hat{B}_x, \hat{B}_y \longrightarrow \hat{\omega} = i k_x \hat{u}_y - i k_y \hat{u}_x\\
!\Downarrow ~~~~~~~~~~~~~~~~~~~~~~&& ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \Downarrow\\
!\vec{\nabla}\cdot \vec{B} \longleftarrow TE,E_y,P, B^2 \longleftarrow u_x, u_y, P &&  \Downarrow IFFT ~~~~~~~~~~~~~~~~~~~~~~~ \omega \rightarrow C^0, C^1, C^2, C^3  \\
!&& \\
!&& \rho, \rho u_x, \rho u_y, E, B_x, B_y \\
!&& \downarrow\\
!&& u_x = \frac{\rho u_x}{\rho}, u_y = \frac{\rho u_y}{\rho}\\
!&& B^2 = B_x^2 + B_y^2\\
!&& \\
!&& \Downarrow FFT\\
!&& \\
!&& \hat{u}_x, \widehat{u}_y\\
!&& \downarrow\\
!&& i k_x \hat{u}_x, i k_y \hat{u}_y\\
!&& i k_y \hat{B}_x, i k_x \hat{B}_y\\
!&& \\
!&& \Downarrow IFFT\\
!&& \\
!&& \frac{du_x}{dx}, \frac{du_y}{dy}\\
!&& \frac{dB_x}{dy}, \frac{dB_y}{dx}\\
!&& \downarrow\\
!&& CurlB = \frac{dB_y}{dx} - \frac{dB_x}{dy}\\
!&& P = \left(\gamma-1\right)\left[E - \frac{1}{2}\left(\rho u_x u_x + \rho u_y u_y + B^2\right)\right]\\
!&& \downarrow\\
!\end{eqnarray*}
!\newpage
!\begin{eqnarray*}
!&& \downarrow\\
!&& Mom_x^1 = \rho u_x u_x + P + \frac{B^2}{2} - B_x B_x\\
!&& Mom_x^2 = \rho u_x u_y - B_x B_y\\
!&& Mom_y^1 = \rho u_x u_y - B_x B_y\\
!&& Mom_y^2 = \rho u_y u_y + P + \frac{B^2}{2} - B_y B_y\\
!&& \downarrow\\
!&& Energy_x = \left(E+P+\frac{B^2}{2}\right)u_x - u_x B_x B_x - u_y B_x B_y - \eta ~ B_y ~ CurlB\\
!&& Energy_y = \left(E+P+\frac{B^2}{2}\right)u_y - u_x B_x B_y - u_y B_y B_y + \eta ~ B_x ~ CurlB\\
!&& \downarrow\\
!&& E_{Visc} = \left(\frac{du_x}{dx}+\frac{du_y}{dy}\right)^2\\
!&& \downarrow\\
!&& Mag_x = u_y B_x - B_y u_x\\
!&& Mag_y = u_x B_y - B_x u_y\\
!&& \\
!&& \Downarrow FFT\\
!&& \\
!&& i k_x \widehat{\rho u_x}, i k_y \widehat{\rho u_y}\\
!&& i k_x \widehat{Mom}_x^1, i k_y \widehat{Mom}_x^2, i k_x \widehat{Mom}_y^1, i k_y \widehat{Mom}_y^2\\
!&& k_x^2 \hat{u}_x, k_y^2 \hat{u}_x, k_x^2 \hat{u}_y, k_y^2 \hat{u}_y\\
!&& i k_x \widehat{Energy}_x, i k_y \widehat{Energy}_y\\
!&& i k_y \widehat{Mag}_x, i k_x \widehat{Mag}_y\\
!&& k_x^2 \hat{B}_x, k_y^2 \hat{B}_x, k_x^2 \hat{B}_y, k_y^2 \hat{B}_y\\
!&& \downarrow\\
!\end{eqnarray*}
!\newpage
!\begin{eqnarray*}
!&& \downarrow\\
!&& \frac{d\hat{\rho}}{dt} = - \left(i k_x \widehat{\rho u_x} + i k_y \widehat{\rho u_y}\right)\\
!&& \\
!&& \frac{d\widehat{\rho u_x}}{dt} = - \left(i k_x \widehat{Mom}_x^1 + i k_y \widehat{Mom}_x^2 \right) - \mu \left(k_x^2 \hat{u}_x + k_y^2 \hat{u}_x \right)\\
!&& \frac{d\widehat{\rho u_y}}{dt} = - \left(i k_x \widehat{Mom}_y^1 + i k_y \widehat{Mom}_y^2 \right) - \mu \left(k_x^2 \hat{u}_y + k_y^2 \hat{u}_y \right)\\
!&& \\
!&& \frac{d\hat{E}}{dt} = - \left(i k_x \widehat{Energy}_x + i k_y \widehat{Energy}_y\right) + \mu \widehat{E_{Visc}}\\
!&& \\
!&& \frac{d\hat{B}_x}{dt} = - i k_y \widehat{Mag}_x - \eta \left(k_x^2 \hat{B}_x + k_y^2 \hat{B}_x\right)\\
!&& \frac{d\hat{B}_y}{dt} = - i k_x \widehat{Mag}_y - \eta \left(k_x^2 \hat{B}_y + k_y^2 \hat{B}_y\right)\\
!&& \downarrow\\
!&& Adams ~ Bashforth\\
!&& \downarrow\\
!&& \hat{\rho}, \widehat{\rho u_x}, \widehat{\rho u_y}, \hat{E}, \hat{B}_x, \hat{B}_y\\
!&& \\
!&& [At~time~ (t+dt)]\\
!\end{eqnarray*}
!
!\end{document}
!___________________________________________________________________________________________________________________________________________
!
! 
