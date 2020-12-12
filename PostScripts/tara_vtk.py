"""
** Transforming PINC output data from HDF5 to VTK format **
Requires: PyEVTK (Python Script to Export VTK)
		  EVTK (Export VTK) package allows exporting data to binary VTK files for visualization
		  and data analysis
@file                HDF52VTK.py
@author              Sayan Adhikari <sayan.adhikari@fys.uio.no>
Instruction for installation:
=============================
git clone https://github.com/paulo-herrera/PyEVTK.git
python3 setup.py install
=============================
Read more about transformation:
https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
https://github.com/paulo-herrera/PyEVTK
"""
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from evtk.hl import gridToVTK
from time import sleep
from tqdm import tqdm

Dir ="data"
outDir = "datavtk"
outFile = "tara"

tmin   = 110
tmax   = 5110
datadt = 1000
#Nt   = round((tmax-tmin)/dt)

Nx = 64
Ny = 64
Nz = 64

Lx = 2 * math.pi
Ly = 2 * math.pi
Lz = 2 * math.pi

x_vals = np.linspace(0, Lx, Nx, dtype='float64')
y_vals = np.linspace(0, Ly, Ny, dtype='float64')
z_vals = np.linspace(0, Lz, Nz, dtype='float64')

#========== TARA Data Directory Setup =============
if os.path.exists(outDir):
    os.system('rm -rf '+outDir)
    os.system('mkdir '+outDir)
else:
    os.system('mkdir '+outDir)
############################################
counter = 0
data_num = np.arange(start=tmin, stop=tmax, step=datadt, dtype=int)
for i in tqdm(range(len(data_num))):
	sleep(0.01)
	file=Dir+'/fort.%d'%data_num[i]
	ux,uy,uz,Bx,By,Bz = np.loadtxt(file, unpack=True)

	uxVTKFlag = ux.reshape(Nx,Ny,Nz)
	uyVTKFlag = uy.reshape(Nx,Ny,Nz)
	uzVTKFlag = uz.reshape(Nx,Ny,Nz)
	BxVTKFlag = Bx.reshape(Nx,Ny,Nz)
	ByVTKFlag = By.reshape(Nx,Ny,Nz)
	BzVTKFlag = Bz.reshape(Nx,Ny,Nz)

	uxVTK = np.require(uxVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	uyVTK = np.require(uyVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	uzVTK = np.require(uzVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	BxVTK = np.require(BxVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	ByVTK = np.require(ByVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	BzVTK = np.require(BzVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])

	gridToVTK("./"+outDir+"/"+outFile+"_%d"%counter, x_vals, y_vals, z_vals, pointData = {"ux" : uxVTK, "uy" : uyVTK, "uz" : uzVTK, "Bx" : BxVTK, "By" : ByVTK, "Bz" : BzVTK})

	counter = counter + 1
