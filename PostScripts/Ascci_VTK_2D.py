"""
** Transforming TARA output data from ASCII to VTK format **
Requires: PyEVTK (Python Script to Export VTK)
		  EVTK (Export VTK) package allows exporting data to binary VTK files for visualization
		  and data analysis
@file                tara_vtk.py
@author              Sayan Adhikari <sayan.adhikari@fys.uio.no>
                     Rupak Mukherjee <rupakm@princeton.edu>
Instruction for installation:
=============================
git clone https://github.com/paulo-herrera/PyEVTK.git
python3 setup.py install
=============================
Read more about transformation:
https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
https://github.com/paulo-herrera/PyEVTK
"""

"""
Instruction for installation [In Shishir Desktop]:
=============================
git clone https://github.com/paulo-herrera/PyEVTK.git
python3 setup.py install --user
pip3 install tqdm
=============================
"""

"""
Instruction for running the script in Antya:
=============================
module load PyEVTK/EVTK 
export PYTHONPATH=/home/application/PyEVTK/lib/python2.7/site-packages/:$PYTHONPATH
** As there is some problem in the module installattion so we have to give the actual path of EVTK in terminal **
=============================
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
outFile = "GIHD2D"

tmin   = 5000
tmax   = 2005000
datadt = 5000
#Nt   = round((tmax-tmin)/dt)

Nx = 2048
Ny = 2048
Nz = 1



Lx = 2 * math.pi
Ly = 2 * math.pi
Lz = 0#2 * math.pi



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
	file='fort.%d'%data_num[i]
	x,y,ux,uy,omega,psi = np.loadtxt(file, unpack=True)
	
	

	uxVTKFlag = ux.reshape(Nx,Ny,Nz)
	uyVTKFlag = uy.reshape(Nx,Ny,Nz)
	omegaVTKFlag = omega.reshape(Nx,Ny,Nz)
	psiVTKFlag = psi.reshape(Nx,Ny,Nz)
	
	

	uxVTK = np.require(uxVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	uyVTK = np.require(uyVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	omegaVTK = np.require(omegaVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	psiVTK = np.require(psiVTKFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
	

	gridToVTK("./"+outDir+"/"+outFile+"_%d"%counter, x_vals, y_vals,  z_vals, pointData = {"ux" : uxVTK, "uy" : uyVTK, "omega" : omegaVTK, "psi" : psiVTK})

	counter = counter + 1

