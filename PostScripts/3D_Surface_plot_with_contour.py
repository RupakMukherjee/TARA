# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This is the python script for the visualization of 3d surface along with the 2D plane projection as well as wired frame view of a 2d data set.
# Developed by: Mr. Shishir Biswas, IPR, India.
# Date: 26 th Sep 2020.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math
Nx=2048
Ny=2048
j=372
for i in range(1860000, 2005000, 5000):
 f=np.loadtxt('fort.'+ str(i))
 xlist=f[:,0]
 ylist=f[:,1]
 xi=np.linspace(np.min(xlist), np.max(xlist), Nx)
 yi=np.linspace(np.min(ylist), np.max(ylist), Ny)
 Y, X = np.meshgrid(yi, xi)
# plt.figure()
# plt.xlabel('x-axis')
# plt.ylabel('y-axis')
 z=f[:,4]
 zi=np.reshape(z,(Nx,Ny))
# plt.ylim([0,2*math.pi])
# plt.xlim([0,2*math.pi])
 levels=900
 
 fig = plt.figure()
 ax = fig.add_subplot(111,projection='3d') #This is needed for 3D plot. The 111 is the ratio of 3 axes. try to keep it 111 for cubic plot.
# ax = plt.axes(projection='3d') # Or you can also use this for making box
 plt.xlabel('X')
 plt.ylabel('Y')
 plt.ylim([0,2*math.pi])
 plt.xlim([0,2*math.pi])
 surf = ax.plot_surface(X, Y, zi,linewidths=0.5,cmap='nipy_spectral',rstride=1, cstride=1) #camp controls the color of the color bar and can be choosen like ['jet','nipy_spectral','gist_ncar','gist_rainbow', 'rainbow'] # [ rstride=1, cstride=1,antialiased=False ] these are only for decoration purpose one can remove these also.
 
# wireframe = ax.plot_wireframe(X, Y, zi,linewidths=0.2,color='black',alpha=0.1,rstride=2, cstride=2)
 
 cset = ax.contour(X, Y, zi, levels, linewidths=0.5, cmap='Greys',zdir='z',offset=np.min(zi)) # This is needed for the contour plot.

 cb = plt.colorbar(surf, shrink=0.7, aspect=12) # shrink aspect controls the dimension of colorbar
# cb = plt.colorbar(cset, shrink=0.7, aspect=12) # shrink aspect controls the dimension of colorbar
 cb.set_label(label='$\omega$', size=15, weight='bold', rotation=0)
# cb.set_label(label='$\omega$ Contour', size=15, weight='bold', rotation=0)
  
# plt.show()

 
 plt.title('Vorticity Evolution' "\n" '[Rn = 228576, $\phi_{m}$ = 0 ]' "\n" 'Time = '+str((i/1000.0)-5.0),color='indigo',fontsize=14, fontweight='bold')


 plt.savefig('w' + str(j) + '.png', format='png', dpi=300)
 j=j+1
 plt.clf()
 plt.close('all')
