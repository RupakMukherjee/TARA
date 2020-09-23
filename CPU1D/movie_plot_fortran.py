#!/usr/bin/ python

import numpy as np
import pylab as plt
import random
import math
import matplotlib.animation as animation

#========= Configuration ===========
show_anim = True
interval = 0.1#in seconds
datadt = 20
# DIR ="weno_fortran"
Lx = 10.0
Ly = 10.0

pi = 3.14159265359;
tmax = 1.5
tmin = 0.0
dt   = 0.001
Nt   = round((tmax-tmin)/dt)
L = 2.0*pi

print(Nt)

data_num = np.arange(start=101, stop=Nt, step=datadt, dtype=int)

if (show_anim == True):
    def animate(i):
        # file=DIR+'/fort.%d'%data_num[i]
        file='fort.%d'%data_num[i]
        datax,datau = np.loadtxt(file, unpack=True)

        ax1.cla()
        img1 = ax1.plot(datax,datau,marker='.',color='b',alpha=1.0)
        ax1.set_title('TimeSteps = %d'%i+'\n Phase Space')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$u(x)$")
        ax1.set_xlim([0, L])
        # ax1.set_ylim([-Ly, Ly])



if (show_anim == True):
    fig,ax1 = plt.subplots(1,1,figsize=(6, 6))
    ani = animation.FuncAnimation(fig,animate,frames=len(data_num),interval=interval*1e+3,blit=False)
    plt.show()
