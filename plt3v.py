import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


for i in range (354100,525101,1000):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #plt.axis('off')
    ax.grid(False)
    ax.set_xlim3d(0, 6.3)
    ax.set_ylim3d(0, 6.3)
    ax.set_zlim3d(0, 6.3)
    ax.view_init(elev=20, azim=30)
    
    
    d0 = np.loadtxt("fort."+str(i), unpack = True)
    d1 = np.loadtxt("fort."+str(i+20), unpack = True)
    d2 = np.loadtxt("fort."+str(i+30), unpack = True)
    d3 = np.loadtxt("fort."+str(i+40), unpack = True)

    ax.scatter(d0[0][:100], d0[1][:100], d0[2][:100], c = 'k', s = 10, lw = 0)
    ax.scatter(d1[0], d1[1], np.zeros(len(d1[0])), c = d1[2], lw = 0)
    ax.scatter(np.zeros(len(d2[0])), d2[0], d2[1], c = d2[2], lw = 0)
    p4 = ax.scatter(d3[1], np.zeros(len(d3[0])), d3[0], c = d3[2], lw = 0)
    plt.colorbar(p4)
    plt.savefig('Tracer_v_'+str(i)+'.png', format='png')
    #plt.show()
    #raw_input('')
    plt.close()




# ffmpeg -framerate 20 -pattern_type glob -i 'Tracer_v_*.png' -c:v libx264 -pix_fmt yuv420p Tracer_v_movie.mp4
    
