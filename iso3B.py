from mayavi import mlab
import numpy as np


for i in range(354100,525101,1000):
    t, t, t, u ,v, w = np.loadtxt("fort."+str(i+10), unpack = True)
    u = np.reshape(u, (64, 64, 64))
    v = np.reshape(v, (64, 64, 64))
    w = np.reshape(w, (64, 64, 64))
    
    src = mlab.pipeline.vector_field(u, v, w)

    magnitude = mlab.pipeline.extract_vector_norm(src)
    mlab.pipeline.iso_surface(magnitude, contours=[0.13,0.16,0.2])
    
    mlab.savefig('Iso_B_'+str(i)+'.png')
    #mlab.show()
    mlab.close()




# ffmpeg -framerate 20 -pattern_type glob -i '*.png' -vcodec libx264 -y -an Iso_B_movie.mp4 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"
    
    
