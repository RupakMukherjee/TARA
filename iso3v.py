from mayavi import mlab
import numpy as np


for i in range(424100,525101,1000):
    u, v, w, t ,t, t = np.loadtxt("fort."+str(i+10), unpack = True)
    u = np.reshape(u, (64, 64, 64))
    v = np.reshape(v, (64, 64, 64))
    w = np.reshape(w, (64, 64, 64))
    
    src = mlab.pipeline.vector_field(u, v, w)

    magnitude = mlab.pipeline.extract_vector_norm(src)
    mlab.pipeline.iso_surface(magnitude, contours=[0.05,0.01,0.001])
    
    mlab.savefig('Iso_v_'+str(i)+'.png')
    #mlab.show()
    mlab.close()




# ffmpeg -framerate 20 -pattern_type glob -i '*.png' -vcodec libx264 -y -an Iso_v_movie.mp4 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"
    
