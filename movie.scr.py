t0 = 354100 
t1 = 525100
dt = 1000
pz = '0.05'

ofl = open('movie', 'w')

ofl.write('''#### gnu.wave ####

set grid

#se xlabel font ",25"
#se ylabel font ",25"
#se zlabel font ",25"
#set tics font ", 25"
#se label font ", 20"

se xlabel "X"
se ylabel "Y"
se zlabel "Z"

#se xr [0:2*pi]
#se yr [0:2*pi]
#se zr [0:2*pi]

#se view 80, 343

xmin = 0
xmax = 2*pi
ymin = 0
ymax = 2*pi
zmin = 0
zmax = 2*pi

#se object 1 polygon from xmin,ymax,zmin to xmin,ymax,zmax to xmax,ymax,zmax to xmax,ymax,zmin fs solid 0.5 fc lt 6 # Back
#se object 2 polygon from xmin,ymin,zmin to xmin,ymax,zmin to xmax,ymax,zmin to xmax,ymin,zmin fs solid 0.15 fc lt 7 # Bottom
#se object 3 polygon from xmin,ymin,zmin to xmin,ymax,zmin to xmin,ymax,zmax to xmin,ymin,zmax fs transparent solid 0.5 fc lt 1 # Left
#se object 4 polygon from xmax,ymin,zmin to xmax,ymax,zmin to xmax,ymax,zmax to xmax,ymin,zmax fs transparent solid 0.5 fc lt 1 # Right
#se object 5 polygon from xmin,ymin,zmax to xmin,ymax,zmax to xmax,ymax,zmax to xmax,ymin,zmax fs transparent solid 0.15 fc lt 7 # Top

sp "<(sed -n '4,4p' ./fort.100)" u 1:2:3 w p pt 7 lt 3 ps 0.1 notitle

#set terminal pngcairo size 1500,1500 enhanced font 'Verdana,9'\n\n''')

for i in range(t0, t1 + 1, dt):
    #ofl.write('''se output "Tracer_p_%i.png"\n'''%(i))
    ofl.write('''rep "<(sed -n '4,4p' ./fort.%i)" u 1:2:3 w p pt 7 lt 3 ps 0.1 notitle\npause %s\n''' % (i, pz))

ofl.close()

# ffmpeg -framerate 20 -pattern_type glob -i 'Tracer_p_*.png' -c:v libx264 -pix_fmt yuv420p Tracer_p_movie.mp4


