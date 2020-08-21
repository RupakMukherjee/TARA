t0 = 100 
t1 = 10100
dt = 10
pz = '0.01'

ofl = open('movie', 'w')

ofl.write('''#### gnu.wave ####

set grid

se xlabel "x"
se ylabel "u(x)"

se xr [0:2*pi]
se yr [-1.1:1.1]

#p "./fort.100" u 2:3 w lp pt 7 lt 3 ps 1 notitle

#set terminal pngcairo size 1500,1500 enhanced font 'Verdana,9'\n\n''')

for i in range(t0, t1 + 1, dt):
    #ofl.write('''se output "mag_field_%i.png"\n'''%(i))
    ofl.write('''p "./fort.%i" u 2:3 w lp lc 1 lt 3 pt 6 ps 0.5 notitle\npause %s\n''' % (i, pz))

ofl.close()

# ffmpeg -framerate 20 -pattern_type glob -i 'mag_field_*.png' -c:v libx264 -pix_fmt yuv420p mag_field_movie.mp4


