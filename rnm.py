import glob
import shutil
fs = glob.glob("*.png")
for f in fs:
    if len(f) == 13:
        shutil.copyfile(f, './new/' + f[:6] + '0000' + f[6:])
    elif len(f) == 14:
        shutil.copyfile(f, './new/' + f[:6] + '000' + f[6:])
    elif len(f) == 15:
        shutil.copyfile(f, './new/' + f[:6] + '00' + f[6:])
    elif len(f) == 16:
        shutil.copyfile(f, './new/' + f[:6] + '0' + f[6:])        
    elif len(f) == 17:
        shutil.copyfile(f, './new/' + f)
        
        
# ffmpeg -framerate 20 -pattern_type glob -i 'Tracer_v_*.png' -c:v libx264 -pix_fmt yuv420p Tracer_v_movie.mp4

# ffmpeg -framerate 20 -pattern_type glob -i '*.png' -vcodec libx264 -y -an Iso_B_movie.mp4 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"

