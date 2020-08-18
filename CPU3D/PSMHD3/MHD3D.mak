#
#include /home/rupak/SourcePath.txt
F95 = /usr/bin/gfortran
SrcPath = /home/rupak/Untitled_Folder/FFTW/PSMHD3

# Definition of the Flags
OMP = -O -fopenmp
FFTWI = -I/usr/include
FFTWL = -lfftw3 -lm
#
##########################
# Object Files for build #
##########################

OBJS = \
PSMHD3.o \
derive.o \
ab.o \

PSMHD3 : $(OBJS)
	 ${F95} $(OMP) $(FFTWI) -o $@ $(OBJS) $(FFTWL)

#######################################
# Object dependencies and compilation #
#######################################
PSMHD3.o : $(SrcPath)/Files/main.f95 \
derive.o \
ab.o
	$(F95) -c $(OMP) $(FFTWI) -o $@ $(SrcPath)/Files/main.f95 $(FFTWL)

derive.o : $(SrcPath)/Files/derive.f95
	$(F95) -c $(OMP) $(FFTWI) -o $@ $(SrcPath)/Files/derive.f95 $(FFTWL)
	
ab.o : $(SrcPath)/Files/Time_Solver/ab.f95
	$(F95) -c $(OMP) $(FFTWI) -o $@ $(SrcPath)/Files/Time_Solver/ab.f95 $(FFTWL)

