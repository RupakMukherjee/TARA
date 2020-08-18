#
#include /home/rupak/SourcePath.txt
F95 = /usr/bin/gfortran
SrcPath = /home/rupak/udbhav/My_PhD/FFTW/3D_MHD/PSMHD3

# Definition of the Flags
OMP = -O -fopenmp
FFTWI = -I/home/rupak/bin/include
FFTWL = -lfftw3 -lm
#
##########################
# Object Files for build #
##########################

OBJS = \
PSMHD3.o \
parameter.o \
input.o \
vorticity.o \
derive.o \
u_f_b.o \
dv_db.o \
real_eval.o \
spectral_eval.o \
dealiaz.o \
nld.o \
ab.o \
reset.o \
div_b.o \
u_b2_p.o \
quantity.o \
e_b_spectra.o \

PSMHD3 : $(OBJS)
	 ${F95} $(OMP) $(FFTWI) -o $@ $(OBJS) $(FFTWL)

#######################################
# Object dependencies and compilation #
#######################################


# main.f95


PSMHD3.o : $(SrcPath)/Files/main.f95 \
parameter.o \
input.o \
vorticity.o \
derive.o \
ab.o \
reset.o \
div_b.o \
u_b2_p.o \
quantity.o \
e_b_spectra.o
	$(F95) -c $(OMP) $(FFTWI) -o $@ $(SrcPath)/Files/main.f95 $(FFTWL)

parameter.o : $(SrcPath)/Files/parameter.f95 
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/parameter.f95

input.o : $(SrcPath)/Files/input.f95
	$(F95) -c -o $@ $(SrcPath)/Files/input.f95

vorticity.o : $(SrcPath)/Files/vorticity.f95
	$(F95) -c $(OMP) $(FFTWI) -o $@ $(SrcPath)/Files/vorticity.f95 $(FFTWL)

reset.o : $(SrcPath)/Files/reset.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/reset.f95

div_b.o : $(SrcPath)/Files/div_b.f95
	$(F95) -c $(OMP) $(FFTWI) -o $@ $(SrcPath)/Files/div_b.f95 $(FFTWL)

u_b2_p.o : $(SrcPath)/Files/u_b2_p.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/u_b2_p.f95

quantity.o : $(SrcPath)/Files/quantity.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/quantity.f95

e_b_spectra.o : $(SrcPath)/Files/e_b_spectra.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/e_b_spectra.f95


# derive.f95


derive.o : $(SrcPath)/Files/N_L_Derivative/derive.f95 \
u_f_b.o \
dv_db.o \
real_eval.o \
spectral_eval.o \
dealiaz.o \
nld.o
	$(F95) -c $(OMP) $(FFTWI) -o $@ $(SrcPath)/Files/N_L_Derivative/derive.f95 $(FFTWL)
	
u_f_b.o : $(SrcPath)/Files/N_L_Derivative/u_f_b.f95	
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/N_L_Derivative/u_f_b.f95	

dv_db.o : $(SrcPath)/Files/N_L_Derivative/dv_db.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/N_L_Derivative/dv_db.f95

real_eval.o : $(SrcPath)/Files/N_L_Derivative/real_eval.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/N_L_Derivative/real_eval.f95

spectral_eval.o : $(SrcPath)/Files/N_L_Derivative/spectral_eval.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/N_L_Derivative/spectral_eval.f95

dealiaz.o : $(SrcPath)/Files/N_L_Derivative/de_aliaz.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/N_L_Derivative/de_aliaz.f95
	
nld.o : $(SrcPath)/Files/N_L_Derivative/nld.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/N_L_Derivative/nld.f95	
		
		
# ab.f95		
		
		
ab.o : $(SrcPath)/Files/Time_Solver/ab.f95
	$(F95) -c $(OMP) -o $@ $(SrcPath)/Files/Time_Solver/ab.f95

