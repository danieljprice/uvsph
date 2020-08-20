FC=gfortran
FFLAGS=-O3
#
# the following files come from splash
#
SPLASH_FILES=asciiutils.f90 timing.f90 kernels.f90 interpolation.f90 interpolate2D.f90
#
# here are the files unique to uvsph
#
SRC=${SPLASH_FILES} read_uv.f90 uvsph.f90
OBJ=${SRC:.f90=.o}
#
# set directory to find splash source files
# by default we assume splash exists as a subdirectory
#
#ifndef SPLASH_DIR
SPLASH_DIR=./splash
#endif

VPATH=${SPLASH_DIR}/src

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

uvsph: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod
