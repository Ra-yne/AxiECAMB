#CAMB Makefile

#Set FISHER=Y to compile bispectrum fisher matrix code
FISHER=

#Edit for your compiler
#Note there are many old ifort versions, some of which behave oddly


#Intel , -openmp toggles mutli-processor:
#note version 10.0 gives wrong result for lensed when compiled with -openmp [fixed in 10.1]
#F90C     = ifort
#FFLAGS= -O3 -qopenmp #-g -warn all -check all # ##-O3  # #RL 08012023 -ip -WB  -fp-model precise -fp-model fast=2 
#(-ip means interprocedural optimizations across files. Potentially 
#FFLAGS = -O2 #3 -ffree-form -x f95-cpp-input
## This is flag is passed to the Fortran compiler allowing it to link C++ if required (not usually):
#F90CRLINK = -cxxlib
#ifneq ($(FISHER),)
#FFLAGS += -mkl
#endif

#Gfortran compiler:
#The options here work in v4.5, delete from RHS in earlier versions (15% slower)
#if pre v4.3 add -D__GFORTRAN__
#With v4.6+ try -Ofast -march=native -fopenmp
#A local machine test indicates that v4.5 can be about 20% slower than ifort
F90C     = gfortran 
#Below are suggested fflags for gfortran users that suppress harmless warnings at runtime but trap serious exceptions  
FFLAGS = -O3 -fopenmp -ffpe-summary=none


#Old Intel ifc, add -openmp for multi-processor (some have bugs):
#F90C     = ifc
#FFLAGS = -O2 -Vaxlib -ip -W0 -WB -quiet -fpp2
#some systems can can also add e.g. -tpp7 -xW

#G95 compiler
#F90C   = g95
#FFLAGS = -O2

#SGI, -mp toggles multi-processor. Use -O2 if -Ofast gives problems.
#F90C     = f90
#FFLAGS  = -Ofast -mp

#Digital/Compaq fortran, -omp toggles multi-processor
#F90C    = f90
#FFLAGS  = -omp -O4 -arch host -math_library fast -tune host -fpe1

#Absoft ProFortran, single processor:
#F90C     = f95
#FFLAGS = -O2 -cpu:athlon -s -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_"

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3

#PGF90
#F90C = pgf90
#FFLAGS = -O2 -DESCAPEBACKSLASH -Mpreprocess

#Sun V880
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign -DMPI

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O2 -xarch=native64 -openmp -ftrap=%none
#try removing -openmp if get bus errors. -03, -04 etc are dodgy.

#IBM XL Fortran, multi-processor (run gmake)
#F90C     = xlf90_r
#FFLAGS  = -DESCAPEBACKSLASH -DIBMXL -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3

#Settings for building camb_fits
#Location of FITSIO and name of library
#FITSDIR       ?= /home/cpac/cpac-tools/lib
#FITSLIB       = cfitsio
#Location of HEALPIX for building camb_fits
#HEALPIXDIR    ?= /home/cpac/cpac-tools/healpix

#ifneq ($(FISHER),)
#FFLAGS += -DFISHER
#EXTCAMBFILES = Matrix_utils.o
#else
#EXTCAMBFILES =
#endif

include ./Makefile_main
