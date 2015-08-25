# Makefile to compile mex/standalone version to Andy Liaw et al.'s C code (used in R package randomForest)
# Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
# License: GPLv2
# Version: 0.02


#  Makefile to generate mex or standalone. will work in cygwin (for windows) or linux
#
#  make mex: generates matlab mex files which can be easily called up
#  make diabetes: generates a standalone file to test on the pima indian
#                 diabetes dataset.
#


#source directory
SRC_FILES=classRF.c classTree.c rfutils.c cokus.c utilities.c classification.c get_args.c

CC=gcc
FORTRAN=gfortran # or g77 whichever is present
INCDIR  = -I. -I$(MATIO_INC) -I$(HDF5INC)
CFLAGS= -fpic -O2 -funroll-loops -march=native  -Wall $(INCDIR)
FFLAGS=-O2 -fpic -march=native#-g -Wall
LDFORTRAN=#-gfortran
MEXFLAGS=-O
INC = classification.h utilities.h 
LIB = -lgfortran -lm -DmxCalloc=calloc -DmxFree=free -L$(MATIO_LIB) -L$(HDF5LIB) -lmatio -lhdf5 -lhdf5_hl 

all:	clean rfsub classification

classification:  clean rfsub 
	$(CC) $(CFLAGS) $(INC) $(SRC_FILES) rfsub.o -o classification $(LIB)

rfsub:	$(SRC)rfsub.f
	echo 'Compiling rfsub.f (fortran subroutines)'
	$(FORTRAN)  $(FFLAGS) -c rfsub.f -o rfsub.o

install:
	cp $(EXE) $(BIN)
	cp ../scripts/* $(BIN)

clean:	
	rm classification -rf
	rm *~ -rf
	rm *.o -rf
