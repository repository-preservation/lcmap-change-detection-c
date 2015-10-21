BIN ?= ../bin
SCRIPTS = ./scripts
EXE = classification

SRC_FILES = classRF.c classTree.c rfutils.c cokus.c utilities.c classification.c get_args.c

CC = gcc
FORTRAN = gfortran # or g77 whichever is present
HDF5INC ?= /usr/include/hdf5
HDF5LIB ?= /usr/lib/x86_64-linux-gnu/hdf5/serial
MATIOLIB ?= /usr/lib/x86_64-linux-gnu
INCDIR = -I. -I$(MATIO_INC) -I$(HDF5INC)
CFLAGS = -fpic -O2 -funroll-loops -march=native  -Wall $(INCDIR)
FFLAGS = -O2 -fpic -march=native#-g -Wall
LDFORTRAN = #-gfortran
MEXFLAGS = -O
INC = classification.h utilities.h rf.h
LIB = -lgfortran -lm -DmxCalloc=calloc -DmxFree=free -L$(MATIOLIB) -L$(HDF5LIB) -lmatio -lhdf5 -lhdf5_hl

all: clean rfsub $(EXE)

$(EXE):  clean rfsub
	$(CC) $(CFLAGS) $(INC) $(SRC_FILES) rfsub.o -o $(EXE) $(LIB)

rfsub: $(SRC)rfsub.f
	echo 'Compiling rfsub.f (fortran subroutines)'
	$(FORTRAN)  $(FFLAGS) -c rfsub.f -o rfsub.o

install:
	mv $(EXE) $(BIN)

clean:
	rm $(BIN)/$(EXE) -rf
	rm *~ -rf
	rm *.o -rf
