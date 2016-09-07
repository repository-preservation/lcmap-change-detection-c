# Set up compile options
CC    = gcc
FORTRAN=gfortran 
RM    = rm -f
MV    = mv
CFLAGS = -Wabi -g
#EXTRA = -Wall -g

# Define the include files
INC = input.h 2d_array.h ccdc.h output.h utilities.h 
INCDIR  = -I. -I$(XML2INC) -I$(ESPAINC) -I$(GSL_SCI_INC) -I$(MATIO_INC) -I$(HDF5INC)
NCFLAGS = $(CFLAGS) $(INCDIR)
FFLAGS=-fdefault-real-8 -g

# Define the source code and object files
SRC = input.c 2d_array.c ccdc.c utilities.c misc.c multirobust.c
OBJ = $(SRC:.c=.o)

# Define the object libraries
LIB = -L$(GSL_SCI_LIB) -lz -lpthread -lrt -lgsl -lgslcblas -lgfortran -lm \
	-DmxCalloc=calloc -DmxFree=free -L$(MATIO_LIB) -L$(HDF5LIB) -lmatio -lhdf5 -lhdf5_hl

# Define the executable
EXE = ccdc

# Target for the executable
all: $(EXE)

ccdc: $(OBJ) glmnet5 $(INC)
	$(CC) $(NCFLAGS) -o ccdc $(OBJ) glmnet5.o $(LIB)

glmnet5: $(SRC) glmnet5.f
	$(FORTRAN) $(FFLAGS) -c glmnet5.f -o glmnet5.o

install:
	cp $(EXE) $(BIN)

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) $(INCDIR) -c $<

