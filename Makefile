# Configuration
SRC_DIR = src

# Set up compile options
CC = gcc
RM = rm -f
MV = mv
EXTRA = -Wall -g

# Define the include files
#INC = input.h 2d_array.h ccdc.h output.h utilities.h
INC = $(wildcard $(SRC_DIR)/*.h)
INCDIR = -I$(SRC_DIR) -I$(XML2INC) -I$(ESPAINC) -I$(GSL_SCI_INC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the source code and object files
#SRC = input.c 2d_array.c ccdc.c utilities.c misc.c
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:.c=.o)

# Define the object libraries
LIB = -L$(GSL_SCI_LIB) -lz -lpthread -lrt -lgsl -lgslcblas -lm

# Define the executable
EXE = ccdc

# Target for the executable
all: $(EXE)

ccdc: $(OBJ) $(INC)
	$(CC) $(NCFLAGS) -o ccdc $(OBJ) $(LIB)

install:
	cp $(EXE) $(BIN)
	cp ../scripts/* $(BIN)

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) $(INCDIR) -c $<

