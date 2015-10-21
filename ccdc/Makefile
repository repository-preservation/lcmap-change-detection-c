# Configuration
SRC_DIR = ./
SCRIPTS = ./scripts
BIN ?= ../bin
XML2INC ?= /usr/include/libxml2/libxml
ESPAINC ?=
GSL_SCI_INC ?= /usr/include/gsl
GSL_SCI_LIB ?= /usr/lib

# Set up compile options
CC = gcc
RM = rm -f
MV = mv
EXTRA = -Wall -Wextra -g

# Define the include files
INC = $(wildcard $(SRC_DIR)/*.h)
INCDIR  = -I$(SRC_DIR) -I$(GSL_SCI_INC) -I$(XML2INC) -I$(ESPAINC) -I$(GSL_SCI_INC)
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

$(BIN):
	mkdir -p $(BIN)

install: $(BIN)
	mv $(EXE) $(BIN)
	cp $(SCRIPTS)/* $(BIN)

clean:
	$(RM) $(BIN)/$(EXE)
	$(RM) $(BIN)/*.r
	$(RM) *.o

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) $(INCDIR) -c $<

