# Configuration
SRC_DIR = src
BIN =? ./bin
XML2INC = /usr/include/libxml2/libxml
ESPAINC =
GSL_SCI_INC = /usr/include/gsl
GSL_SCI_LIB = /usr/lib

# Set up compile options
CC = gcc
RM = rm -f
MV = mv
EXTRA = -Wall -Wextra -g

# Define the include files
#INC = input.h 2d_array.h ccdc.h output.h utilities.h
#INC = $(SRC_DIR)/input.h $(SRC_DIR)/2d_array.h $(SRC_DIR)/ccdc.h $(SRC_DIR)/output.h $(SRC_DIR)/utilities.h
INC = $(wildcard $(SRC_DIR)/*.h)
#INCDIR = -I$(SRC_DIR) -I$(XML2INC) -I$(ESPAINC) -I$(GSL_SCI_INC)
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
	cp $(EXE) $(BIN)
	cp scripts/* $(BIN)

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) $(INCDIR) -c $<

