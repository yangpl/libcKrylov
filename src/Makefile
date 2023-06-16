CC = gcc
CFLAGS = -O3 -g -Wall 

# variables for the libraries
AR = ar
ARFUNCT = cruv
AREXTRA = xo
RANLIB = ranlib
RM = rm -rf
#==================================================================

LIB = -lm -fopenmp 
INC = -I../include 
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: lib clean

lib : $(OBJ) 
	$(AR) $(ARFUNCT) ../lib/libcKrylov.a $(OBJ)

%.o : %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

clean :	
	$(RM) *.o  *~
