CC = mpicc
FC = mpif90
CFLAGS = -O3 -g -Wall 
FFLAGS = -O3 -Wall -cpp -ffree-line-length-none 

BIN = .
LIB = -lm -lmpi -lfftw3
INC = -I.
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)


all: clean main

main:	$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN)/main $(OBJ) $(LIB)

%.o: %.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) main *.mod


