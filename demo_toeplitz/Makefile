CC = gcc
CFLAGS = -O3 -g -Wall 

libcKrylov_INC = ../include
libcKrylov_LIB = ../lib


BIN = .
LIB = -L$(libcKrylov_LIB) -lcKrylov -lm 
INC = -I$(libcKrylov_INC) -I. 
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
	rm -f $(OBJ) main


