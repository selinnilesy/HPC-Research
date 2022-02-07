CC = g++

matrix:  matrixGenerator.o mmio.o
    $(CC) matrixGenerator.o mmio.o -o matrix.out

%.o: %.c %.h
    $(CC) -ansi -Wall  -c $<