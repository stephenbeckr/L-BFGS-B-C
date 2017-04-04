# Feb 17 2015, converted fortran to C
# Stephen Becker, stephen.becker@colorado.edu


FC = gfortran
#CC = gcc
CC = g++
#CC = gcc-4.9

#FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 
#FFLAGS += -fdefault-integer-8
CFLAGS = -O3 -Wall -Wno-uninitialized 
#CFLAGS = -g -Wall -lm -DDEBUG -Wuninitialized -Wextra -Wno-unused-parameter
POSTFLAGS=-lm

# my f2c.h is here:
INCLUDES = -I/usr/include -I. -I/usr/local/include
#CFLAGS += $(INCLUDES)

DRIVER1_77 = driver1.c
DRIVER2_77 = driver2.c
DRIVER3_77 = driver3.c

LBFGSB  = lbfgsb.c linesearch.c subalgorithms.c print.c

LINPACK = linpack.c

BLAS 	= miniCBLAS.c
#POSTFLAGS += -D_USE_OPTIMIZED_BLAS -lblas

TIMER   = timer.c

SRC = $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(F2CFILES)

default: all test_1
all :  lbfgsb_77_1 lbfgsb_77_2 lbfgsb_77_3

lbfgsb_77_1 : $(DRIVER1_77) $(SRC)
	$(CC) -o x.lbfgsb_77_1 $(CFLAGS) $(DRIVER1_77) $(SRC) $(POSTFLAGS)

lbfgsb_77_2 : $(DRIVER2_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) -o x.lbfgsb_77_2 $(CFLAGS) $(DRIVER2_77) $(SRC) $(POSTFLAGS)

lbfgsb_77_3 : $(DRIVER2_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) -o x.lbfgsb_77_3 $(CFLAGS) $(DRIVER3_77) $(SRC) $(POSTFLAGS)

lbfgsb_77_1_nonstatic : $(SRC)
	$(CC) -o x.lbfgsb_77_1_nonstatic $(CFLAGS) driver_nonstatic.c $(SRC) $(POSTFLAGS)

test_1 : x.lbfgsb_77_1
	./x.lbfgsb_77_1
