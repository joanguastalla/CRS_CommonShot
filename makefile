# Prefix used to recognize rules in recipes
.RECIPEPREFIX= 

#Define C Compiler
CC=gcc

# Compiler flags:
# -g Add debuging information to the executable file
# -Wall turn on most of compiler warnings
CFLAGS= -c -g -Wall -pedantic

# Define Library paths to include
LFLAGS=-L/usr/lib/x86_64-linux-gnu

# Define libraries to be used
LIBS= -lm -lblas -llapack

# Define C source files
SRC=cubicinterpol.c

# Header files for SRC
HCUBIC=mathutils.h lapack.h
 
# Construction of .o files
OBJCUBIC=cubicinterpol.o transpose.o sgn.o

#Targets not constructed as files
.PHONY=all 

# Name of executable file
MAIN=cubicinterpol


all: $(OBJCUBIC) $(MAIN)
	@echo Compiling executable $(MAIN)


$(MAIN): $(OBJCUBIC) 
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS)

%.o: %.c $(HCUBIC)
	$(CC) -o $@ $< $(CFLAGS)

clean:
	rm -f *.o
