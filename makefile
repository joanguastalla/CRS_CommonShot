# Prefix used to recognize rules in recipes
.RECIPEPREFIX= 

#Define C Compiler
CC=gcc

# Compiler flags:
# -g Add debuging information to the executable file
# -Wall turn on most of compiler warnings
CFLAGS=-Wall  

# Define Library paths to include
LFLAGS=-L/usr/lib/x86_64-linux-gnu 
OPENMP=-fopenmp
# Define libraries to be used
LIBS= -lm -lblas -llapack 


# Header files for SRC
HCUBIC=mathutils.h	seismicutils.h omp.h 
 
# Construction of .o files
OBJCUBIC=commonreceiver.o mathutils.o cubicinterpol.o  crs_commmonshot.o

# Joining file object
SRC=testcrs.o

#Targets not constructed as files
.PHONY=all 

# Name of executable file
MAIN=testcrs


all: $(MAIN)
	@echo Compiling executable $(MAIN)


$(MAIN): $(SRC) $(OBJCUBIC)  
	$(CC) $(OPENMP) -o $@ $^ $(LFLAGS) $(LIBS) 

$(SRC):	testcrs.c
	$(CC) $(CFLAGS) $(OPENMP) -c -o $@ $< 

%.o: %.c $(HCUBIC)
	$(CC) -o $@ $< $(CFLAGS) 

clean:
	rm -f *.o
