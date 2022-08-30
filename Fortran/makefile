#
# Try to determine the operating system:
#
ifeq '$(findstring ;,$(PATH))' ';'
    UNAME := Windows
else
    UNAME := $(shell uname 2>/dev/null || echo Unknown)
    ifeq ($(UNAME), Unknown) 
        UNAME := Windows # Assuming windows if we cannot find a uname
        $(warning    Cannot determine system; assuming Windows)
    endif
endif

$(info    System is $(UNAME))

# Check for other windows name
ifeq ($(findstring NT-,$(UNAME)),NT-)
    # Found
    UNAME=Windows
endif
#
# Define lib-names depending on the system:
#
ifeq ($(UNAME), Darwin)
 libMatlab = ../lib/libFEISTY.dylib
 libR = ../lib/libFEISTY.dylib
endif

ifeq ($(UNAME), Linux)
 libMatlab = ../lib/libFEISTY_linux_matlab.so
 libR = ../lib/libFEISTY.so
endif

ifeq ($(UNAME), Windows)
 libMatlab = ../lib/libFEISTY.dll
 libR = ../lib/libFEISTY.dll
endif

#
# iFortran:
#
#FORT=ifort
#FLAGS=-fPIC#-O2
#LINK= #-l
#
# gfortran:
#
FORT=gfortran
#FLAGS=-g3
FLAGS= -fPIC -O2
#LINK=
#
# nag fortran
#
#FORT = nagfor



OBJ = input.o globals.o spectrum.o fish.o FEISTY.o FEISTYtest.o FEISTY_wrap_R.o FEISTY_wrap_MATLAB.o

all: FEISTY libMatlab libR

debugging: FLAGS=-fPIC -g3 -fbacktrace#-debug extended # -check all -fpe0 -warn -traceback 

debugging:  all

%.o: %.f90
	$(FORT) -c  $(FLAGS) -o $@ $<

FEISTY: $(OBJ)  FEISTYtest.o
	$(FORT) $(LINK) -fPIC $(FLAGS) -o $@ $^

libMatlab: $(OBJ) FEISTY_wrap_MATLAB.o
	$(FORT) -shared $(LINK) $(FLAGS) -o $(libMatlab) $^

libR: $(OBJ) FEISTY_wrap_R.o
	$(FORT) -shared $(LINK) $(FLAGS) -o $(libR) $^

clean:
	rm -f *.o FEISTY $(libMatlab) $(libR)
