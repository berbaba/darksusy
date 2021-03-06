#FF=$(FC)
FF=g77

### Compiler options ###

# Options for linux
FOPT = -funroll-loops -fno-emulate-complex -ffast-math -O -ffixed-line-length-none

### Setups for the DarkSUSY install directory ###

# Determine where to install stuff (prefix is set in configure)
prefix=/sps/km3net/users/ctoennis/Macros/WimpSim/v3.01/darksusy-5.0.6
# DS_INSTALL is where the library and data files will be installed
DS_INSTALL=${prefix}

LIB=$(DS_INSTALL)/lib
INC=$(DS_INSTALL)/include
cfitsio=.

### you must set GALPROP_LIBS if you have compiled galprop

all: dsmain dstest dstest-isasugra dstest-msugra-scan dstest-mssm7-scan dstest-mssm25-scan

dsmain : dsmain.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dsmain dsmain.f \
        $(LIB)/aldata.o -ldarksusy -lFH -lHB

dstest : dstest.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest dstest.f \
        -ldarksusy -lFH -lHB

dstest-galprop : dstest-galprop.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-galprop \
	dstest-galprop.f -ldarksusy -lstdc++ -L${cfitsio} -lcfitsio -lFH -lHB

dstest-galprop-one : dstest-galprop-one.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-galprop-one \
	dstest-galprop-one.f -ldarksusy -lstdc++ -L${cfitsio} -lcfitsio -lFH -lHB

dstest-isasugra : dstest-isasugra.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-isasugra \
        dstest-isasugra.f $(LIB)/aldata.o \
        -ldarksusy -lFH -lHB

dstest-msugra-scan : dstest-msugra-scan.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-msugra-scan \
        dstest-msugra-scan.f $(LIB)/aldata.o \
        -ldarksusy -lFH -lHB
	
dstest-msugra-fast-scan : dstest-msugra-fast-scan.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-msugra-fast-scan \
        dstest-msugra-fast-scan.f $(LIB)/aldata.o \
        -ldarksusy -lFH -lHB		
	
dstest-mssm7-scan : dstest-mssm7-scan.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-mssm7-scan \
        dstest-mssm7-scan.f $(LIB)/aldata.o \
        -ldarksusy -lFH -lHB	
	
dstest-mssm7-fast-scan : dstest-mssm7-fast-scan.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-mssm7-fast-scan \
        dstest-mssm7-fast-scan.f $(LIB)/aldata.o \
        -ldarksusy -lFH -lHB		
	
dstest-mssm25-scan : dstest-mssm25-scan.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-mssm25-scan \
        dstest-mssm25-scan.f $(LIB)/aldata.o \
        -ldarksusy -lFH -lHB		
	
dstest-mssm25-fast-scan : dstest-mssm25-fast-scan.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dstest-mssm25-fast-scan \
        dstest-mssm25-fast-scan.f $(LIB)/aldata.o \
        -ldarksusy -lFH -lHB		

dsslha2slha : dsslha2slha.f $(LIB)/libdarksusy.a
	$(FF) $(FOPT) -I$(INC) -L$(LIB) -o dsslha2slha dsslha2slha.f \
        $(LIB)/aldata.o -ldarksusy -lFH -lHB
