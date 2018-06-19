# !! Define Environment variable FARGO_ARCH to reflect your architecture
# for ex.: setenv FARGO_ARCH LINUX (or export FARGO_ARCH=LINUX depending on your shell)
# Possible choices : undefined, LINUX, INTEL, MACOSX, EV68, EVMP, SPARC

# Generic Linux platform
# FARGO_ARCH must be set to LINUX
CC_LINUX  = gcc
OPT_LINUX = -O4 -march=native -Wall -ffast-math -Wno-unknown-pragmas
OPTSEQ_LINUX =
PARAOPT_LINUX =
PARACC_LINUX = mpicc

# AMD Opteron platform
# FARGO_ARCH must be set to OPTERON
CC_OPTERON = gcc
OPT_OPTERON = -O3 -march=opteron -Wall -ffast-math -Wno-unknown-pragmas
OPTSEQ_OPTERON =
PARAOPT_OPTERON =
PARACC_OPTERON = mpicc

# Intel Linux platform
# optimized for pentium 4 (otherwise edit accordingly)
# FARGO_ARCH must be set to INTEL
# (sequential only, MPI identical to generic IA32)
CC_INTEL  = icc
OPT_INTEL =
OPTSEQ_INTEL = -O3 -tpp7 -ipo -unroll -w2 -Wall -wd1418,310,266,810,1419,981,161
PARAOPT_INTEL = -O4 -march=pentium3 -Wall -ffast-math
PARACC_INTEL = mpicc

# Macintosh MAC OS X platform **POWERPC ONLY** (with gcc)
# FARGO_ARCH must be set to MACOSX
CC_MACOSX  = gcc
OPT_MACOSX = -Wall -O3 -ffast-math -mcpu=G4 -Wno-unknown-pragmas
OPTSEQ_MACOSX =
PARAOPT_MACOSX =
PARACC_MACOSX = mpicc

# Macintosh MAC OS X platform **Mac Intel  ONLY** (with gcc)
# FARGO_ARCH must be set to MacIntel
CC_MacIntel = gcc
#Below : OpenMP support; experimental feature: code not recently tested
#with OpenMP parallelization. Uncomment this line at your own risk.
#OPT_MacIntel = -O3 -Wall -arch i386 -fast -ftree-vectorize -fopenmp
OPT_MacIntel = -O3 -Wall -arch i386 -fast -Wno-unknown-pragmas
OPTSEQ_MacIntel =
PARAOPT_MacIntel =
PARACC_MacIntel = mpicc

# Generic options valid for any platfom
# which has 'cc' in the path
# These are used if you leave FARGO_ARCH undefined
CC_  = cc
OPT_ = -O
OPTSEQ_ =
PARAOPT_ =
PARACC_ = cc

# Setup for Compaq HP/SC (Digital EV68)
# with no OpenMP support
# Set the FARGO_ARCH variable to EV68
CC_EV68  = cc
OPT_EV68 = -arch ev6 -fast -O3 -inline all -msg_disable ompdirignored
OPTSEQ_EV68 =
PARAOPT_EV68 =
PARACC_EV68 = cc

# Setup for Compaq HP/SC (Digital EV68)
# with OpenMP support
# Set the FARGO_ARCH variable to EVMP
CC_EVMP  = cc
OPT_EVMP = -arch ev6 -fast -O3 -inline all -mp
OPTSEQ_EVMP =
PARAOPT_EVMP =
PARACC_EVMP = cc

# Generic Sparc (v7) SunStation platform w/ gcc
# FARGO_ARCH must be set to SPARC
CC_SPARC  = gcc
OPT_SPARC = -O3 -Wall
OPTSEQ_SPARC =
PARAOPT_SPARC =
PARACC_SPARC = mpicc
#
#
#
#--------------------No Changes Needed after this line------------------------
#
#
#
SHELL		=  /bin/sh

MAINOBJ         = LowTasks.o SideEuler.o Output.o Init.o main.o Theo.o\
		  Interpret.o SourceEuler.o TransportEuler.o Stockholm.o\
		  Planet.o RungeKunta.o Viscosity.o Psys.o Force.o\
		  var.o Pframeforce.o split.o merge.o commbound.o fpe.o rebin.o

DUMMY		= mpi_dummy.o

COMP        = $(CC_$(FARGO_ARCH))
OPT         = $(OPT_$(FARGO_ARCH))
OPTSEQ      = $(OPTSEQ_$(FARGO_ARCH))
PARAOPT     = $(PARAOPT_$(FARGO_ARCH)) -D_PARALLEL
PARACOMP    = $(PARACC_$(FARGO_ARCH))
#Uncomment the following line for OpenMP support (with gcc) and comment the
#subsequent line
#LIBS        = -lm -lgomp
LIBS        = -lm
PARALIBS    = -lmpi
AUTOINCL    = param.h param_noex.h global_ex.h

include	.config

EXENAME        = ../fargo

ARCHIVE		= $(EXENAME:../%=%.tar)
ARCHIVECOMP	= $(EXENAME:../%=%.tar.gz)

SRC = *.c
INCLUDE = *.h

ifeq ($(BUILD),parallel)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS)
OPTIONS		= $(OPT) $(PARAOPT)
OBJ		= $(MAINOBJ)
else
COMPILER	= $(COMP)
LIBRARIES	= $(LIBS)
OPTIONS		= $(OPT) $(OPTSEQ)
OBJ		= $(MAINOBJ) $(DUMMY)
endif

all: conditionalrebuild $(AUTOINCL) $(OBJ) $(EXENAME) archive
	@echo ""
	@echo ""
	@echo "      NOTE"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "This built is PARALLEL (MPI)."
else
	@echo "This built is SEQUENTIAL."
endif
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "gmake BUILD=sequential"
else
	@echo "gmake BUILD=parallel"
endif
	@echo ""

$(EXENAME): $(OBJ)
	$(COMPILER) $(OBJ) $(OPTIONS) -o $(EXENAME) $(LIBRARIES)

.PHONY: conditionalrebuild
ifneq ($(BUILD),$(OLDBUILD))
conditionalrebuild: clean
	@echo "OLDBUILD = $(BUILD)" > .config
	@echo "BUILD = $(BUILD)" >> .config
else
conditionalrebuild:
endif

.oldconfig:
.config:

archive : $(SRC) $(INCL) makefile varparser.pl
	@echo "Creating ../.source.tar.bz2"
	@tar cf ../.source.tar *.c
	@tar rf ../.source.tar *.h
	@tar rf ../.source.tar makefile
	@tar rf ../.source.tar varparser.pl
	@bzip2 -9 -f ../.source.tar

para:
	@gmake BUILD=parallel

seq:
	@gmake BUILD=sequential

$(AUTOINCL) : var.c global.h makefile varparser.pl
	@./varparser.pl

$(OBJ): fargo.h fondam.h param.h param_noex.h types.h makefile

.PHONY: clean mrproper package

mrproper:
	rm -f *.o *~ *.s *.il $(AUTOINCL) $(EXENAME) ../core.*\
	*.tex *.dvi *.pdf *.ps *.log *.aux *.lint $(ARCHIVE)\
	$(ARCHIVECOMP)

clean:
	rm -f *.o *~ *.s *.il

package: $(ARCHIVECOMP)

# The Doxyfile provided here is relatively general and should allow
# you to build your own documentation if you modify the source. You
# should have the required executable files in your path
# (doxygen, dot, latex, epstopdf, dvips, gs). On Mac Os X I made soft
# links in /usr/bin to these executables (either in /sw/bin or
# /Applications/Doxygen.app/Content/Resources).
doc:
	doxygen Doxyfile

$(ARCHIVECOMP) : $(SRC) $(INCL) makefile varparser.pl	Doxyfile
	@echo "Creating $(ARCHIVECOMP)"
	@cd ../..;tar cf $(ARCHIVE) fargo/src/*.c
	@cd ../..;tar rf $(ARCHIVE) fargo/src/*.h
	@cd ../..;tar rf $(ARCHIVE) fargo/src/varparser.pl
	@cd ../..;tar rf $(ARCHIVE) fargo/src/makefile
	@cd ../..;tar rf $(ARCHIVE) fargo/src/Doxyfile
	@cd ../..;tar rf $(ARCHIVE) fargo/in/template.par
	@cd ../..;tar rf $(ARCHIVE) fargo/in/Jup.cfg
	@cd ../..;tar rf $(ARCHIVE) fargo/in/Nep.cfg
	@cd ../..;tar rf $(ARCHIVE) fargo/in/stockholm*par
	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out1/
	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out2/
	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out3/
	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out4/
	@cd ../..;tar rf $(ARCHIVE) fargo/mp.pro
	@cd ../..;tar rf $(ARCHIVE) fargo/INSTALL
	@cd ../..;gzip -9 $(ARCHIVE)
	@cd ../..;mv $(ARCHIVECOMP) fargo/

.c.o  :
	$(COMPILER) $*.c -c $(OPTIONS)
