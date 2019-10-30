# Track git commit when in git repository
GIT_COMMIT = "$(shell git rev-parse HEAD)"
GIT_CHANGED = "$(shell git diff-index --name-only HEAD | sed -e '{:q;N;s/\n/, /g;t q}')"

ifeq (,$(subst "",,$(GIT_COMMIT)$(GIT_CHANGED)))
	GIT_COMMIT = "not in git repo"
	GIT_CHANGED = "not in git repo"
endif

# !! Define Environment variable FARGO_ARCH to reflect your architecture
# for ex.: setenv FARGO_ARCH LINUX (or export FARGO_ARCH=LINUX depending on your shell)
# Possible choices : undefined, LINUX, BINAC

# Generic Linux platform
# FARGO_ARCH must be set to LINUX
CC_LINUX = mpicc
CXX_LINUX = mpic++
CFLAGS_LINUX = -Ofast -std=c++11 -Wall -march=native -Wno-unknown-pragmas
#CFLAGS_LINUX = -g -Wall -march=native -Wno-unknown-pragmas
MPI_PREFIX_LINUX = /usr
FFTW_PREFIX_LINUX = /usr
GSL_PREFIX_LINUX = /usr

# Generic Linux platform
# FARGO_ARCH must be set to LINUXENV
CC_LINUXENV = mpicc
CXX_LINUXENV = mpicxx
CFLAGS_LINUXENV = -Ofast -std=c++11 -Wall -march=native -Wno-unknown-pragmas
MPI_PREFIX_LINUXENV = $(MPI_HOME)
FFTW_PREFIX_LINUXENV = $(FFTW_HOME)
GSL_PREFIX_LINUXENV = /usr/gsl/2.1

# HLR2
# FARGO_ARCH must be set to HLR2
CC_HLR2 = mpicc
CXX_HLR2 = mpic++
CFLAGS_HLR2 = -Ofast -Wall -march=native -Wno-unknown-pragmas
MPI_PREFIX_HLR2 = /opt/openmpi/1.10
FFTW_PREFIX_HLR2 = /opt/fftw/3.3_openmpi_1.10
GSL_PREFIX_HLR2 = /opt/gsl/2.1

# BINAC
# FARGO_ARCH must be set to BINAC
CC_BINAC = mpicc
CXX_BINAC = mpicxx
CFLAGS_BINAC = -Ofast -Wall -march=native -Wno-unknown-pragmas -std=c++11
#MPI_PREFIX_BINAC = /opt/bwhpc/common/mpi/openmpi/1.10.3-gnu-6.1
MPI_PREFIX_BINAC = $(MPI_HOME)
#FFTW_PREFIX_BINAC = /opt/bwhpc/common/numlib/fftw/3.3.5-gnu-6.1
FFTW_PREFIX_BINAC = $(FFTW_HOME)
GSL_PREFIX_BINAC = /opt/bwhpc/common/numlib/gsl/2.1

# Default Arch is LINUX
CC_ = $(CC_LINUX)
CXX_ = $(CXX_LINUX)
CFLAGS_ =  $(CFLAGS_LINUX)
MPI_PREFIX_ = $(MPI_PREFIX_LINUX)
FFTW_PREFIX_ = $(FFTW_PREFIX_LINUX)
GSL_PREFIX_ = $(GSL_PREFIX_LINUX)

# # current git commit
# GIT_COMMIT = "$(shell git rev-parse HEAD)"
# GIT_CHANGED = "$(shell git diff-index --name-only HEAD | sed -e '{:q;N;s/\n/, /g;t q}')"

# options to compile
OPTIONS := -D_TRAP_FPE -D_GNU_SOURCE -DGIT_COMMIT='$(GIT_COMMIT)' -DGIT_CHANGED='$(GIT_CHANGED)'
ifneq "$(FARGO_DEBUG)" "1"
	OPTIONS := $(OPTIONS) -DNDEBUG
endif

# (MPI)-C compiler to use
CC = $(CC_$(FARGO_ARCH))
# (MPI)-C++ compiler to use
CXX = $(CXX_$(FARGO_ARCH))

# C Flags
CFLAGS := $(CFLAGS_$(FARGO_ARCH))
ifeq "$(FARGO_DEBUG)" "1"
	CFLAGS := $(CFLAGS) -g
endif
$(info CFLAGS=$(CFLAGS))
# C++ Flags
CXXFLAGS = $(CFLAGS)

# use special MPI/FFTW libs
MPI_PREFIX = $(MPI_PREFIX_$(FARGO_ARCH))
FFTW_PREFIX= $(FFTW_PREFIX_$(FARGO_ARCH))
GSL_PREFIX= $(GSL_$(FARGO_ARCH))

_CFLAGS = $(OPTIONS) $(CFLAGS) -I$(FFTW_PREFIX)/include -I$(MPI_PREFIX)/include -I$(GSL_PREFIX)/include
_CXXFLAGS = $(OPTIONS) $(CXXFLAGS) -I$(FFTW_PREFIX)/include -I$(MPI_PREFIX)/include -I$(GSL_PREFIX)/include

SRCS_C = $(wildcard *.c)
SRCS_CXX = $(wildcard *.cpp)
SRCS = $(SRCS_C) $(SRCS_CXX)
OBJS_C = $(SRCS_C:%.c=%.o)
OBJS_CXX = $(SRCS_CXX:%.cpp=%.o)
OBJS = $(OBJS_C) $(OBJS_CXX)
HDRS = $(wildcard *.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

MAKEFILE_DIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

LIBS = -lm -L$(MPI_PREFIX)/lib -lmpi -L$(FFTW_PREFIX)/lib -I$(FFTW_PREFIX)/include -lfftw3_mpi -lfftw3 -L$(GSL_PREFIX)/lib64 -I$(GSL_PREFIX)/include -lgsl -lgslcblas -L$(MAKEFILE_DIR)rebound -lrebound
# -Wl,--rpath -Wl,$(FFTW_PREFIX)/lib

EXENAME        = ../fargo

$(EXENAME): $(OBJS)
	$(MAKE) -C rebound
	$(CXX) $(OBJS) $(_CFLAGS) -o $(EXENAME) $(LIBS) 

$(OBJS): $(HDRS) makefile

main.o: $(HDRS) makefile $(SRCS_C) $(SRCS_CXX)

.PHONY: clean

clean:
	rm -f $(OBJS) $(TEMP)
	$(MAKE) -C rebound clean

%.o : %.c
	$(CC) $*.c -c $(_CFLAGS)

%.o : %.cpp
	$(CXX) $*.cpp -c $(_CXXFLAGS)
