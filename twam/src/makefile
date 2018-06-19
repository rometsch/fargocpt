# !! Define Environment variable FARGO_ARCH to reflect your architecture
# for ex.: setenv FARGO_ARCH LINUX (or export FARGO_ARCH=LINUX depending on your shell)
# Possible choices : undefined, LINUX, BINAC

# Generic Linux platform
# FARGO_ARCH must be set to LINUX
CC_LINUX = mpicc
CXX_LINUX = mpic++
CFLAGS_LINUX = -O2 -Wall -march=native -Wno-unknown-pragmas -ggdb
MPI_PREFIX_LINUX = /usr
FFTW_PREFIX_LINUX = /usr
GSL_PREFIX_LINUX = /usr

# BINAC
# FARGO_ARCH must be set to BINAC
CC_BINAC = /opt/bwhpc/common/mpi/openmpi/2.1.2-gnu-5.2/bin/mpicc
CXX_BINAC = /opt/bwhpc/common/mpi/openmpi/2.1.2-gnu-5.2/bin/mpic++
CFLAGS_BINAC = -O2 -Wall -march=native -Wno-unknown-pragmas -ggdb
MPI_PREFIX_BINAC = /opt/bwhpc/common/mpi/openmpi/2.1.2-gnu-5.2
FFTW_PREFIX_BINAC = /opt/bwhpc/common/numlib/fftw/3.3.5-gnu-5.2
GSL_PREFIX_LINUX = /opt/bwhpc/common/numlib/gsl/2.1

# Default Arch is LINUX
CC_ = $(CC_LINUX)
CXX_ = $(CXX_LINUX)
CFLAGS_ =  $(CFLAGS_LINUX)
MPI_PREFIX_ = $(MPI_PREFIX_LINUX)
FFTW_PREFIX_ = $(FTTW_PREFIX_LINUX)
GSL_PREFIX_ = $(GSL_PREFIX_LINUX)

# # current git commit
# GIT_COMMIT = "$(shell git rev-parse HEAD)"
# GIT_CHANGED = "$(shell git diff-index --name-only HEAD | sed -e '{:q;N;s/\n/, /g;t q}')"

# options to compile
OPTIONS = -D_TRAP_FPE -D_GNU_SOURCE -DGIT_COMMIT='$(GIT_COMMIT)' -DGIT_CHANGED='$(GIT_CHANGED)'

# (MPI)-C compiler to use
CC = $(CC_$(FARGO_ARCH))
# (MPI)-C++ compiler to use
CXX = $(CXX_$(FARGO_ARCH))

# C Flags
CFLAGS = $(CFLAGS_$(FARGO_ARCH))
# C++ Flags
CXXFLAGS = $(CFLAGS)

# use special MPI/FFTW libs
MPI_PREFIX = $(MPI_PREFIX_$(FARGO_ARCH))
FFTW_PREFIX= $(FTTW_PREFIX_$(FARGO_ARCH))
GSL_PREFIX= $(FTTW_GSL_$(FARGO_ARCH))

_CFLAGS = $(OPTIONS) $(CFLAGS) -I$(FFTW_INC_DIR) -I$(MPI_INC_DIR) -I$(GSL_PREFIX)/include
_CXXFLAGS = $(OPTIONS) $(CXXFLAGS) -I$(FFTW_INC_DIR) -I$(MPI_INC_DIR) -I$(GSL_PREFIX)/include

SRCS_C = $(wildcard *.c)
SRCS_CXX = $(wildcard *.cpp)
SRCS = $(SRCS_C) $(SRCS_CXX)
OBJS_C = $(SRCS_C:%.c=%.o)
OBJS_CXX = $(SRCS_CXX:%.cpp=%.o)
OBJS = $(OBJS_C) $(OBJS_CXX)
HDRS = $(wildcard *.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

LIBS = -lm -L$(MPI_LIB_DIR) -lmpi -L$(FFTW_LIB_DIR) -I$(FFTW_INC_DIR) -lfftw3_mpi -lfftw3 -L$(GSL_PREFIX)/lib64 -I$(GSL_PREFIX)/include -lgsl -lgslcblas
# -Wl,--rpath -Wl,$(FFTW_PREFIX)/lib

EXENAME        = ../fargo

$(EXENAME): $(OBJS)
	$(CXX) $(OBJS) $(_CFLAGS) -o $(EXENAME) $(LIBS)

$(OBJS): $(HDRS) makefile

main.o: $(HDRS) makefile $(SRCS_C) $(SRCS_CXX)

.PHONY: clean

clean:
	rm -f $(OBJS) $(TEMP)

%.o : %.c
	$(CC) $*.c -c $(_CFLAGS)

%.o : %.cpp
	$(CXX) $*.cpp -c $(_CXXFLAGS)
