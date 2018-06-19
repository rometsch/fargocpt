# current git commit
GIT_COMMIT = "$(shell git rev-parse HEAD)"
GIT_CHANGED = "$(shell git diff-index --name-only HEAD | sed -e '{:q;N;s/\n/, /g;t q}')"

# options to compile
OPTIONS = -D_TRAP_FPE -D_GNU_SOURCE -DGIT_COMMIT='$(GIT_COMMIT)' -DGIT_CHANGED='$(GIT_CHANGED)'

# (MPI)-C compiler to use
CC = mpicc
# (MPI)-C++ compiler to use
CXX = mpic++

# C Flags
CFLAGS = -O3 -Wall -march=native -Wno-unknown-pragmas -ggdb
# C++ Flags
CXXFLAGS = $(CFLAGS)

# use special MPI/FFTW libs
MPI_PREFIX = /usr
FFTW_PREFIX= /usr/local

_CFLAGS = $(OPTIONS) $(CFLAGS) -I$(FFTW_PREFIX)/include -I$(MPI_PREFIX)/include
_CXXFLAGS = $(OPTIONS) $(CXXFLAGS) -I$(FFTW_PREFIX)/include -I$(MPI_PREFIX)/include

SRCS_C = $(wildcard *.c)
SRCS_CXX = $(wildcard *.cpp)
SRCS = $(SRCS_C) $(SRCS_CXX)
OBJS_C = $(SRCS_C:%.c=%.o)
OBJS_CXX = $(SRCS_CXX:%.cpp=%.o)
OBJS = $(OBJS_C) $(OBJS_CXX)
HDRS = $(wildcard *.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

LIBS = -lm -L$(MPI_PREFIX)/lib -lmpi -L$(FFTW_PREFIX)/lib -I$(FFTW_PREFIX)/include -lfftw3_mpi -lfftw3 -lgsl -lgslcblas
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
