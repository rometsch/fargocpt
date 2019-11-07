include makefile.defs

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
