include makefile.defs

SRCS_C = $(wildcard *.c)
SRCS_CXX = $(wildcard *.cpp)
SRCS = $(SRCS_C) $(SRCS_CXX)
OBJS_C = $(SRCS_C:%.c=%.o)
OBJS_CXX = $(SRCS_CXX:%.cpp=%.o)
OBJS = $(OBJS_C) $(OBJS_CXX)
HDRS = $(wildcard *.h) $(wildcard rebound/*.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

MAKEFILE_DIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

LIBS = -lm -lstdc++fs -L$(MPI_PREFIX)/lib -lmpi -L$(FFTW_PREFIX)/lib -I$(FFTW_PREFIX)/include -lfftw3_mpi -lfftw3 -L$(GSL_PREFIX)/lib64 -I$(GSL_PREFIX)/include -lgsl -lgslcblas

#statically link rebound
LIBS += -Lrebound -lrebound

EXENAME        = ../fargo

$(EXENAME): rebound/librebound.a $(OBJS)
	$(CXX) $(OBJS) $(_CFLAGS) -o $(EXENAME) $(LIBS)

rebound/librebound.a:
	$(MAKE) -C rebound

$(OBJS): $(HDRS) makefile

main.o: $(HDRS) makefile $(SRCS_C) $(SRCS_CXX)

.PHONY: clean

clean:
	rm -f $(OBJS) $(TEMP)
	cd rebound && $(MAKE) clean

%.o : %.c
	$(CC) $(_CFLAGS) -c -o $@ $<

%.o : %.cpp
	$(CXX) $*.cpp -c $(_CXXFLAGS)
