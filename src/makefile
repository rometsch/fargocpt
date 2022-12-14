include makefile.defs

SRCS_C = $(wildcard *.c)
SRCS_CXX = $(wildcard *.cpp particles/*.cpp random/*.cpp nbody/*.cpp viscosity/*.cpp)
SRCS = $(SRCS_C) $(SRCS_CXX)
OBJS_C = $(SRCS_C:%.c=%.o)
OBJS_CXX = $(patsubst %.cpp,%.o,$(SRCS_CXX))
OBJS = $(OBJS_C) $(OBJS_CXX)
HDRS = $(wildcard *.h rebound/*.h particles/*.h random/*.h nbody/*.h viscosity/*.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

MAKEFILE_DIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

$(EXENAME): rebound/librebound.a yaml-cpp/libyamlcpp.a units/libunits.a $(OBJS)
	$(CXX) $(OBJS) $(_CFLAGS) -o $(EXENAME) $(LIBS)

rebound/librebound.a:
	$(MAKE) -C rebound

yaml-cpp/libyamlcpp.a:
	$(MAKE) -C yaml-cpp

units/libunits.a:
	$(MAKE) -C units

$(OBJS): $(HDRS) makefile

main.o: $(HDRS) makefile $(SRCS_C) $(SRCS_CXX)

.PHONY: clean

clean:
	rm -f $(OBJS) $(TEMP)
	cd rebound && $(MAKE) clean
	cd yaml-cpp && $(MAKE) clean
	cd units && $(MAKE) clean

%.o : %.c
	$(CC) $(_CFLAGS) -c -o $@ $<

%.o : %.cpp
	$(CXX) $*.cpp -c $(_CXXFLAGS) -o $(patsubst %.cpp,%.o,$*.cpp)
