# C compiler
CC = gcc
# C++ compiler
CXX = g++
# C flags
CFLAGS = -O2 -Wall -march=native -Wno-unknown-pragmas -Wall
# C++ flags
CXXFLAGS = $(CFLAGS)

# all target
all:	polargrid_generator modeanalysis read1d jeans fargo2gnuplot merge1d extract_particle extract_profile

polargrid_generator: polargrid_generator.o
	$(CXX) $(CXXFLAGS) -o polargrid_generator polargrid_generator.o

modeanalysis: modeanalysis.o
	$(CXX) $(CXXFLAGS) -o modeanalysis modeanalysis.o

read1d: read1d.o
	$(CXX) $(CXXFLAGS) -o read1d read1d.o

jeans: jeans.o
	$(CXX) $(CXXFLAGS) -o jeans jeans.o

merge1d: merge1d.o
	$(CXX) $(CXXFLAGS) -o merge1d merge1d.o

extract_profile: extract_profile.o
	$(CXX) $(CXXFLAGS) -o extract_profile extract_profile.o

extract_particle: extract_particle.o
	$(CXX) $(CXXFLAGS) -o extract_particle extract_particle.o

fargo2gnuplot: fargo2gnuplot.o
	$(CXX) $(CXXFLAGS) -o fargo2gnuplot fargo2gnuplot.o

# clean target
.PHONY: clean

clean:
	rm -f *.o *~ *.il modeanalysis polargrid_generator read1d

%.o : %.c
	$(CC) $*.c -c $(CFLAGS)

%.o : %.cpp
	$(CXX) $*.cpp -c $(CXXFLAGS)

