# inspired by the rebound makefile https://github.com/hannorein/rebound/blob/main/src/Makefile

SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
HEADERS=$(wildcard *.h)

CXX = g++ -I..

CXX_DEFINES = -DYAML_CPP_STATIC_DEFINE
CXX_FLAGS = -fPIC -Wall -Wextra -Wshadow  -Wno-long-long -pedantic -pedantic-errors -std=gnu++17

default_target: all

all: $(SOURCES) libunits.a 

%.o : %.cpp $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CXX) $*.cpp -c $(CXX_FLAGS) -o $*.o


libunits.a: $(OBJECTS)
	@echo ""        
	@echo "Linking static library $@ ..."
	
	ar rcs $@ $(OBJECTS)
	
	@echo ""        
	@echo "The static library $@ has been created successfully."
	
	
clean:
	@echo "Removing object files *.o ..."
	@-rm -f *.o
	@echo "Removing static libraries *.a ..."
	@-rm -f *.a
