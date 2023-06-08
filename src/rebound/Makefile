include Makefile.defs

SOURCES=rebound.c tree.c particle.c gravity.c integrator.c integrator_whfast.c integrator_saba.c integrator_ias15.c integrator_sei.c integrator_bs.c integrator_leapfrog.c integrator_mercurius.c integrator_eos.c integrator_tes.c boundary.c input.c binarydiff.c output.c collision.c communication_mpi.c display.c tools.c rotations.c derivatives.c simulationarchive.c glad.c integrator_janus.c transformations.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) librebound.a

%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) -o $@ $<

librebound.a: $(OBJECTS)
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
	@echo "Removing coverage files ..."
	@-rm -f *.gcda
	@-rm -f *.gcno
	
