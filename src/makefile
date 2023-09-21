include makefile.defs

# CFLAGS += -DDO_USE_PLUTO_UNITS

SRCS = $(wildcard *.cpp particles/*.cpp random/*.cpp nbody/*.cpp viscosity/*.cpp boundary_conditions/*.cpp)
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
HDRS = $(wildcard *.h rebound/*.h particles/*.h random/*.h nbody/*.h viscosity/*.h boundary_conditions/*.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

MAKEFILE_DIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

build: $(EXENAME)

$(EXENAME): rebound/librebound.a yaml-cpp/libyamlcpp.a units/libunits.a $(OBJS)
	@mkdir -p $(@D)
	@echo
	@echo "Compiling and linking $@"
	$(CXX) $(OBJS) $(CFLAGS) $(INCLUDES) $(LIBS) $(GIT_INFO) $(OPTIONS) -o $(EXENAME) $(LIBS)
	@echo
	@echo "Build successful!"
	@$(MAKE) --no-print-directory info

rebound/librebound.a:
	$(MAKE) -C rebound

yaml-cpp/libyamlcpp.a:
	$(MAKE) -C yaml-cpp

units/libunits.a:
	$(MAKE) -C units

main.o: $(HDRS) makefile $(SRCS)

$(OBJS): $(HDRS) makefile

# http://www.gnu.org/software/make/manual/make.html#Static-Pattern
$(OBJS): $(BUILD_DIR)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CXX) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) -o $@ $<

info:
	@echo
	@echo Build parameters were:
	@echo
	@echo Compiler: 
	$(CC) --version
	@echo CFLAGS:
	@echo $(CFLAGS)
	@echo
	@echo OPTIONS:
	@echo $(OPTIONS)
	@echo
	@echo INCLUDES:
	@echo $(INCLUDES)
	@echo
	@echo LIBS:
	@echo $(LIBS)
	@echo
	@echo GIT_INFO:
	@echo $(GIT_INFO)
	@echo
	@echo Include path MPI:  INCLUDE_MPI:  $(INCLUDE_MPI)
	@echo Include path FFTW: INCLUDE_FFTW: $(INCLUDE_FFTW)
	@echo Include path GSL:  INCLUDE_GSL:  $(INCLUDE_GSL)
	@echo Include path OMP:  INCLUDE_OMP:  $(INCLUDE_OMP)

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean info

clean:
	rm -f $(OBJS) $(TEMP)
	cd rebound && $(MAKE) clean
	cd yaml-cpp && $(MAKE) clean
	cd units && $(MAKE) clean
	rm -rf $(BUILD_DIR)

cleanfargo:
	rm -f $(OBJS) $(TEMP)
	rm -rf $(BUILD_DIR)
