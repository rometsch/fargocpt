include makefile.defs

BUILD_DIR = ../build

SRCS = $(wildcard *.cpp particles/*.cpp random/*.cpp nbody/*.cpp viscosity/*.cpp)
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
HDRS = $(wildcard *.h rebound/*.h particles/*.h random/*.h nbody/*.h viscosity/*.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

MAKEFILE_DIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

$(info Build parameters:)
$(info CFLAGS: $(CFLAGS))
$(info INCLUDES: $(INCLUDES))
$(info LIBS: $(LIBS))
$(info GIT_INFO: $(GIT_INFO))

build: $(EXENAME)

$(EXENAME): rebound/librebound.a yaml-cpp/libyamlcpp.a units/libunits.a $(OBJS)
	@mkdir -p $(@D)
	@$(CXX) $(OBJS) $(CFLAGS) $(INCLUDES) $(LIBS) $(GIT_INFO) -o $(EXENAME) $(LIBS)
	@echo "Compiling and linking $@"
	@echo "Build successful!"

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
	@$(CXX) -c $(CFLAGS) $(INCLUDES) $(LIBS) $(OPTIONS) -o $@ $<

list:
	@echo $(OBJS)

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean list

clean:
	rm -f $(OBJS) $(TEMP)
	cd rebound && $(MAKE) clean
	cd yaml-cpp && $(MAKE) clean
	cd units && $(MAKE) clean
	rm -rf $(BUILD_DIR)
