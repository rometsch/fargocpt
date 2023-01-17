include makefile.defs

BUILD_DIR = ../build

SRCS = $(wildcard *.cpp particles/*.cpp random/*.cpp nbody/*.cpp viscosity/*.cpp)
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
HDRS = $(wildcard *.h rebound/*.h particles/*.h random/*.h nbody/*.h viscosity/*.h)
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

MAKEFILE_DIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

$(info $(OBJS))

build: $(EXENAME)

# $(EXENAME): rebound/librebound.a yaml-cpp/libyamlcpp.a units/libunits.a $(OBJS)
$(EXENAME): $(OBJS)
	# mkdir -p $(@D)
	@$(CXX) $(OBJS) $(_CFLAGS) -o $(EXENAME) $(LIBS)
	@echo Building $@ with the following arguments:
	@echo $(_CFLAGS)
	@echo and with libraries:
	@echo $(LIBS)
	@echo "Build successful!"

# rebound/librebound.a:
# 	$(MAKE) -C rebound

# yaml-cpp/libyamlcpp.a:
# 	$(MAKE) -C yaml-cpp

# units/libunits.a:
# 	$(MAKE) -C units

$(OBJS): $(HDRS) makefile

main.o: $(HDRS) makefile $(SRCS)

list:
	@echo $(OBJS)

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean list

clean:
	rm -f $(OBJS) $(TEMP)
	cd rebound && $(MAKE) clean
	cd yaml-cpp && $(MAKE) clean
	cd units && $(MAKE) clean

$(OBJS): $(BUILD_DIR)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CXX) -c $(_CXXFLAGS) -o $@ $<