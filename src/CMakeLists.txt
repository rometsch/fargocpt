#########
# Tested with following modules on the BINAC cluster
#
# module load numlib/gsl/2.5-gnu-9.2
# module load mpi/openmpi/3.1-gnu-9.2
# module load numlib/fftw/3.3.8-openmpi-3-gnu-9.2
# module load devel/cmake/3.27.5
#
#########
cmake_minimum_required(VERSION 3.1)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

PROJECT(fargocpt)

# where the compiled program is located
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# find MPI, needed to check whether we are on Binac
FIND_PATH(MPI_INCLUDE_DIR
        NAME
        mpi.h
        PATHS
        ${MPI_INCLUDE_DIRS}
        /opt/bwhpc/common/mpi/openmpi/3.1.6-gnu-9.2/include/
        DOC "The directory where mpi.h resides")
        set(MPI_INCLUDE_DIRS ${MPI_INCLUDE_DIR})
# MPI dir to check if we are on BINAC
string(SUBSTRING ${MPI_INCLUDE_DIR} 0 11 MPI_SUBSTR)

find_package(MPI REQUIRED)

if("${MPI_SUBSTR}" STREQUAL "/opt/bwhpc/")

        set(BINAC_FFTW_PATH "/opt/bwhpc/common/numlib/fftw/3.3.8-openmpi-3.1-gnu-9.2")
        set(BINAC_GASL_PATH "/opt/bwhpc/common/numlib/gsl/2.5-gnu-9.2")

        include_directories("${BINAC_FFTW_PATH}/include/" "${BINAC_GSL_PATH}/include/")
	
	find_library(FFTW3_LIBRARY fftw3 HINTS "${BINAC_FFTW_PATH}/lib/")

	find_library(FFTW3_MPI_LIBRARY fftw3_mpi HINTS "${BINAC_FFTW_PATH}/lib/")
	find_library(FFTW3_OMP_LIBRARY fftw3_omp HINTS "${BINAC_FFTW_PATH}/lib/")
	find_library(GSL_LIBRARY gsl HINTS "${BINAC_GSL_PATH}/lib/")
	find_library(GSLCBLAS_LIBRARY gslcblas HINTS "${BINAC_GSL_PATH}/lib/")

    	find_package(GSL REQUIRED)
else()  
	
	find_library(FFTW3_LIBRARY fftw3 HINTS "")
	find_library(FFTW3_MPI_LIBRARY fftw3_mpi HINTS "")
	find_library(FFTW3_OMP_LIBRARY fftw3_omp HINTS "")
	find_library(GSL_LIBRARY gsl HINTS "")
	find_library(GSLCBLAS_LIBRARY gslcblas HINTS "")

endif()


########## Include OpenMP ############################################################################################################
find_package(OpenMP REQUIRED)

# Copy & Pasted from https://cliutils.gitlab.io/modern-cmake/chapters/packages/OpenMP.html
# needed to link OpenMP lib
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()
########## END Include OpenMP ########################################################################################################

############ Build Units ################################################################################################################
FILE(GLOB UNITS_CPP_SOURCES "units/*.cpp"            )
FILE(GLOB UNITS_H_SOURCES   "units/*.h" "units/*.hpp")
ADD_LIBRARY(LocalUnitsLib STATIC ${UNITS_CPP_SOURCES} ${UNITS_H_SOURCES})
target_compile_options(LocalUnitsLib PRIVATE -fPIC -O3 -Wall -Wextra -Wshadow -Wno-long-long -pedantic -pedantic-errors -std=gnu++17)
############ END Buld Units #############################################################################################################

############ Build Rebound ################################################################################################################
FILE(GLOB REB_CPP_SOURCES "rebound/*.c")
FILE(GLOB REB_H_SOURCES   "rebound/*.h")
ADD_LIBRARY(LocalReboundLib STATIC ${REB_CPP_SOURCES} ${REB_H_SOURCES})
target_compile_options(LocalReboundLib PRIVATE -std=c99 -fPIC -Wpointer-arith -O3 -Wall -Wno-unused-result -lfftw3 -fopenmp)
target_compile_definitions(LocalReboundLib PRIVATE _GNU_SOURCE OPENMP FFTW)
target_link_libraries(LocalReboundLib -lm -lrt)
############ End Build Rebound ############################################################################################################


############ Build Yaml-cpp ###############################################################################################################
FILE(GLOB YAML_CPP_SOURCES "yaml-cpp/*.cpp" "yaml-cpp/contrib/*.cpp" "yaml-cpp/node/*.cpp" "yaml-cpp/node/detail/*.cpp")
FILE(GLOB YAML_H_SOURCES   "yaml-cpp/*.h" "yaml-cpp/contrib/*.h" "yaml-cpp/node/*.h" "yaml-cpp/node/detail/*.h")

ADD_LIBRARY(LocalYamlLib STATIC ${YAML_CPP_SOURCES} ${YAML_H_SOURCES})
target_compile_options(LocalYamlLib PRIVATE -fPIC -O3 -Wall -Wextra -Wshadow -Weffc++ -Wno-long-long -pedantic -pedantic-errors -std=gnu++11)
target_compile_definitions(LocalYamlLib PRIVATE YAML_CPP_STATIC_DEFINE)
target_include_directories(LocalYamlLib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
FILE(GLOB CPP_SOURCES "./*.cpp" "boundary_conditions/*.cpp" "nbody/*.cpp" "particles/*.cpp" "random/*.cpp"            "viscosity/*.cpp")
FILE(GLOB H_SOURCES   "./*.h"   "boundary_conditions/*.h"   "nbody/*.h"   "particles/*.h"   "random/.h" "random/.hpp" "viscosity/*.h"  )
############ END Build Yaml-cpp ###########################################################################################################


############ Build fargocpt ###############################################################################################################
# Use imported targets would make things much eazier. Thanks Levi for pointing it out.
add_executable(${PROJECT_NAME} ./main.cpp ${CPP_SOURCES} ${H_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_CXX_INCLUDE_PATH})




target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wextra -Wmissing-declarations -Werror=redundant-decls -Werror=format-security -Wno-unknown-pragmas)
target_compile_options(${PROJECT_NAME} PUBLIC -std=c++17 -Ofast -march=native -flto=auto -fopenmp)
# add symbols for backtrace
if (UNIX)
    target_compile_options(${PROJECT_NAME} PUBLIC -rdynamic)
endif (UNIX)
target_compile_definitions(${PROJECT_NAME} PUBLIC OMPI_SKIP_MPICXX NDEBUG _TRAP_FPE _GNU_SOURCE OPENMP)

target_link_libraries(${PROJECT_NAME} -lm ${MPI_CXX_LIBRARIES} ${FFTW3_LIBRARY} ${FFTW3_MPI_LIBRARY} ${FFTW3_OMP_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} OpenMP::OpenMP_CXX LocalReboundLib LocalUnitsLib LocalYamlLib)
target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
