#########
# Tested with following modules on the BINAC cluster
#
# module load numlib/gsl/2.5-gnu-9.2
# module load mpi/openmpi/3.1-gnu-9.2
# module load numlib/fftw/3.3.8-openmpi-3-gnu-9.2
#
#########
cmake_minimum_required(VERSION 3.9)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

PROJECT(fargocpt)

# where the compiled program is located
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# find MPI, needed to check whether we are on Binac
FIND_PATH( MPI_INCLUDE_DIR
        NAME
        mpi.h
        PATHS
        ${MPI_INCLUDE_DIRS}
        /opt/bwhpc/common/mpi/openmpi/3.1-gnu-9.2/include/
        DOC "The directory where mpi.h resides")
        set(MPI_INCLUDE_DIRS ${MPI_INCLUDE_DIR})
# MPI dir to check if we are on BINAC
string(SUBSTRING ${MPI_INCLUDE_DIR} 0 11 MPI_SUBSTR)

if("${MPI_SUBSTR}" STREQUAL "/opt/bwhpc/")
        set(BINAC_FFTW_PATH "/opt/bwhpc/common/numlib/fftw/3.3.8-openmpi-3.1-gnu-9.2")
        set(BINAC_GASL_PATH "/opt/bwhpc/common/numlib/gsl/2.5-gnu-9.2")
        include_directories("${BINAC_FFTW_PATH}/include/" "${BINAC_GSL_PATH}/include/")
        #link_directories("${BINAC_GSL_PATH}/lib/" "${BINAC_FFTW_PATH}/lib/")
	
	find_library(FFTW3_LIBRARY fftw3 HINTS "${BINAC_FFTW_PATH}/lib/")
	find_library(FFTW3_MPI_LIBRARY fftw3_mpi HINTS "${BINAC_FFTW_PATH}/lib/")
	find_library(FFTW3_OMP_LIBRARY fftw3_omp HINTS "${BINAC_FFTW_PATH}/lib/")
	find_library(GSL_LIBRARY gsl HINTS "${BINAC_GSL_PATH}/lib/")
	find_library(GSLCBLAS_LIBRARY gslcblas HINTS "${BINAC_GSL_PATH}/lib/")

    # otherwise cmake finds gsl 1.5 on binac
    find_package(GSL REQUIRED)
else()  
	
	find_library(FFTW3_LIBRARY fftw3 HINTS "")
	find_library(FFTW3_MPI_LIBRARY fftw3_mpi HINTS "")
	find_library(FFTW3_OMP_LIBRARY fftw3_omp HINTS "")
	find_library(GSL_LIBRARY gsl HINTS "")
	find_library(GSLCBLAS_LIBRARY gslcblas HINTS "")

endif()

find_package(MPI REQUIRED)

# Include OpenMP
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)

FILE(GLOB UNITS_CPP_SOURCES "units/*.cpp"            )
FILE(GLOB UNITS_H_SOURCES   "units/*.h" "units/*.hpp")
ADD_LIBRARY(LocalUnitsLib STATIC ${UNITS_CPP_SOURCES} ${UNITS_H_SOURCES})
target_compile_options(LocalUnitsLib PRIVATE -fPIC -Wall -Wextra -Wshadow -Wno-long-long -pedantic -pedantic-errors -std=gnu++17)

FILE(GLOB REB_CPP_SOURCES "rebound/*.c")
FILE(GLOB REB_H_SOURCES   "rebound/*.h")
ADD_LIBRARY(LocalReboundLib STATIC ${REB_CPP_SOURCES} ${REB_H_SOURCES})
target_compile_options(LocalReboundLib PRIVATE -std=c99 -Wpointer-arith -O3 -Wall -Wno-unused-result -lfftw3 -fopenmp)
target_compile_definitions(LocalReboundLib PRIVATE _GNU_SOURCE OPENMP FFTW)
target_link_libraries(LocalReboundLib m rt)

FILE(GLOB YAML_CPP_SOURCES "yaml-cpp/*.cpp")
FILE(GLOB YAML_H_SOURCES   "yaml-cpp/*.h"  )
ADD_LIBRARY(LocalYamlLib STATIC ${YAML_CPP_SOURCES} ${YAML_H_SOURCES})
target_compile_options(LocalUnitsLib PRIVATE -fPIC -Wall -Wextra -Wshadow -Weffc++ -Wno-long-long -pedantic -pedantic-errors -std=gnu++11)

FILE(GLOB CPP_SOURCES "./*.cpp" "boundary_conditions/*.cpp" "nbody/*.cpp" "particles/*.cpp" "random/*.cpp"            "viscosity/*.cpp")
FILE(GLOB H_SOURCES   "./*.h"   "boundary_conditions/*.h"   "nbody/*.h"   "particles/*.h"   "random/.h" "random/.hpp" "viscosity/*.h"  )

# Use imported targets would make things much eazier. Thanks Levi for pointing it out.
add_executable(${PROJECT_NAME} ./main.cpp ${CPP_SOURCES} ${H_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_CXX_INCLUDE_PATH})

# add symbols for backtrace
if (UNIX)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -rdynamic")
endif (UNIX)


target_compile_options(${PROJECT_NAME} PRIVATE -std=c++17 -Ofast -march=native -flto=auto -fopenmp)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wmissing-declarations -Werror=redundant-decls -Werror=format-security -Wno-unknown-pragmas)
target_compile_definitions(${PROJECT_NAME} PRIVATE OMPI_SKIP_MPICXX NDEBUG _TRAP_FPE _GNU_SOURCE OPENMP)

#target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX ${MPI_CXX_LIBRARIES} ${FFTW3_LIBRARY} ${FFTW3_MPI_LIBRARY} ${FFTW3_OMP_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} LocalReboundLib LocalUnitsLib LocalYamlLib)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX ${MPI_CXX_LIBRARIES} ${FFTW3_LIBRARY} ${FFTW3_MPI_LIBRARY} ${FFTW3_OMP_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} LocalReboundLib LocalUnitsLib LocalYamlLib)
#target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} fftw3 fftw3_mpi gsl gslcblas)
