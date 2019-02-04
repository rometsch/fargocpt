#########
# Tested with following modules on the BINAC cluster
#
# module load numlib/gsl/2.1
# module load mpi/openmpi/3.1-gnu-7.3
# module load numlib/fftw/3.3.8-openmpi-3.1-gnu-7.3
#
#########
cmake_minimum_required(VERSION 3.1)

PROJECT(fargo)

# find MPI 
FIND_PATH( MPI_INCLUDE_DIR
        NAME
        mpi.h
        PATHS
        ${MPI_INCLUDE_DIRS}
        /opt/bwhpc/common/mpi/openmpi/3.1.2-gnu-7.3/include/	/usr/local/openmpi-4.0/include/
        DOC "The directory where mpi.h resides")
        set(MPI_INCLUDE_DIRS ${MPI_INCLUDE_DIR})

FIND_LIBRARY( MPI_LIBRARY
        NAMES MPI mpi
        PATHS /opt/bwhpc/common/mpi/openmpi/3.1.2-gnu-7.3/	/usr/local/openmpi-4.0/lib/
        DOC "The MPI library")
        SET(MPI_LIBRARIES ${MPI_LIBRARY})
# where the compiled program is located

# find fftw3-mpi.h
FIND_PATH( FFTW_MPI_INCLUDE_DIR
        NAME
	fftw3-mpi.h
        PATHS
	${FFTW_MPI_INCLUDE_DIRS}
	/usr/local/fftw-3.3.8/include/
        DOC "The directory where fftw3-mpi.h resides")
	set(FFTW_MPI_INCLUDE_DIRS ${FFTW_MPI_INCLUDE_DIR})




set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# include files
FILE(GLOB CPP_SOURCES "./*.cpp")
FILE(GLOB H_SOURCES "./*.h")

# MPI dir to check if we are on BINAC
string(SUBSTRING ${MPI_INCLUDE_DIR} 0 11 MPI_SUBSTR)

if("${MPI_SUBSTR}" STREQUAL "/opt/bwhpc/")
        set(BINAC_FFTW_PATH "/opt/bwhpc/common/numlib/fftw/3.3.8-openmpi-3.1-gnu-7.3")
        set(BINAC_GASL_PATH  "/opt/bwhpc/common/numlib/gsl/2.1")
        include_directories("${BINAC_FFTW_PATH}/include/" "${BINAC_GSL_PATH}/include/")
        #link_directories("${BINAC_GSL_PATH}/lib/" "${BINAC_FFTW_PATH}/lib/")
	
	find_library(FFTW3_LIBRARY fftw3 HINTS "${BINAC_FFTW_PATH}/lib/")
	find_library(FFTW3_MPI_LIBRARY fftw3_mpi HINTS "${BINAC_FFTW_PATH}/lib/")
	find_library(GSL_LIBRARY gsl HINTS "${BINAC_GSL_PATH}/lib/")
	find_library(GSLCBLAS_LIBRARY gslcblas HINTS "${BINAC_GSL_PATH}/lib/")
else()  
        # otherwise cmake finds gsl 1.5 on binac
        find_package(GSL REQUIRED)
	
	find_library(FFTW3_LIBRARY fftw3 HINTS "/usr/local/fftw-3.3.8/lib/")
	find_library(FFTW3_MPI_LIBRARY fftw3_mpi HINTS "/usr/local/fftw-3.3.8/lib/")
	find_library(GSL_LIBRARY gsl HINTS "")
	find_library(GSLCBLAS_LIBRARY gslcblas HINTS "")

endif()

find_package(MPI REQUIRED)


add_compile_options(-std=c++17 -O3 -Wall -Wextra -Wno-old-style-cast -march=native -DNDEBUG -D_TRAP_FPE -D_GNU_SOURCE) 
# add_compile_options(-std=c++17 -g -Wall -Wextra -Wno-old-style-cast -march=native -D_TRAP_FPE -D_GNU_SOURCE) 

# Use imported targets would make things much eazier. Thanks Levi for pointing it out.
add_executable(${PROJECT_NAME} ./main.cpp ${CPP_SOURCES} ${H_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_CXX_INCLUDE_PATH} ${FFTW_MPI_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} ${FFTW3_MPI_LIBRARY} ${FFTW3_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
# target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} fftw3 fftw3_mpi gsl gslcblas)
