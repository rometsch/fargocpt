#########
# Tested with following modules on the BINAC cluster
#
# module load numlib/gsl/2.1
# module load mpi/openmpi/3.1-gnu-7.3
# module load numlib/fftw/3.3.8-openmpi-3.1-gnu-7.3
#
#########
cmake_minimum_required(VERSION 3.1)

PROJECT(rebound)
# include files
FILE(GLOB CPP_SOURCES "./*.c")
FILE(GLOB H_SOURCES "./*.h")

#add_library(rebound SHARED  ${CPP_SOURCES} ${H_SOURCES})
