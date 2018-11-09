cmake_minimum_required(VERSION 3.1)

PROJECT(fargo)
find_package(MPI REQUIRED)
# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)

# set(MPI_EXECUTABLE_SUFFIX ".openmpi")

FILE(GLOB CPP_SOURCES "./*.cpp")
FILE(GLOB H_SOURCES "./*.h")

include_directories(${MPI_INCLUDE_PATH})

add_compile_options(-std=c++11 -O3 -Wall -march=native -Wno-unknown-pragmas)

# Use imported targets would make things much eazier. Thanks Levi for pointing it out.
add_executable(main ./main.cpp ${CPP_SOURCES} ${H_SOURCES})
target_link_libraries(main
  PRIVATE
  MPI_C)

