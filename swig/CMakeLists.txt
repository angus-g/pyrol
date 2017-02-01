# This is a CMake example for Python

cmake_minimum_required(VERSION 2.8.9)
project(ROL)
set(CMAKE_BUILD_TYPE Release)

message(STATUS "Checking for Trilinos")
find_package(Trilinos PATHS ${TRILINOS_DIR} ${Trilinos_DIR}
  $ENV{TRILINOS_DIR})

if ("${Trilinos_VERSION}" VERSION_LESS "12.10.0")
  set(Trilinos_FOUND FALSE)
  message(FATAL_ERROR "Unable to find Trilinos (>= 12.10.0)")
endif()
message(STATUS "Found Trilinos")

find_library(ROL PATHS ${Trilinos_DIR})
message(STATUS "Found ROL: " ${ROL_INCLUDE_DIRS})

find_package(PythonLibs)
find_package(MPI)

#Bring the headers, such as Student.h into the project
include_directories(${PYTHON_INCLUDE_DIRS} ${ROL_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})
add_compile_options(-std=c++11)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# SET(CMAKE_SWIG_FLAGS "-python")

SET_SOURCE_FILES_PROPERTIES(ROL.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(example.i PROPERTIES SWIG_FLAGS "")
SWIG_ADD_MODULE(ROL python ROL.i)
SWIG_LINK_LIBRARIES(ROL ${PYTHON_LIBRARIES} ${MPI_LIBRARY} ${ROL_LIBRARIES})


#Generate the shared library from the sources
# add_library(_ROL SHARED ROL_wrap.cxx)
# target_link_libraries(_ROL ${MPI_LIBRARY} ${ROL_LIBRARIES})
# set_target_properties(_ROL PROPERTIES PREFIX "")
