cmake_minimum_required(VERSION 2.8.9)
project(ROL)
set(CMAKE_BUILD_TYPE Release)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "Checking for Trilinos")
find_package(Trilinos PATHS ${TRILINOS_DIR} ${Trilinos_DIR}
  $ENV{TRILINOS_DIR})

if ("${Trilinos_VERSION}" VERSION_LESS "12.10.0")
  set(Trilinos_FOUND FALSE)
  message(FATAL_ERROR "Unable to find Trilinos (>= 12.10.0)")
endif()
message(STATUS "Found Trilinos")

find_library(ROL PATHS ${Trilinos_DIR})
if (DEFINED ROL_LIBRARIES)
  set(ROL_FOUND TRUE)
  message(STATUS "Found ROL: " ${ROL_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "ROL library not found - make sure Trilinos is compiled with ROL support")
endif()

find_package(PythonLibs)
find_package(MPI)

find_path(EIGEN_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
  HINTS
  ${EIGEN_DIR}
  ENV EIGEN_DIR
  ENV EIGEN3_ROOT
  ENV EIGEN3_ROOT_DIR
  PATHS
  /usr/include/eigen3
  ${CMAKE_INSTALL_PREFIX}/include
  ${KDE4_INCLUDE_DIR}
  PATH_SUFFIXES eigen3 eigen
  )

#Bring the headers, such as Student.h into the project
include_directories(pybind11/include ${EIGEN_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS} ${ROL_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})

add_compile_options(-std=c++11 -g)

#Generate the shared library from the sources
add_library(ROL MODULE ROL.cpp)
target_link_libraries(ROL ${MPI_LIBRARY} ${ROL_LIBRARIES} ${PYTHON_LIBRARY})
set_target_properties(ROL PROPERTIES PREFIX "")

install(TARGETS ROL LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
