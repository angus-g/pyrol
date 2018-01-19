cmake_minimum_required(VERSION 2.8.9)
project(_ROL)
set(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
  message(STATUS "Found ROL: " ${ROL_LIBRARIES})
else()
  message(FATAL_ERROR "ROL library not found - make sure Trilinos is compiled with ROL support")
endif()

##find_package(PythonLibs)
##find_package(PythonInterp)
##message(STATUS "Python version found:" ${PYTHON_VERSION_STRING})
find_package(MPI)

add_subdirectory(pybind11)
#find_package(pybind11 REQUIRED) # if not including pybind11 in subdirectory

pybind11_add_module(_ROL MODULE ROL/ROL.cpp
  ROL/algorithm.cpp 
  ROL/algorithmstate.cpp 
  ROL/augmentedlagrangian.cpp 
  ROL/bounds.cpp 
  ROL/constraint.cpp 
  ROL/moreauyosidapenalty.cpp 
  ROL/objective.cpp 
  ROL/optimizationproblem.cpp 
  ROL/optimizationsolver.cpp 
  ROL/parameterlist.cpp 
  ROL/stdvector.cpp
  ROL/vector.cpp)

target_include_directories(_ROL PRIVATE ${ROL_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(_ROL PRIVATE ${MPI_LIBRARY} ${ROL_LIBRARIES})
set_target_properties(_ROL PROPERTIES COMPILE_DEFINITIONS "ENABLE_PYROL")
install(TARGETS _ROL DESTINATION ${CMAKE_INSTALL_PREFIX})
