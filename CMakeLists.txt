cmake_minimum_required(VERSION 3.14)
project(_ROL VERSION 0.3.0)

# minimum cpp for using filesystem 
set(CMAKE_CXX_STANDARD 17)


set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# making sure rpath is set correctly
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
   set(MACOSX_RPATH TRUE)
else()
   set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()

# options for Trilinos
set(BUILD_SHARED_LIBS ON CACHE BOOL "")
set(TPL_ENABLE_MPI OFF CACHE BOOL "")
set(Trilinos_INSTALL_LIBRARIES_AND_HEADERS OFF CACHE BOOL "")
set(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "")
set(Trilinos_ENABLE_TESTS OFF CACHE BOOL "")
set(Trilinos_ENABLE_ROL ON CACHE BOOL "")
set(Trilinos_ENABLE_EXPLICIT_INSTANTIATION ON CACHE BOOL "")
set(ROL_ParameterList boost::property_tree CACHE STRING "")
set(ROL_Ptr std::shared_ptr CACHE STRING "")
set(ROL_stacktrace noop CACHE STRING "")
set(ROL_ENABLE_EXAMPLES OFF CACHE BOOL "")
set(ROL_ENABLE_TESTS OFF CACHE BOOL "")
set(ROL_ENABLE_Sacado OFF CACHE BOOL "")
set(TPL_BLAS_LIBRARIES "" CACHE STRING "")
set(TPL_LAPACK_LIBRARIES "" CACHE STRING "")


include(FetchContent)
FetchContent_Declare(Trilinos
  GIT_REPOSITORY https://github.com/angus-g/Trilinos
  GIT_SHALLOW TRUE
  GIT_TAG "pyrol-2.0-checkpointing")
FetchContent_MakeAvailable(Trilinos)

find_package(Boost)

add_subdirectory(pybind11)
pybind11_add_module(_ROL MODULE ROL/ROL.cpp
  ROL/algorithmstate.cpp
  ROL/augmentedlagrangian.cpp
  ROL/bounds.cpp
  ROL/constraint.cpp
  ROL/objective.cpp
  ROL/statustest.cpp
  ROL/optimizationproblem.cpp
  ROL/optimizationsolver.cpp
  ROL/secant.cpp
  # ROL/step.cpp
  ROL/algorithm.cpp
  ROL/problem.cpp
  ROL/solver.cpp
  ROL/parameterlist.cpp
  ROL/stdvector.cpp
  ROL/vector.cpp)

# Testing if we have filesystem
try_compile(HAS_FS "${CMAKE_BINARY_DIR}/temp" 
"${CMAKE_SOURCE_DIR}/tests/has_filesystem.cc" 
            CMAKE_FLAGS -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON
            LINK_LIBRARIES stdc++fs)

target_include_directories(_ROL PRIVATE ${ROL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} include)
target_link_libraries(_ROL PRIVATE rol)
if (HAS_FS)
   target_link_libraries(_ROL PRIVATE -lstdc++fs)
endif(HAS_FS)
set_target_properties(_ROL PROPERTIES COMPILE_DEFINITIONS "ENABLE_PYROL")



