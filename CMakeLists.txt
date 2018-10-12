# CMake entry point
cmake_minimum_required (VERSION 2.6)

project (feature-mser)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-g -Wall -Wextra -O2 -std=c++14")
set(CMAKE_C_FLAGS "-Wall -O2")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

cmake_policy(SET CMP0054 NEW) #FindOpenMP.cmake
#find_package(OpenMP)
#    if (OPENMP_FOUND)
#        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#    endif()

message( "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}" )
message( "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}" )
message( "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}" )
message( "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}" )

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory." )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
    message( "Your Source Directory contains spaces." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
    message( "Your Build Directory contains spaces." )
endif()

# external libs
set(klib "external/klib")
INCLUDE_DIRECTORIES(
    src
    ${klib}
)

message( "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}" )

ADD_EXECUTABLE(feature-mser
    src/utils.c
    src/main.c
)

TARGET_LINK_LIBRARIES(feature-mser
    ${ALL_LIBS}
)
