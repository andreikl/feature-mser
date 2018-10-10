# CMake entry point
cmake_minimum_required (VERSION 2.6)

project (feature-mser)

if(POLICY CMP0026)
    cmake_policy(SET CMP0026 OLD)
endif()

if(POLICY CMP0012)
    cmake_policy(SET CMP0012 NEW)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})

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

include_directories(
     src
)

add_executable(feature-mser
    src/main.c
)

target_link_libraries(feature-mser
    ${ALL_LIBS}
)
