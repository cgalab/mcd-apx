cmake_minimum_required (VERSION 2.6)
project (min_convex_decomp)

set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-gdwarf-2 -g3 -DDEBUG -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# warning flags for the C++ compiler
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -pedantic" )

if( WITH_DEBUGSYMBOLS )
   message("Adding debugging symbols")
   set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g3 -gdwarf-2")
endif()

set(PACKAGE_STRING min_convex_decomp)
set(PACKAGE_BUGREPORT held@cs.sbg.ac.at)
set(PACKAGE_VERSION 1.0)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
# More verbose output by make
set(CMAKE_VERBOSE_MAKEFILE "OFF")
add_subdirectory(src)
