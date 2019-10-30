cmake_minimum_required (VERSION 2.6)
project (min_convex_decomp)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_STANDARD 17)

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall" )

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -W -Wfloat-equal -Wundef -Wpointer-arith" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wwrite-strings" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wredundant-decls -Wchar-subscripts -Wcomment -Wformat=2" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wwrite-strings -Wmissing-declarations -Wredundant-decls" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wswitch-enum" )

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Winit-self -Wmissing-field-initializers" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Waddress -Wmissing-noreturn -Wstrict-overflow=1" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wextra -Warray-bounds" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pedantic" )

if (POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif ()


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
