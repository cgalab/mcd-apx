cmake_minimum_required(VERSION 2.8)

project (min_convex_decomp)

set( min_convex_decomp_SRC
	io.cc arg_eval.cc onion.cc decomp.cc
	"${PROJECT_BINARY_DIR}/gitversion.h"
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_executable(mcd_apx main.cc ${min_convex_decomp_SRC})
target_link_libraries(mcd_apx ${min_convex_decomp_LIBS})

add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/gitversion.h"
  DEPENDS "${CMAKE_SOURCE_DIR}/.git/logs/HEAD"
  COMMAND ${CMAKE_COMMAND}
    -Dlocal_dir="${CMAKE_SOURCE_DIR}"
    -Doutput_dir="${PROJECT_BINARY_DIR}"
    -P "${CMAKE_SOURCE_DIR}/cmake/gitversion.cmake"
)


get_directory_property( DirDefs COMPILE_DEFINITIONS )
message( "COMPILE_DEFINITIONS = ${DirDefs}" )
