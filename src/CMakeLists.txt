cmake_minimum_required(VERSION 2.8)

project (min_convex_decomp)

set( min_convex_decomp_SRC
	io.cc arg_eval.cc onion.cc decomp.cc data.cpp
)

include(CheckFunctionExists) 
check_function_exists(random HAVE_RANDOM)

if( NOT HAVE_RANDOM )
	message("Adding custom random function")
	add_definitions(-DRAND)
endif( NOT HAVE_RANDOM )

add_executable(mcd_apx main.cc ${min_convex_decomp_SRC}) 
target_link_libraries(mcd_apx ${min_convex_decomp_LIBS})


get_directory_property( DirDefs COMPILE_DEFINITIONS ) 
message( "COMPILE_DEFINITIONS = ${DirDefs}" )
