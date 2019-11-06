#project (min_convex_decomp VERSION "1.99" DESCRIPTION "min_convex_decomp")
cmake_minimum_required(VERSION 2.8)
project (min_convex_decomp)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(trianglelib
                       triangle.c
           )
set_target_properties(trianglelib PROPERTIES PUBLIC_HEADER triangle.h)
target_compile_definitions(trianglelib PRIVATE TRILIBRARY=1 ANSI_DECLARATORS=1)

add_executable(triangle
                       triangle.c
              )


target_link_libraries(triangle m)
target_link_libraries(trianglelib m)
