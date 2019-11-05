project (mcd-apx VERSION "1.99" DESCRIPTION "mcd-apx")

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
