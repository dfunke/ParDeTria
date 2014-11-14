set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

add_library(delaunay delaunay.c delaunay.h predicates.c)
