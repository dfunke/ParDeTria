set(CMAKE_C_FLAGS "-g -O3 -march=native -fPIC")

add_library(delaunay delaunay.c delaunay.h predicates.c)
