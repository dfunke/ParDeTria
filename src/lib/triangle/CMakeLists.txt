set(CMAKE_CXX_FLAGS "-g -O3 -march=native -DTRILIBRARY -DANSI_DECLARATORS -DLINUX -DNO_TIMER -I/usr/X11R6/include -L/usr/X11R6/lib -fPIC")

add_library(triangle triangle.cpp triangle.h)
