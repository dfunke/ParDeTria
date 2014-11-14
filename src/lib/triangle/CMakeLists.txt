set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRILIBRARY -DANSI_DECLARATORS -DLINUX -DNO_TIMER -fPIC")

add_library(triangle triangle.cpp triangle.h)
