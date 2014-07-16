CC = g++
CSWITCHES = -g -O3 -march=native -fPIC

RM = /bin/rm
# The action starts here.

all: libvoronoi.a test

test: voronoi_main.cpp voronoi.o
	$(CC) $(CSWITCHES) -o test voronoi.o  voronoi_main.cpp

libvoronoi.a: voronoi.o
	ar rcs libvoronoi.a voronoi.o

voronoi.o: VoronoiDiagramGenerator.cpp VoronoiDiagramGenerator.h
	$(CC) $(CSWITCHES) -c -o voronoi.o VoronoiDiagramGenerator.cpp

distclean:
	$(RM) libvoronoi.a voronoi.o test
