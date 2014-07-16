CC = cc
CSWITCHES = -g -O3 -march=native -fPIC

RM = /bin/rm
# The action starts here.

all: libdelaunay.a

libdelaunay.a: delaunay.o predicates.o
	ar rcs libdelaunay.a delaunay.o predicates.o

delaunay.o: delaunay.c predicates.o
	$(CC) $(CSWITCHES) -c -o delaunay.o delaunay.c

predicates.o: predicates.c
	$(CC) $(CSWITCHES) -c -o predicates.o predicates.c

distclean:
	$(RM) libdelaunay.a delaunay.o predicates.o
