all: circle

circle: circle.c Makefile
	scorep mpicc -std=c99 -ggdb -Wall -Wextra circle.c -o circle

clean:
	rm -f circle 