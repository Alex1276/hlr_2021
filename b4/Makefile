# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff
OBJS = partdiff.o askparams.o

# Targets ...
all: partdiff

partdiff: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff.o: partdiff.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
	rm -f *.o
	rm -f partdiff-openmp
	rm -f partdiff-openmp-element
	rm -f partdiff-openmp-zeilen
	rm -f partdiff-openmp-spalten

partdiff-openmp:
	partdiff-openmp: $(OBJS) Makefile
	$(CC) $(CFLAGS) -fopenmp -o $@ $(OBJS) $(LIBS)
	partdiff.o: partdiff.c Makefile
	%.o: %.c
		$(CC) -c $(CFLAGS) -fopenmp $*.c

partdiff-openmp-element:
	NEWOBJS = partdiff-openmp-element.o askparams.o
	partdiff-openmp-element: $(NEWOBJS) Makefile
	$(CC) $(CFLAGS) -fopenmp -o $@ $(NEWOBJS) $(LIBS)
	partdiff.o: partdiff.c Makefile
	%.o: %.c
	  $(CC) -c $(CFLAGS) -fopenmp $*.c

partdiff-openmp-zeilen:
	NEWOBJS = partdiff-openmp-zeilen.o askparams.o
	partdiff-openmp-zeilen: $(NEWOBJS) Makefile
	$(CC) $(CFLAGS) -fopenmp -o $@ $(NEWOBJS) $(LIBS)
	partdiff.o: partdiff.c Makefile
	%.o: %.c
	  $(CC) -c $(CFLAGS) -fopenmp $*.c

partdiff-openmp-spalten:
	NEWOBJS = partdiff-openmp-spalten.o askparams.o
	partdiff-openmp-spalten: $(NEWOBJS) Makefile
	$(CC) $(CFLAGS) -fopenmp -o $@ $(NEWOBJS) $(LIBS)
	partdiff.o: partdiff.c Makefile
	%.o: %.c
		$(CC) -c $(CFLAGS) -fopenmp $*.c
