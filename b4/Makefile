# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -fopenmp -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff
OBJS = partdiff.o askparams.o

# Targets ...
all: partdiff

partdiff: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-openmp: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-openmp-element: partdiff-openmp-element.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-openmp-element.o askparams.o $(LIBS)

partdiff-openmp-zeilen: partdiff-openmp-zeilen.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-openmp-zeilen.o askparams.o $(LIBS)

partdiff-openmp-spalten: partdiff-openmp-spalten.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-openmp-spalten.o askparams.o $(LIBS)

partdiff.o: partdiff.c Makefile

partdiff-openmp-element.o: partdiff-openmp-element.c Makefile
partdiff-openmp-zeilen.o: partdiff-openmp-zeilen.c Makefile
partdiff-openmp-spalten.o: partdiff-openmp-spalten.c Makefile

askparams.o: askparams.c Makefile


# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)









