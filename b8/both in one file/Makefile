# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O0
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-par partdiff-par-hybrid 
PAR-OBJS = partdiff-par.o askparams.o
HYBRID-OBJS = partdiff-par-hybrid.o askparams.o

# Targets ...
all: partdiff-par partdiff-par-hybrid

partdiff-par: $(PAR-OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(PAR-OBJS) $(LIBS)

partdiff-par.o: partdiff-par.c Makefile

partdiff-par-hybrid: $(HYBRID-OBJS) Makefile
	$(CC) $(LFLAGS) -fopenmp -DOPENMP -o $@ $(HYBRID-OBJS) $(LIBS)

partdiff-par-hybrid.o: partdiff-par.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

partdiff-par-hybrid.o: partdiff-par.c
	$(CC) -c $(CFLAGS) -fopenmp -DOPENMP partdiff-par.c


clean:
	$(RM) $(PAR-OBJS)
	$(RM) $(HYBRID-OBJS)
	$(RM) $(TGTS)
