all: timempi timempi2
	
timempi: timempi.o
	mpicc -o timempi timempi.o

timempi.o: timempi.c
	mpicc -c timempi.c

timempi2: timempi2.o
	mpicc -o timempi2 timempi2.o

timempi2.o: timempi2.c
	mpicc -c timempi2.c

clean:
	rm -f timempi
	rm -f timempi2
	rm -f *.o
	rm -f *~