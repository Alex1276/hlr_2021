all: pointer map
$all: %.o
	gcc -o $@ $^
pointer.o: pointer.c
	gcc -c pointer.c
map.o: map.c
	gcc -c map.c
clean:
	rm -f map
	rm -f pointer
	rm -f *.o
	rm -f *~