var = pointer map
all: $(var)
$(all): %.o
	gcc -o $@ $^
$(all): %.o: %.c
	gcc -c $^
clean:
	rm -f $(var)
	rm -f *.o
	rm -f *~