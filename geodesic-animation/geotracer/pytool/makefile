
CC = gcc
CFLAGS = -Wall -Werror -fpic -O3

geotrace.so: libagg.o
	$(CC) -shared -o libgeotrace.so mathfunc.o geometry.o model.o libgeotrace.o 
	-rm *.o

libagg.o: libgeotrace.o mathfunc.o model.o geometry.o 
	ld -r mathfunc.o geometry.o model.o libgeotrace.o -o libagg.o

clean:
	-rm libgeotrace.so
	-rm *.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
