CFLAGS = -O -Wall $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default clean

default: out.txt
	cat out.txt

out.txt: main
	./$< > $@

main: main.o
main.o: main.c


clean:
	$(RM) main *.o *.dat *.png
