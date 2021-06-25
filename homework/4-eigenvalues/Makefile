CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt plot.png
	cat out.txt

out.txt: main
	./$< > $@

main: jacobi_diag.o

plot.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Eigenfunctions to particle in a box";\
	set ylabel "Wavefunction";\
	set xlabel "x";\
	plot \
	"$<" using 1:2 with lines title "psi0",\
	"$<" using 1:3 with lines title "psi1",\
	"$<" using 1:4 with lines title "psi2"\
	'|pyxplot

data.txt: main

clean:
	$(RM) main *.txt
