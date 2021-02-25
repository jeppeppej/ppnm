CFLAGS = -O -Wall -std=gnu11
#CFLAGS += `gsl-config --cflags`
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: erf.pyxplot.png gamma.pyxplot.png

erf.pyxplot.png: error.txt Makefile input.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "myerf"\
	,"input.data.txt" using 1:2 with points pointtype 2 title "tabulated data"\'| pyxplot


gamma.pyxplot.png: gamma.txt Makefile tabulatedGamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "tgamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "mygamma"\
	,"tabulatedGamma.txt" using 1:2 with points pointtype 2 title "tabulated data"\'| pyxplot


error.txt gamma.txt: main
	./$< > $@

main: main.o

.PHONEY:clean
clean:
	$(RM) *.o main *.png error.txt gamma.txt
