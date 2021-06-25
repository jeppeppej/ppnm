CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


default: plot.png

plot.png: out.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set key outside;\
	set title "ANN interpolation";\
	f(x) = x*exp(-x);\
	plot \
	"$<" index 0 using 1:2 with points title "points"\
	,"$<" index 1 using 1:2 with line color black title "ANN spline"\
	,"$<" index 1 using 1:3 with line color red title "ANN integral"\
	,"$<" index 1 using 1:4 with line color green title "ANN derivative"\
	' | pyxplot

out.txt: main
	./$< > $@

main: main.c ann.c functions.c

clean:
	$(RM) out* main *.png
