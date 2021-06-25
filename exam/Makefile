CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: default clean

default: out.txt plot.png
	cat $<

plot.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Cubic (sub-)spline interpolation";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [-4:4];\
	set yrange [-1.25:1.25];\
	plot \
	"$<" index 0 using 1:2 with points title "points",\
	"$<" index 1 using 1:3 with line color black title "generating function",\
	"$<" index 1 using 1:2 with line color red title "akima spline",\
	"$<" index 1 using 1:4 with line color green title "regular cubic spline"\
	' | pyxplot

out.txt error.txt data.txt:main
	./$< 1>out.txt 2>error.txt

main: main.o aspline.o cspline.o

clean:
	$(RM) main *.o *.txt *.png
