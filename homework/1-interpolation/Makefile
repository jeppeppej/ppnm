CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: fig-linear.png fig-quadratic.png fig-cubic.png

fig-linear.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Linear spline with deriv and int";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:10];\
	set yrange [-10:30];\
	plot \
	"points.txt" using 1:2 with points title "points",\
	"$<" index 0 using 1:2 with line color red title "linear spline",\
	"$<" index 0 using 1:3 with line color red title "integral of linear spline",\
	"$<" index 0 using 1:4 with line color blue title "GSLs linear spline",\
	"$<" index 0 using 1:5 with line color blue title "GSLs integral of spline"\
	' | pyxplot

fig-quadratic.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Quadratic spline with deriv and int";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:10];\
	set yrange [-10:30];\
	plot \
	"points.txt" using 1:2 with points title "points",\
	 "$<" index 1 using 1:2 with line color black title "Quadratic spline",\
	 "$<" index 1 using 1:3 with line color red title "Integral of spline",\
	 "$<" index 1 using 1:4 with line color green title "Derivative of spline"\
	' | pyxplot

fig-cubic.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Cubic spline with deriv and int";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:10];\
	set yrange [-10:30];\
	plot \
	"points.txt" using 1:2 with points title "points",\
	 "$<" index 2 using 1:2 with line color black title "Cubic spline",\
	 "$<" index 2 using 1:3 with line color red title "Integral of spline",\
	 "$<" index 2 using 1:4 with line color green title "Derivaite spline"\
	' | pyxplot

out.txt: main
	./$< > $@

main: main.o linear_spline.o quad_spline.o cubic_spline.o

clean:
	$(RM) ms.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex *.png *.out
