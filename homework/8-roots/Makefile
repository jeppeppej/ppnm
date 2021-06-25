CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default: out.txt errors.txt wavefunction.png
	cat out.txt

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt

main: functions.o

wavefunction.png: odeoutput.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Radial wavefunction";\
	set key top right;\
	set tics out;\
	set xlabel "y";\
	set ylabel "x";\
	set xrange [0:8];\
	set yrange [0:0.4];\
	plot \
	 "$<" index 0 using 1:2 with line color black title "Calculated",\
	 "$<" index 1 using 1:2 with line color red title "Exact"\
	' | pyxplot



odeoutput.txt: main

.PHONEY: default clean

clean:
	$(RM) errors.txt *.o *.txt *.png main
