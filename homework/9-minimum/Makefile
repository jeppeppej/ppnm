CFLAGS = -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


default: out.txt wave.png
	cat out.txt

out.txt: main
	./$< > $@

wave.png: plot.txt data.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "Energy(GeV)";\
	set ylabel "Cross Section";\
	set title "Higgs fit";\
	plot \
	"plot.txt" using 1:2 with lines color blue title "fit" \
	,"data.txt" using 1:2 with points color red title "data" \
	' | pyxplot

main: main.c functions.c testfunctions.c

.PHONEY:clean
clean:
	$(RM) out.txt plot.txt main *png
