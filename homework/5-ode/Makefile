CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default : SIR.png

SIR.png: s.txt i.txt r.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set title "SIR-model, recovery time=14 days";\
	set key outside;\
	set tics out;\
	set xlabel "t (days)";\
	set ylabel "Population";\
	set xrange [0:50];\
	set yrange [0:60000000];\
	plot \
	 "s.txt" using 1:2 with line color black title "S(contact time=1)",\
	 "i.txt" using 1:2 with line color red title "I",\
	 "r.txt" using 1:2 with line color green title "R",\
	 "s.txt" using 1:3 with line lt 2 color black title "S(contact time=2)",\
	 "i.txt" using 1:3 with line lt 2 color red title "I",\
	 "r.txt" using 1:3 with line lt 2 color green title "R",\
	 "s.txt" using 1:4 with line lt 3 color black title "S(contact time=4)",\
	 "i.txt" using 1:4 with line lt 3 color red title "I",\
	 "r.txt" using 1:4 with line lt 3 color green title "R"\
	' | pyxplot

s.txt i.txt r.txt : main
	@./$< > $@

main : main.o

.PHONEY : clean

clean :
	@$(RM) main *.o *.txt *.png
