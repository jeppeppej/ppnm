CFLAGS = -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: sqrt.png given.png
	cat out.txt

sqrt.png: main Makefile plainSqrt.txt quasiSqrt.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Error";\
	set title "Plain MC and using quasi-random sequences on sqrt(x)";\
	set logscale y;\
	a=1; b=-0.5; f(x)=a*x**b ;\
	c=1; d=-1; g(x)=c*x**d ;\
	fit f(x) withouterrors "plainSqrt.txt" using 1:3 via a,b ;\
	fit g(x) withouterrors "quasiSqrt.txt" using 1:3 via c,d ;\
	plot \
	"plainSqrt.txt" using 1:3 with points title "Plain error ",\
	"quasiSqrt.txt" using 1:3 with points title "Quasi-random error",\
	 f(x) with line title "$$%.2fn^{%.2f}$$"%(a,b),\
	 g(x) with line title "$$%.2fn^{%.2f}$$"%(c,d)\
	' |pyxplot

given.png: main Makefile plainGiven.txt quasiGiven.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Plain MC and using quasi-random sequences on a given integrand";\
	set key outside;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Error";\
	set logscale y;\
	e=1;f=-0.5;\
	i(x)=e*x**f;\
	fit i(x) withouterrors "plainGiven.txt" using 1:3 via e,f;\
	g=1;h=-1;\
	j(x)=g*x**h;\
	fit j(x) withouterrors "quasiGiven.txt" using 1:3 via g,h;\
	plot \
	 "plainGiven.txt" using 1:3 with points title "pseudo",\
	 "quasiGiven.txt" using 1:3 with points title "quasi",\
	 i(x) with line title "$$%.2fn^{%.2f}$$"%(e,f),\
	 j(x) with line title "$$%.2fn^{%.2f}$$"%(g,h)\
	' | pyxplot


plainSqrt.txt:main

quasiSqrt.txt:main

plainGiven.txt:main

qusaiGiven.txt:main

out.txt: main
	./$< > $@

main: functions.c main.c

clean:
	$(RM) main main.o functions.o *.txt *.png *.gpi *.log
