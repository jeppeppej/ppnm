CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

ms.pdf: ms.tex  fig-pyxplot.pdf #fig-gpl.tex
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<

out.data.txt: main
	./$< > $@

main: main.o myexp.o


#fig-gpl.tex : fig-gpl.gpi out.data.txt
#	gnuplot $<


fig-pyxplot.pdf: fig-pyxplot.gpi out.data.txt
	pyxplot $<

clean:
	$(RM) ms.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex
