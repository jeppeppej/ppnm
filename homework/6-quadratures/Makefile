CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: out.txt

out.txt: main
	./$< > $@
	cat out.txt

main: main.o linAlg.o linAlg.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c linAlg.h
	$(CC) $(CFLAGS) -c $< -o $@

linAlg.o: linAlg.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) main.o main functions.o out.txt *.png
