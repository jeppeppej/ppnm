CFLAGS = -O -Wall -Werror  -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp -lgomp
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm -pthread
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

clean:
	$(RM) *.txt *.o main
