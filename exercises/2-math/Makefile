CC = gcc # the C compiler
CFLAGS = -O -std=gnu11 -Wall # options for the C compiler
LDLIBS = -lm # libraries to link

default: out.txt    # default target is to make out.txt
	cat out.txt # show out.txt on the screen

out.txt: math            # out.txt needs math
	./math > out.txt # run math, output goes into out.txt

math: math.o                           # math needs math.o
	$(CC) -o math math.o $(LDLIBS) # link math.o into math

math.o: math.c                   # math.o needs math.c
	$(CC) $(CFLAGS) -c math.c # compile math.c

.PHONEY: clean
clean:                              # clean is a phoney target
	$(RM) math.o math out.txt # clean up the directory

.PHONEY: test
test:                      # test target used for debugging
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
