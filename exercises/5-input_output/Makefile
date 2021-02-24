CFLAGS = -Wall -std=gnu11 -O
LDLIBS = -lm


default: out.txt
	cat out.txt

out.txt: main
	./main $(cat input.txt) > out.txt

main: main.o
main.o: main.c


.PHONEY:clean
clean:
	$(RM) main *.o out.txt 
