CFLAGS = -Wall -std=gnu11 -O
LDLIBS = -lm


default: out.txt input.txt
	cat out.txt

out.txt: main
	./$< CMD $(shell cat input.txt) > $@
	./$< STDIN < input.txt >> $@
	./$< FILE input.txt $@

#input.txt: read
#	read > input.txt

.PHONEY:clean
clean:
	$(RM) main *.o out.txt
