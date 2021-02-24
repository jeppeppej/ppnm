CC = gcc # a driver that runs linker and compileren. Bruges som variabel, så man nemt kan udskifte.
CFLAGS = -gnu99 -0 -Wall # nogle indstillinger, -Wall er all warnings
LDLIBS = -lm # mathemathical library - standard for C is called M

out.txt: hello  # kør programmet hello i filen out.txt
	./hello > out.txt	# ./ bruges fordi man bruger program fra egen folder. 
				#Virker overflødigt, men vigtigt at gøre alting uafhængigt af terminalen.
hell.o: hello.c hello.h
	$(CC) $(CFLAGS) -c hello.c -o hello.o

hello: hello.o world.o	# der er intet link til libraries, det gøres her.
	$(CC) $(LDFLAGS) hello.o world.o -o hello $(LDLIBS) # long options kom senere, de korte virker stadig. Libraries skal helst komme til sidst i denne command.

world.o: world.c hello.h # good practice
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o out* hello
