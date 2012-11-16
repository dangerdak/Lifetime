lifetime: main.o lifetime.o
	g++ -g -Wall -static main.o lifetime.o -o lifetime -lgsl -lgslcblas -lm

main.o: main.cpp lifetime.h
	g++ -g -Wall -c main.cpp

lifetime.o: lifetime.cpp lifetime.h
	g++ -g -Wall -I/usr/local/include -c lifetime.cpp

clean:
	rm -rf *.o lifetime
