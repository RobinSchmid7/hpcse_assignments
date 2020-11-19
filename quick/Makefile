CC=g++
CFLAGS=-O3 -Wall -Wextra -Wpedantic

all: calculation

calculation.o: calculation.cpp
	$(CC) -c -o calculation.o calculation.cpp $(CFLAGS)

calculation: calculation.o
	$(CC) -o calculation calculation.o

clean:
	rm -f *.o *~ calculation

.PHONY: all
.PHONY: clean
