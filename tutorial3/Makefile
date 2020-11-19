all: Burgers1D debug sanitizer profiler

Burgers1D: Burgers1D.cpp
	g++ -g -DNDEBUG -Wall -o $@ $<

debug: Burgers1D.cpp
	g++ -g -Wall -o $@ $<

sanitizer: Burgers1D.cpp
	g++ -g -DNDEBUG -Wall -fsanitize=address -o $@ $<

profiler: Burgers1D.cpp
	g++ -g -DNDEBUG -pg -Wall -o $@ $<

clean:
	rm -f Burgers1D debug sanitizer profiler example.txt *.png *.out
