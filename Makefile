CC=g++
CFLAGS=-std=c++11 -pthread -g

common = common.h tools/debug.h
file = file.o file.h
myGraph = myGraph.o myGraph.h
decomp = decomp.o decomp.h
myAlgorithm = myAlgorithm.o myAlgorithm.h
main: main.o $(common) $(file) $(myGraph) $(decomp) $(myAlgorithm)
	$(CC) $(CFLAGS) -rdynamic -o main.out main.o file.o myGraph.o decomp.o myAlgorithm.o 
main.o :main.cpp $(common) file.h myGraph.h decomp.h
	$(CC) $(CFLAGS) -rdynamic -c main.cpp
file.o: file.cpp $(common) file.h myGraph.h
	$(CC) $(CFLAGS) -rdynamic -c file.cpp
myGraph.o: myGraph.cpp myGraph.h $(common)
	$(CC) $(CFLAGS) -rdynamic -c myGraph.cpp
decomp.o: decomp.cpp decomp.h $(common) myGraph.h
	$(CC) $(CFLAGS) -rdynamic -c decomp.cpp
myAlgorithm.o: myAlgorithm.cpp  myAlgorithm.h $(common)  myGraph.h 
	$(CC) $(CFLAGS) -rdynamic -c myAlgorithm.cpp

clean:
	rm file.o
