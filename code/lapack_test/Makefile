COMPILER=g++
FLAGS=-llapacke -lblas -pedantic

test: test.o
	$(COMPILER) -o test test.o $(FLAGS)

test.o: test.cpp
	$(COMPILER) -c test.cpp $(FLAGS)

clean:
	rm *.o test
