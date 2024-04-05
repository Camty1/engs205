COMPILER=g++
FLAGS=-pedantic -lblas -llapacke

out: hw1.o utils_205.o
	$(COMPILER) -o out hw1.o utils_205.o $(FLAGS)

hw1.o: hw1.cpp
	$(COMPILER) -c hw1.cpp $(FLAGS)

utils_205.o: utils_205.cpp
	$(COMPILER) -c utils_205.cpp $(FLAGS)

clean:
	rm *.o out
