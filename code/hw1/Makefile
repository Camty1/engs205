COMPILER=g++
FLAGS=-pedantic -lblas -llapacke

out: hw1.o hw1_redo.o utils_205.o
	$(COMPILER) -o out hw1_redo.o utils_205.o $(FLAGS)

hw1.o: hw1.cpp
	$(COMPILER) -c hw1.cpp $(FLAGS)

hw1_redo.o: hw1_redo.cpp
	$(COMPILER) -c hw1_redo.cpp $(FLAGS)

utils_205.o: utils_205.cpp
	$(COMPILER) -c utils_205.cpp $(FLAGS)

clean:
	rm *.o out
