COMPILER=g++
FLAGS=-pedantic -lblas -llapacke -Wall -g
OBJDIR=objects/

all: hw3

hw3: $(OBJDIR)hw3.o
	$(COMPILER) -o hw3 $(OBJDIR)hw3.o $(FLAGS)

$(OBJDIR)hw3.o: hw3.cpp
	$(COMPILER) -c hw3.cpp $(FLAGS)
	mv hw3.o $(OBJDIR)

clean:
	rm $(OBJDIR)/*
	rm hw3
