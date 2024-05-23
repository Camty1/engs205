COMPILER=g++
FLAGS=-pedantic -lblas -llapacke -Wall -g
OBJDIR=objects/

all: hw3 simple

hw3: $(OBJDIR)hw3.o
	$(COMPILER) -o hw3 $(OBJDIR)hw3.o $(FLAGS)

simple: $(OBJDIR)simple.o
	$(COMPILER) -o simple $(OBJDIR)simple.o $(FLAGS)

$(OBJDIR)hw3.o: hw3.cpp
	$(COMPILER) -c hw3.cpp $(FLAGS)
	mv hw3.o $(OBJDIR)

$(OBJDIR)simple.o: simple.cpp
	$(COMPILER) -c simple.cpp $(FLAGS)
	mv simple.o $(OBJDIR)

clean:
	rm $(OBJDIR)/*
	rm hw3 simple
