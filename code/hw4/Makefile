COMPILER=g++
FLAGS=-pedantic -lblas -llapacke -lfftw3 -lm -Wall -g
OBJDIR=objects/

all: hw4

hw4: $(OBJDIR)hw4.o
	$(COMPILER) -o hw4 $(OBJDIR)hw4.o $(FLAGS)

$(OBJDIR)hw4.o: hw4.cpp
	$(COMPILER) -c hw4.cpp $(FLAGS)
	mv hw4.o $(OBJDIR)

clean:
	rm $(OBJDIR)/*
	rm hw4
