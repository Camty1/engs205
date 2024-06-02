COMPILER=g++
FLAGS=-pedantic -lblas -llapacke -lfftw3 -lm -Wall -g
OBJDIR=objects/

all: project

project: $(OBJDIR)project.o
	$(COMPILER) -o project $(OBJDIR)project.o $(FLAGS)

$(OBJDIR)project.o: project.cpp
	$(COMPILER) -c project.cpp $(FLAGS)
	mv project.o $(OBJDIR)

clean:
	rm $(OBJDIR)/*
	rm project
