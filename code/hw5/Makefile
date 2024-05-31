COMPILER=g++
FLAGS=-pedantic -lblas -llapacke -lfftw3 -lm -Wall -g
OBJDIR=objects/

all: hw5_pseudo hw5_coll

hw5_pseudo: $(OBJDIR)hw5_pseudo.o
	$(COMPILER) -o hw5_pseudo $(OBJDIR)hw5_pseudo.o $(FLAGS)

$(OBJDIR)hw5_pseudo.o: hw5_pseudo.cpp
	$(COMPILER) -c hw5_pseudo.cpp $(FLAGS)
	mv hw5_pseudo.o $(OBJDIR)

hw5_coll: $(OBJDIR)hw5_coll.o
	$(COMPILER) -o hw5_coll $(OBJDIR)hw5_coll.o $(FLAGS)

$(OBJDIR)hw5_coll.o: hw5_coll.cpp
	$(COMPILER) -c hw5_coll.cpp $(FLAGS)
	mv hw5_coll.o $(OBJDIR)

clean:
	rm $(OBJDIR)/*
	rm hw5_pseudo
	rm hw5_coll
