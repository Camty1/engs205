COMPILER=g++
FLAGS=-pedantic -lblas -llapacke
OBJDIR=objects/

all: hw1 hw2

hw1: $(OBJDIR)hw1.o $(OBJDIR)utils_205.o
	$(COMPILER) -o hw1 $(OBJDIR)hw1.o $(OBJDIR)utils_205.o $(FLAGS)

hw2: $(OBJDIR)hw2.o $(OBJDIR)utils_205.o
	$(COMPILER) -o hw2 $(OBJDIR)hw2.o $(OBJDIR)utils_205.o $(FLAGS)

$(OBJDIR)hw1.o: hw1.cpp
	$(COMPILER) -c hw1.cpp $(FLAGS)
	mv hw1.o $(OBJDIR)

$(OBJDIR)hw2.o: hw2.cpp
	$(COMPILER) -c hw2.cpp $(FLAGS)
	mv hw2.o $(OBJDIR)

$(OBJDIR)utils_205.o: utils_205.cpp
	$(COMPILER) -c utils_205.cpp $(FLAGS)
	mv utils_205.o $(OBJDIR)

clean:
	rm $(OBJDIR)*.o hw1 hw2
