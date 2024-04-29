CC=g++
CFLAGS=-pedantic -lblas -llapacke
OBJDIR=objects/

out: $(OBJDIR)hw1.o $(OBJDIR)hw1_elem.o $(OBJDIR)utils_205.o
	$(CC) -o out $(OBJDIR)hw1.o $(OBJDIR)utils_205.o $(CFLAGS)

$(OBJDIR)hw1.o: hw1.cpp
	$(CC) -c hw1.cpp $(CFLAGS)
	mv hw1.o $(OBJDIR)

$(OBJDIR)hw1_elem.o: hw1_elem.cpp
	$(CC) -c hw1_elem.cpp $(CFLAGS)
	mv hw1_elem.o $(OBJDIR)

$(OBJDIR)utils_205.o: utils_205.cpp
	$(CC) -c utils_205.cpp $(CFLAGS)
	mv utils_205.o $(OBJDIR)

clean:
	rm *.o out
