FC = gfortran
CC=g++
CFLAGS=-pedantic -lblas -llapacke
OBJDIR=objects/
OBJs=105.o fortran_utils.o dsolve.o
FFLAGS=-ffree-line-length-none -ffixed-line-length-none

all: bem fd

bem: $(OBJDIR)hw1.o $(OBJDIR)hw1_elem.o $(OBJDIR)utils_205.o
	$(CC) -o bem $(OBJDIR)hw1.o $(OBJDIR)utils_205.o $(CFLAGS)

fd: $(OBJDIR)105.o $(OBJDIR)dsolve.o $(OBJDIR)fortran_utils.o
	$(FC) -o fd $(OBJDIR)105.o $(OBJDIR)dsolve.o $(OBJDIR)fortran_utils.o $(FFLAGS)

$(OBJDIR)hw1.o: hw1.cpp
	$(CC) -c hw1.cpp $(CFLAGS)
	mv hw1.o $(OBJDIR)

$(OBJDIR)hw1_elem.o: hw1_elem.cpp
	$(CC) -c hw1_elem.cpp $(CFLAGS)
	mv hw1_elem.o $(OBJDIR)

$(OBJDIR)utils_205.o: utils_205.cpp
	$(CC) -c utils_205.cpp $(CFLAGS)
	mv utils_205.o $(OBJDIR)

$(OBJDIR)105.o: 105.f
	$(FC) -c 105.f $(FFLAGS)
	mv 105.o $(OBJDIR)

$(OBJDIR)dsolve.o: dsolve.f
	$(FC) -c dsolve.f $(FFLAGS)
	mv dsolve.o $(OBJDIR)

$(OBJDIR)fortran_utils.o: dsolve.f
	$(FC) -c fortran_utils.f $(FFLAGS)
	mv fortran_utils.o $(OBJDIR)

clean:
	rm -rf *.o *.mod fd bem
	rm $(OBJDIR)*
