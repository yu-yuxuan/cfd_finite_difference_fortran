OBJECTS_2D_LP = boundary_conditions.o scheme.o initialize.o two_d_laplace.o


FFLAGS = -g

.PHONY:   clean clean_doc test

two_d_laplace: two_d_laplace.exe
	./two_d_laplace.exe
two_d_laplace.exe: $(OBJECTS_2D_LP)
	gfortran $(FFLAGS) $(OBJECTS_2D_LP) -o two_d_laplace.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
