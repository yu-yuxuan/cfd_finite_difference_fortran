OBJECTS_2D_FTBS = boundary_conditions.o scheme.o initialize.o linear_advection_equation.o


FFLAGS = -g

.PHONY:   clean clean_doc test

linear_advection_equation: linear_advection_equation.exe
	./linear_advection_equation.exe
linear_advection_equation.exe: $(OBJECTS_2D_FTBS)
	gfortran $(FFLAGS) $(OBJECTS_2D_FTBS) -o linear_advection_equation.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
