OBJECTS_2D_CF = boundary_conditions.o scheme.o initialize.o cavity_flow.o


FFLAGS = -g

.PHONY:   clean clean_doc test

cavity_flow: cavity_flow.exe
	./cavity_flow.exe
cavity_flow.exe: $(OBJECTS_2D_CF)
	gfortran $(FFLAGS) $(OBJECTS_2D_CF) -o cavity_flow.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
