OBJECTS_2D_BG = boundary_conditions.o scheme.o initialize.o two_d_burgers_ftbs.o


FFLAGS = -g

.PHONY:   clean clean_doc test

two_d_burgers_ftbs: two_d_burgers_ftbs.exe
	./two_d_burgers_ftbs.exe
two_d_burgers_ftbs.exe: $(OBJECTS_2D_BG)
	gfortran $(FFLAGS) $(OBJECTS_2D_BG) -o two_d_burgers_ftbs.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
