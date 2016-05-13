OBJECTS_1D_BG = boundary_conditions.o scheme.o initialize.o one_d_burger_FTBS.o


FFLAGS = -g

.PHONY:   clean clean_doc test

one_d_burger_FTBS: one_d_burger_FTBS.exe
	./one_d_burger_FTBS.exe
one_d_burger_FTBS.exe: $(OBJECTS_1D_BG)
	gfortran $(FFLAGS) $(OBJECTS_1D_BG) -o one_d_burger_FTBS.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
