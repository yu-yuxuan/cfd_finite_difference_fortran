OBJECTS_2D_FTBS = boundary_conditions.o scheme.o initialize.o two_d_advt_burgers_FTBS.o

FFLAGS = -g

.PHONY:   clean clean_doc test

two_d_advt_burgers_FTBS: two_d_advt_burgers_FTBS.exe
	./two_d_advt_burgers_FTBS.exe
two_d_advt_burgers_FTBS.exe: $(OBJECTS_2D_FTBS)
	gfortran $(FFLAGS) $(OBJECTS_2D_FTBS) -o two_d_advt_burgers_FTBS.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
