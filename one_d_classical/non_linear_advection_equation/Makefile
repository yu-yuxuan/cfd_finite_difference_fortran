OBJECTS_1D_FTBS = boundary_conditions.o scheme.o initialize.o one_d_advt_burgers_FTBS.o

FFLAGS = -g

.PHONY:   clean clean_doc test

one_d_advt_burgers_FTBS: one_d_advt_burgers_FTBS.exe
	./one_d_advt_burgers_FTBS.exe
one_d_advt_burgers_FTBS.exe: $(OBJECTS_1D_FTBS)
	gfortran $(FFLAGS) $(OBJECTS_1D_FTBS) -o one_d_advt_burgers_FTBS.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
