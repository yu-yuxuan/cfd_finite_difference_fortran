OBJECTS_2D_FTCS = boundary_conditions.o scheme.o initialize.o two_d_diffusion_FTCS.o


FFLAGS = -g

.PHONY:   clean clean_doc test

two_d_diffusion_FTCS: two_d_diffusion_FTCS.exe
	./two_d_diffusion_FTCS.exe
two_d_diffusion_FTCS.exe: $(OBJECTS_2D_FTCS)
	gfortran $(FFLAGS) $(OBJECTS_2D_FTCS) -o two_d_diffusion_FTCS.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
