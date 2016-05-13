OBJECTS_1D_FTCS = boundary_conditions.o scheme.o initialize.o one_d_diffusion_FTCS.o
OBJECTS_1D_BTCS = boundary_conditions.o scheme.o initialize.o one_d_diffusion_BTCS.o

FFLAGS = -g

.PHONY:   clean clean_doc test

one_d_diffusion_FTCS: one_d_diffusion_FTCS.exe
	./one_d_diffusion_FTCS.exe
one_d_diffusion_FTCS.exe: $(OBJECTS_1D_FTBS)
	gfortran $(FFLAGS) $(OBJECTS_1D_FTBS) -o one_d_diffusion_FTCS.exe
one_d_diffusion_BTCS: one_d_diffusion_BTCS.exe
	./one_d_diffusion_BTCS.exe
one_d_diffusion_BTCS.exe: $(OBJECTS_1D_BTCS)
	gfortran $(FFLAGS) $(OBJECTS_1D_BTCS) -o one_d_diffusion_BTCS.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
