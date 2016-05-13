OBJECTS_2D_CF = boundary_conditions.o scheme.o initialize.o channel_flow.o


FFLAGS = -g

.PHONY:   clean clean_doc test
channel_flow: channel_flow.exe
	./channel_flow.exe
channel_flow.exe: $(OBJECTS_2D_CF)
	gfortran $(FFLAGS) $(OBJECTS_2D_CF) -o channel_flow.exe
%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
