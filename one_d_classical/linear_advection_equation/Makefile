OBJECTS_1D_FTBS = boundary_conditions.o scheme.o initialize.o one_d_linear_advection_FTBS.o
OBJECTS_1D_FTCS = boundary_conditions.o scheme.o initialize.o one_d_linear_advection_FTCS.o
OBJECTS_1D_BTCS = boundary_conditions.o scheme.o initialize.o one_d_linear_advection_BTCS.o
OBJECTS_1D_laxf = boundary_conditions.o scheme.o initialize.o one_d_linear_advt_lax_fr.o
OBJECTS_1D_laxw = boundary_conditions.o scheme.o initialize.o one_d_linear_advt_lax_wen.o
OBJECTS_1D_leap = boundary_conditions.o scheme.o initialize.o one_d_linear_advt_leapfrg_ftbs.o


FFLAGS = -g

.PHONY:   clean clean_doc test

one_d_linear_advection_FTBS: one_d_linear_advection_FTBS.exe
	./one_d_linear_advection_FTBS.exe
one_d_linear_advection_FTBS.exe: $(OBJECTS_1D_FTBS)
	gfortran $(FFLAGS) $(OBJECTS_1D_FTBS) -o one_d_linear_advection_FTBS.exe
one_d_linear_advection_FTCS: one_d_linear_advection_FTCS.exe
	./one_d_linear_advection_FTCS.exe
one_d_linear_advection_FTCS.exe: $(OBJECTS_1D_FTCS)
	gfortran $(FFLAGS) $(OBJECTS_1D_FTCS) -o one_d_linear_advection_FTCS.exe
one_d_linear_advection_BTCS: one_d_linear_advection_BTCS.exe
	./one_d_linear_advection_BTCS.exe
one_d_linear_advection_BTCS.exe: $(OBJECTS_1D_BTCS)
	gfortran $(FFLAGS) $(OBJECTS_1D_BTCS) -o one_d_linear_advection_BTCS.exe
one_d_linear_advt_lax_fr: one_d_linear_advt_lax_fr.exe
	./one_d_linear_advt_lax_fr.exe
one_d_linear_advt_lax_fr.exe: $(OBJECTS_1D_laxf)
	gfortran $(FFLAGS) $(OBJECTS_1D_laxf) -o one_d_linear_advt_lax_fr.exe
one_d_linear_advt_lax_wen: one_d_linear_advt_lax_wen.exe
	./one_d_linear_advt_lax_wen.exe
one_d_linear_advt_lax_wen.exe: $(OBJECTS_1D_laxw)
	gfortran $(FFLAGS) $(OBJECTS_1D_laxw) -o one_d_linear_advt_lax_wen.exe
one_d_linear_advt_leapfrg_ftbs: one_d_linear_advt_leapfrg_ftbs.exe
	./one_d_linear_advt_leapfrg_ftbs.exe
one_d_linear_advt_leapfrg_ftbs.exe: $(OBJECTS_1D_leap)
	gfortran $(FFLAGS) $(OBJECTS_1D_leap) -o one_d_linear_advt_leapfrg_ftbs.exe




test: test.exe
	./test.exe

test.exe: test_mode.o test.o
	gfortran $(FFLAGS) test_mode.o test.o  -o test.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $<
clean_all: clean clean_doc

clean:
	rm -f *.o *.exe *.mod
clean_doc:
	rm -f *.txt *.png
