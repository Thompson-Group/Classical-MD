
%.o : %.f90
	pgf90 -O3 -c $<

OBJECTS = kinds.o common_variables.o forces.o calculate_bond.o calc_angle.o calc_vdw.o coulomb.o\
      coulomb_sf.o coulomb_dsf.o fcheck.o dists.o integrator.o read_data.o thermostat.o read_input.o\
      initvel.o dump.o thermo_dump.o main.o 


main: $(OBJECTS)
	pgf90 -O3 $(OBJECTS) -o main -lm 

