FLAGS = -fopenmp -O3 -march=native -std=f2003
#FLAGS = -O3 -march=native -std=f2003
#FLAGS = -pg -g -fopenmp -O3 -march=native 

# HDF5 external libraries
# HDF5 = /nfs/cptg/u4/broberts/usr/hdf5-1.8.14
# FORTRANLIB = -I$(HDF5)/include -L$(HDF5)/lib/ $(HDF5)/lib/libhdf5hl_fortran.a $(HDF5)/lib/libhdf5_hl.a $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -lrt -ldl -lm -Wl,-rpath -Wl,$(HDF5)/lib
 
FORTRANLIB = -I${SCINET_HDF5_INC} -L${SCINET_HDF5_LIB} -lhdf5_fortran -lhdf5 -lsz -lz

all: main_rotation

main_rotation: mt95.f90 global_mod.f95 gauss_rand_mod.f95 sort_hilbert_mod.f95 neighbour_list_mod.f95 forces_mod.f95 output_mod.f95 initial_mod.f95 restart_mod.f95 evolve_mod.f95 distribute_mod.f95 mpc_rotation_mod.f95 reaction_simple_mod.f95 main_rotation.f95
	gfortran $(FLAGS) $^ -o $@ $(FORTRANLIB)

clean:
	$(RM) main_rotation *.mod *.dat *.h5 *.pdf *.mp4

