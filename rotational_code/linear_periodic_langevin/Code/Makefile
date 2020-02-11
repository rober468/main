FLAGS = -fopenmp -O3 -march=native #-std=f2003
#FLAGS = -O3 -march=native -std=f2003
#FLAGS = -pg -g -fopenmp -O3 -march=native 

# HDF5 external libraries
# HDF5 = /nfs/cptg/u4/broberts/usr/hdf5-1.8.14
# FORTRANLIB = -I$(HDF5)/include -L$(HDF5)/lib/ $(HDF5)/lib/libhdf5hl_fortran.a $(HDF5)/lib/libhdf5_hl.a $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -lrt -ldl -lm -Wl,-rpath -Wl,$(HDF5)/lib
 
FORTRANLIB = -I${SCINET_HDF5_INC} -L${SCINET_HDF5_LIB} -lhdf5_fortran -lhdf5 -L${SCINET_SZIP_LIB}-lsz -lz

#FORTRANLIB = -I/scinet/gpc/Libraries/HDF5-1.8.7/v18-serial-gcc/include -L/scinet/gpc/Libraries/HDF5-1.8.7/v18-serial-gcc/lib /scinet/gpc/Libraries/HDF5-1.8.7/v18-serial-gcc/lib/libhdf5hl_fortran.a /scinet/gpc/Libraries/HDF5-1.8.7/v18-serial-gcc/lib/libhdf5_hl.a /scinet/gpc/Libraries/HDF5-1.8.7/v18-serial-gcc/lib/libhdf5_fortran.a /scinet/gpc/Libraries/HDF5-1.8.7/v18-serial-gcc/lib/libhdf5.a -L/scinet/gpc/Libraries/HDF5-1.8.7/zlib-gcc/lib -L/scinet/gpc/Libraries/HDF5-1.8.7/szip-gcc/lib -lsz -lz -lm -lgpfs -Wl,-rpath -Wl,/scinet/gpc/Libraries/HDF5-1.8.7/v18-serial-gcc/lib

all: main_langevin_rotors

main_langevin_rotors: mt95.f90 dqng_mod.f90 global_mod.f90 random_mod.f90 forces_mod.f90 derived_quantities_mod.f90 output_mod.f90 initial_mod.f90 evolve_mod.f90 main_langevin_rotors.f90
	gfortran $(FLAGS) $^ -o $@ $(FORTRANLIB)

clean:
	$(RM) main_langevin_rotors *.mod *.dat *.h5 *.pdf *.mp4

