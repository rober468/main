all: Main_SelkovDimer

Main_SelkovDimer: mt95.f90 Global_SelkovDimer.f90 Affected_Cell.o Distribute_SelkovDimer.o Initial_SelkovDimer.o Restart.o Motion_SelkovDimer.o ForcesMD_SelkovDimer.o RMPC_SelkovDimer.f90 Output_SelkovDimer.f90 Gaussian_Random.o Sort.o Output_Concentration_Catalytic.o Main_SelkovDimer.o
	gfortran -g -pg -fopenmp -march=native -o $@ $^

Main_SelkovDimer.o: mt95.f90 Global_SelkovDimer.f90 Main_SelkovDimer.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Initial_SelkovDimer.o: mt95.f90 Global_SelkovDimer.f90 Initial_SelkovDimer.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Affected_Cell.o: mt95.f90 Global_SelkovDimer.f90 Affected_Cell.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Distribute_SelkovDimer.o: mt95.f90 Global_SelkovDimer.f90 Distribute_SelkovDimer.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Restart.o: mt95.f90 Global_SelkovDimer.f90 Restart.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Motion_SelkovDimer.o: mt95.f90 Global_SelkovDimer.f90 Motion_SelkovDimer.f90
	gfortran -fopenmp -O3 -march=native $^ -c

ForcesMD_SelkovDimer.o: mt95.f90 Global_SelkovDimer.f90 ForcesMD_SelkovDimer.f90
	gfortran -fopenmp -O3 -march=native $^ -c

RMPC_SelkovDimer.o: mt95.f90 Global_SelkovDimer.f90 RMPC_SelkovDimer.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Output_SelkovDimer.o: mt95.f90 Global_SelkovDimer.f90 Output_SelkovDimer.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Gaussian_Random.o: mt95.f90 Global_SelkovDimer.f90 Gaussian_Random.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Sort.o: mt95.f90 Global_SelkovDimer.f90 Sort.f90
	gfortran -fopenmp -O3 -march=native $^ -c

Output_Concentration_Catalytic.o: mt95.f90 Global_SelkovDimer.f90 Output_Concentration_Catalytic.f90
	gfortran -fopenmp -O3 -march=native $^ -c

clean:
	$(RM) Main_SelkovDimer *.o *.mod
