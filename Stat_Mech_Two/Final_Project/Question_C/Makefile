CC = gfortran
FLAGS = -fopenmp -O3 -march=native

all: Main_Reactive_Flux

Main_Reactive_Flux: mt95.f90 Global.f95 Initialize.o Output_Time_Dep_Rate.o Velocity_Verlet_Evolve.o Gaussian_Random.o Main_Reactive_Flux.o
	$(CC) $(FLAGS) -o $@ $^

Main_Reactive_Flux.o: Global.f95 Main_Reactive_Flux.f95
	$(CC) $(FLAGS) $^ -c

Initialize.o: Global.f95 Initialize.f95
	$(CC) $(FLAGS) $^ -c

Output_Time_Dep_Rate.o: Global.f95 Output_Time_Dep_Rate.f95
	$(CC) $(FLAGS) $^ -c

Velocity_Verlet_Evolve.o: Global.f95 Velocity_Verlet_Evolve.f95
	$(CC) $(FLAGS) $^ -c

Gaussian_Random.o: mt95.f90 Gaussian_Random.f95
	$(CC) $(FLAGS) $^ -c

clean:
	rm -f Main_Reactive_Flux *.mod *.o momenta*.dat position*.dat
