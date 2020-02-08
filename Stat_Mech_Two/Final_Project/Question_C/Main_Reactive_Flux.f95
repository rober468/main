program Main_Reactive_Flux
use mt95
use Global
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!										        !
! Main program to calculate the rate constant via the reactive flux method of a quartic !
! bistable system coupled to a harmonic bath with Hamiltonian:                          !
!											!
!   H = P_0**2/2m + V_0 + sum_{i=1}^N [P_i**2/(2M) + 0.5Momega_i**2R_i**2 - R_0c_iR_i]	!
!											!
! where V_0 = aR_0**4 / 4 - bR_0**2 / 2 . Numerical calculation performed by first	!
! using Monte-Carlo sampling (via Gaussian distributions) at t=0 for initial states,	!
! then evolving the system in time with MD using the velocity Verlet algorithm.		!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(R_i(1:runs,1:int(N)))
allocate(P_i(1:runs,1:int(N)))
allocate(F_i(1:runs,1:int(N)))
allocate(R_0(1:runs))
allocate(P_0(1:runs))
allocate(P_00(1:runs))
allocate(F_0(1:runs))
allocate(k_f(1:time+1))
allocate(st_err(1:time+1))
allocate(theta_0(1:runs))

! Initialize PRNG
seed = 24106587
CALL genrand_init( put=seed )

! Initialize using Monte-Carlo sampling
CALL Initialize

! Time evolution of quantities by velocity Verlet algorithm for many runs
cstep = 1.
do while (int(cstep) < time)
        CALL Velocity_Verlet_Evolve
	CALL Output_Time_Dep_Rate
        cstep = cstep + 1.
        if (mod(cstep,200) == 0) then
                write(*,*) cstep
        end if
end do

deallocate(R_i)
deallocate(P_i)
deallocate(F_i)
deallocate(R_0)
deallocate(P_0)
deallocate(P_00)
deallocate(F_0)
deallocate(k_f)
deallocate(st_err)
deallocate(theta_0)

end program
