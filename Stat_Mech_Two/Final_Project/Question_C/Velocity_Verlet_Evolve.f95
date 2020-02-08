subroutine Velocity_Verlet_Evolve
use Global
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
! Evolves the system according to the velocity Verlet algorithm.			!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j = 1,runs,1

        ! Evolve reaction coordinate and oscillator up to force recalculation
        P_0(j) = P_0(j) + (tstep/2.)*F_0(j)
        R_0(j) = R_0(j) + P_0(j)*tstep
        do k = 1,int(N),1
	        P_i(j,k) = P_i(j,k) + (tstep/2.)*F_i(j,k)
	        R_i(j,k) = R_i(j,k) + P_i(j,k)*tstep
        end do

        ! Recalculate forces
        F_0(j) = -a*R_0(j)**3 + b*R_0(j)
        do k = 1,int(N),1
	        F_0(j) = F_0(j) + sqrt( xi/N*(1.-exp(-omega_max)) ) * R_i(j,k)
        end do
        do k = 1,int(N),1
                F_i(j,k) = -R_i(j,k) + sqrt( xi/N*(1.-exp(-omega_max)) ) * R_0(j)
        end do

        ! Final momentum displacement
        P_0(j) = P_0(j) + (tstep/2.)*F_0(j)
        do k = 1,int(N),1
	        P_i(j,k) = P_i(j,k) + (tstep/2.)*F_i(j,k)
        end do

end do

end subroutine
