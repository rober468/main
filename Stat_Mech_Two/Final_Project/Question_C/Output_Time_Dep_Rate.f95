subroutine Output_Time_Dep_Rate
use Global
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									!
! This program calculates the time-dependent rate constant.		!
!									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate average quantity of A at equilibrium
R_A = -sqrt( b/a + xi/a*(1.-exp(-omega_max)) )
d2UdR02 = 3.*a*R_A**2 - (b+xi*(1.-exp(-omega_max)))
W_R_A = a/4.*R_A**4 - b/2.*R_A**2 - xi/2.*(1.-exp(-omega_max))*R_A**2
M = sqrt(2.*pi/(beta*d2UdR02)) * exp(-beta*W_R_A)

! Calculate average over all runs for at each time.
av = 0.
do j = 1,runs,1
	if ( R_0(j) >= 0. ) then
		theta_0(j) = 1.
	else if ( R_0(j) < 0. ) then
		theta_0(j) = 0.
	end if
	av = av + theta_0(j) * P_00(j) / M
end do
k_f(cstep) = av / runs

! Calculate standard error from estimator
st_err(cstep) = 0.
do i = 1,runs,1
	st_err(cstep) = st_err(cstep) + ( ( theta_0(j) * P_00(j) / K ) - k_f(cstep) )**2
end do
st_err(cstep) = sqrt( (st_err(cstep)/float(runs)) ) / sqrt(float(runs)-1.)

! Output rate constant data
open(unit=10,file="Time_Dep_Rate.dat",action="write",status="replace",access="append")
write(10,*) k_f(cstep) , st_err(cstep)
close(10)

! Output trajectory of x
open(unit=10,file="trajectory_R_0.dat",action="write",status="replace",access="append")
write(10,*) ( R_0(i) , i=1,10 )
close(10)

end subroutine Output_Time_Dep_Rate
