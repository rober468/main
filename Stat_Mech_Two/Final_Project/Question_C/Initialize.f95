subroutine Initialize
use Global
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
! This subroutine initializes the oscillator positions and all momenta at t=0.		!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate mean and variances
meanoscpos = 0.
varoscpos = 1. / beta 
meanoscmom = 0.
varoscmom = 1. / beta
meanrxnmom = 0.
varrxnmom = 1. / beta


! Pick initial values
do i = 1,runs,1
	! Oscillator positions and momenta
	do j = 1,int(N),1
		CALL Gaussian_Random(meanoscpos,varoscpos,R_i(i,j))
		CALL Gaussian_Random(meanoscmom,varoscmom,P_i(i,j))
	end do
	! Rxn coordinate momentum
	CALL Gaussian_Random(meanrxnmom,varrxnmom,P_0(i))
	P_00(i) = P_0(i)
end do

! Calculate initial force
do i = 1,runs,1
	! Rxn coordinate force
	F_0(i) = 0.
	do j = 1,int(N),1
		F_0(i) = F_0(i) + sqrt( xi/N*(1.-exp(-omega_max)) ) * R_i(i,j)
	end do
	do j = 1,int(N),1
		F_i(i,j) = -R_i(i,j) 
	end do
end do

open(unit=10,file="position_osc.dat",action='write',status='new')
open(unit=11,file="momenta_osc.dat",action='write',status='new')
open(unit=12,file="momenta_rxn.dat",action='write',status='new')
do i = 1,runs,1
        write(10,*) (R_i(i,j) , j = 1,int(N),1)
        write(11,*) (P_i(i,j) , j = 1,int(N),1)
        write(12,*) P_0(i)
end do
close(10)
close(11)
close(12)

end subroutine Initialize
