module reaction_simple_mod
implicit none

contains

        subroutine reaction_simple(MPCt,numcell,k1,flag,type_B,tBCellN,sBCellN)
use mt95
implicit none

!This subroutine calculates the probability for the reactions in the
!multiparticle collision cells and performs the reaction.

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=dp) , intent(in) :: numcell                        ! Number of subcells
real(kind=dp) , intent(in) :: MPCt                              ! MPC collision time
real(kind=dp) , intent(in) :: k1                                ! Reaction rate for B -> A
integer(kind=sp) , dimension(:), intent(inout) :: flag          ! Identity of solvent
integer(kind=sp) , dimension(:), intent(inout) :: type_B          ! Identity of solvent
integer(kind=dp) , dimension(:), intent(in) :: tBCellN          ! Total number of B solvent in cell N
integer(kind=dp) , dimension(:,:), intent(in) :: sBCellN        ! Label I in Cell A with i

! Definitions for internal variables
integer(kind=dp) :: i                                           ! Loop variable
real(kind=dp) :: x                                              ! Random variable
real(kind=dp) :: a00                                            ! a00 = a1+a2+a3+a4+a5+a6
real(kind=dp) :: a1                                             ! Propensity of reaction 1
real(kind=dp) :: prob1                                          ! Probability of reaction 1

!***********************************************************************************!
!******************************Reactions********************************************!
!***********************************************************************************!

do i = 1,numcell
 
!!!!!!Probability calculations
 !Probability factor a_{}
 if ( tBCellN(i) == 0 ) then
  a1 = 0.
 else if ( tBCellN(i) > 0 ) then
  a1 = k1 * tBCellN(i)
 end if

 a00 = a1

 !Probability for back reaction2
 prob1 = 1.d0 - dexp(-a00*MPCt)

 call genrand_real1(x)
 if ( 0. < x .and. x <= prob1 ) then
  flag(sBCellN(i,1)) = 0
  type_B(sBCellN(i,1)) = 0
 end if

end do

end subroutine reaction_simple

end module reaction_simple_mod
