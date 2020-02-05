module reaction_mod
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reaction_simple(MPCt,numcell,k1,flag,tBCellN,sBCellN,mt,mti)
use mtmod
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
integer(kind=dp) , dimension(:), intent(in) :: tBCellN          ! Total number of B solvent in cell N
integer(kind=dp) , dimension(:,:), intent(in) :: sBCellN        ! Label I in Cell A with i
integer , dimension(:,:) , intent(inout) :: mt
integer , dimension(:) , intent(inout) :: mti

! Definitions for internal variables
integer(kind=dp) :: i                                           ! Loop variable
real(kind=dp) :: random                                         ! Random variable
real(kind=dp) :: a00                                            ! a00 = a1+a2+a3+a4+a5+a6
real(kind=dp) :: a1                                             ! Propensity of reaction 1
real(kind=dp) :: prob1                                          ! Probability of reaction 1
integer :: thread_num,omp_get_thread_num

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numcell,tbcelln,sbcelln,k1,MPCt,flag,mt,mti) &
!$OMP PRIVATE(i,a1,a00,prob1,thread_num,random)

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

 thread_num = omp_get_thread_num()+1
 random = grnd(mt(:,thread_num),mti(thread_num))
 if ( 0. < random .and. random <= prob1 ) then
  flag(sBCellN(i,1)) = 0
 end if

end do

!$OMP END PARALLEL DO

end subroutine reaction_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reaction_depletion(MPCt,numcell,k1,k2,flag,tACellN,sACellN,tBCellN,sBCellN,mt,mti)
use mtmod
implicit none

!This subroutine calculates the probability for the reactions in the
!multiparticle collision cells and performs the reaction.
!
! Reaction: B -> S irreversible, A -> S irreversible
!

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=dp) , intent(in) :: numcell                        ! Number of subcells
real(kind=dp) , intent(in) :: MPCt                              ! MPC collision time
real(kind=dp) , intent(in) :: k1                                ! Reaction rate for A -> S
real(kind=dp) , intent(in) :: k2                                ! Reaction rate for B -> S
integer(kind=sp) , dimension(:), intent(inout) :: flag          ! Identity of solvent
integer(kind=dp) , dimension(:), intent(in) :: tACellN          ! Total number of A solvent in cell N
integer(kind=dp) , dimension(:,:), intent(in) :: sACellN        ! Label I in Cell A with i
integer(kind=dp) , dimension(:), intent(in) :: tBCellN          ! Total number of B solvent in cell N
integer(kind=dp) , dimension(:,:), intent(in) :: sBCellN        ! Label I in Cell B with i
integer , dimension(:,:) , intent(inout) :: mt
integer , dimension(:) , intent(inout) :: mti

! Definitions for internal variables
integer(kind=dp) :: i                                           ! Loop variable
real(kind=dp) :: random                                              ! Random variable
real(kind=dp) :: a00                                            ! a00 = a1+a2
real(kind=dp) :: a1                                             ! Propensity of reaction 1
real(kind=dp) :: a2                                             ! Propensity of reaction 1
real(kind=dp) :: prob_rxn                                       ! Probability of reaction
integer :: thread_num,omp_get_thread_num

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numcell,tacelln,sacelln,tbcelln,sbcelln,k1,k2,MPCt,flag,mt,mti) &
!$OMP PRIVATE(i,a1,a2,a00,prob_rxn,thread_num,random)

do i = 1,numcell

!!!!!!Probability calculations
 !Probability factor a_{}
 if ( tACellN(i) == 0 ) then
  a1 = 0.
 else if ( tACellN(i) > 0 ) then
  a1 = k1 * tACellN(i)
 end if
 if ( tBCellN(i) == 0 ) then
  a2 = 0.
 else if ( tBCellN(i) > 0 ) then
  a2 = k2 * tBCellN(i)
 end if

 a00 = a1+a2

 !Probability for a reaction
 prob_rxn = 1.d0 - dexp(-a00*MPCt)

 thread_num = omp_get_thread_num()+1
 random = grnd(mt(:,thread_num),mti(thread_num))
 if ( 0. < random .and. random <= prob_rxn ) then

  random = grnd(mt(:,thread_num),mti(thread_num))
  random = a00*random
  if ( 0. < random .and. random < a1 ) then
   flag(sACellN(i,1)) = 2
  else if ( a1 <= random .and. random < a00 ) then
   flag(sBCellN(i,1)) = 2
  end if

 end if

end do

!$OMP END PARALLEL DO

end subroutine reaction_depletion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module reaction_mod
