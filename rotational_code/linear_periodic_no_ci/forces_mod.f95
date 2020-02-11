module forces_mod
implicit none

contains

        subroutine LJ_force(motor_number,sphere,NT,Box,flag,type_B,Xsol,Ysol,Zsol,Xdim,Ydim,Zdim, &
                    FpairX,FpairY,FpairZ,Potential,RPab,RPba,cutoff2,Ea,Eb,R6,R12,r12R12,r6R6,random, &
                    neighbour_num,neighbour)
use mt95
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer :: motor_number
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
real(kind=dp) , dimension(:) , intent(in) :: Box                ! Length of simulation box
integer(kind=sp) , dimension(:), intent(inout) :: flag          ! Identity of solvent
integer(kind=sp) , dimension(:), intent(inout) :: type_B          ! Identity of solvent
real(kind=dp) , dimension(:) , intent(in) :: Xsol               ! Force on solvent in x direction
real(kind=dp) , dimension(:) , intent(in) :: Ysol               ! Force on solvent in y direction
real(kind=dp) , dimension(:) , intent(in) :: Zsol               ! Force on solvent in z direction
real(kind=dp) , intent(in) :: Xdim             ! Force on solvent in x direction
real(kind=dp) , intent(in) :: Ydim             ! Force on solvent in y direction
real(kind=dp) , intent(in) :: Zdim             ! Force on solvent in z direction
real(kind=dp) , dimension(:) , intent(inout) :: FpairX      ! LJ force in x direction
real(kind=dp) , dimension(:) , intent(inout) :: FpairY      ! LJ force in y direction
real(kind=dp) , dimension(:) , intent(inout) :: FpairZ      ! LJ force in z direction
real(kind=dp) , dimension(:) , intent(inout) :: Potential   ! Potential of each particle
integer , intent(in) :: sphere                                  ! Incoming sphere identity
real(kind=dp) , intent(in) :: RPab                              ! Probability of reaction from F to G
real(kind=dp) , intent(in) :: RPba                              ! Probability of reaction from G to F
real(kind=dp) , intent(in) :: cutoff2                           ! LJ cutoff for monomer
real(kind=dp) , intent(in) :: Ea                                ! LJ pot. with monomer
real(kind=dp) , intent(in) :: Eb                                ! LJ pot. with monomer
real(kind=dp) , intent(in) :: R6                                ! R^6
real(kind=dp) , intent(in) :: R12                               ! R^12
real(kind=dp) , intent(in) :: r12R12                            ! 12*R^12
real(kind=dp) , intent(in) :: r6R6                              ! 6*R^6
real(kind=dp) , dimension(:) , intent(inout) :: random          ! Vector of random numbers
integer(kind=dp) , intent(in) :: neighbour_num ! Number of solvent particles in list
integer(kind=dp) , dimension(:) , intent(in) :: neighbour   ! Identify solvent number in list

! Definitions for internal variables
integer(kind=dp) :: i,j                                         ! Loop variable
real(kind=dp) :: dRsq                                           ! Solvent-monomer distance, squared
real(kind=dp) , dimension(3) :: dR                              ! Solvent-monomer distance
real(kind=dp) :: IdRsq                                          ! Inverse of the square of dR
real(kind=dp) :: IdR6                                           ! Inverse of (dR)**6
real(kind=dp) :: IdR8                                           ! Inverse of (dR)**8
real(kind=dp) :: IdR12                                          ! Inverse of (dR)**12
real(kind=dp) :: IdR14                                          ! Inverse of (dR)**14
real(kind=dp) :: Fbrack                                         ! LJ force calc., part with radii
real(kind=dp) :: Ubrack                                         ! LJ pot. calc., part with radii

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i,j,dR,dRsq,IdRsq,IdR6,IdR8,IdR12,IdR14,Fbrack,Ubrack) &
!$OMP SHARED(Xdim,Ydim,Zdim,Xsol,Ysol,Zsol,Box,cutoff2,r12R12,r6R6,R12) &
!$OMP SHARED(R6,Ea,Eb,RPab,RPba,random,sphere,flag,type_B,neighbour_num,neighbour) &
!$OMP SHARED(FpairX,FpairY,FpairZ,Potential,motor_number) &
!$OMP SCHEDULE(static)

do i = 1,neighbour_num

 dR(1) = Xsol(neighbour(i)) - Xdim
 dR(2) = Ysol(neighbour(i)) - Ydim
 dR(3) = Zsol(neighbour(i)) - Zdim

 do j = 1,2
  if ( dR(j) > 0.5d0*Box(j) ) then
   dR(j) = dR(j) - Box(j)
  else if ( dR(j) < -0.5d0*Box(j) ) then
   dR(j) = dR(j) + Box(j)
  end if
 end do
 dRsq = dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) 
 
 if ( dRsq <= cutoff2 ) then                                    ! If within cutoff2 region of LJ pot.
  if ( sphere == 1 ) then
   flag(neighbour(i)) = 1                            ! Convert A -> B for catalytic sphere
   type_B(neighbour(i)) = motor_number
  end if
  if ( flag(neighbour(i)) == 0 ) then                  ! For A molecules:
   IdRsq = 1.d0/dRsq                                            ! Inverse of the square of dR
   IdR6 = IdRsq**3                                              ! Inverse of (dR)**6
   IdR8 = IdR6 * IdRsq                                          ! Inverse of (dR)**8
   IdR12 = IdR6 * IdR6                                          ! Inverse of (dR)**12
   IdR14 = IdR12 * IdRsq                                        ! Inverse of (dR)**14
   Fbrack = 4.d0 * ( r12R12*IdR14-r6R6*IdR8 )                   ! LJ force calc., part with radii
   Ubrack = 4.d0*(R12*IdR12-R6*IdR6)+1.d0                       ! LJ pot. calc., part with radii
   FpairX(i) = Ea * Fbrack * dR(1)                     ! LJ force in x direction
   FpairY(i) = Ea * Fbrack * dR(2)                     ! LJ force in y direction
   FpairZ(i) = Ea * Fbrack * dR(3)                     ! LJ force in z direction
   Potential(i) = Potential(i) + Ea * Ubrack  ! LJ Potential

  else if ( flag(neighbour(i)) == 1 ) then             ! If molecule is B, then:
   if ( type_B(neighbour(i)) == motor_number ) then
    IdRsq = 1.d0/dRsq                                            ! Inverse of the square of dR
    IdR6 = IdRsq**3                                              ! Inverse of (dR)**6
    IdR8 = IdR6 * IdRsq                                          ! Inverse of (dR)**8
    IdR12 = IdR6 * IdR6                                          ! Inverse of (dR)**12
    IdR14 = IdR12 * IdRsq                                        ! Inverse of (dR)**14
    Fbrack = 4.d0 * ( r12R12*IdR14-r6R6*IdR8 )                   ! LJ force calc., part with radii
    Ubrack = 4.d0*(R12*IdR12-R6*IdR6)+1.d0                       ! LJ pot. calc., part with radii
    FpairX(i) = Eb * Fbrack * dR(1)                              ! LJ force in x direction
    FpairY(i) = Eb * Fbrack * dR(2)                              ! LJ force in y direction
    FpairZ(i) = Eb * Fbrack * dR(3)                              ! LJ force in z direction
    Potential(i) = Potential(i) + Eb * Ubrack                    ! LJ Potential	
   else
    IdRsq = 1.d0/dRsq                                            ! Inverse of the square of dR
    IdR6 = IdRsq**3                                              ! Inverse of (dR)**6
    IdR8 = IdR6 * IdRsq                                          ! Inverse of (dR)**8
    IdR12 = IdR6 * IdR6                                          ! Inverse of (dR)**12
    IdR14 = IdR12 * IdRsq                                        ! Inverse of (dR)**14
    Fbrack = 4.d0 * ( r12R12*IdR14-r6R6*IdR8 )                   ! LJ force calc., part with radii
    Ubrack = 4.d0*(R12*IdR12-R6*IdR6)+1.d0                       ! LJ pot. calc., part with radii
    FpairX(i) = Ea * Fbrack * dR(1)                              ! LJ force in x direction
    FpairY(i) = Ea * Fbrack * dR(2)                              ! LJ force in y direction
    FpairZ(i) = Ea * Fbrack * dR(3)                              ! LJ force in z direction
    Potential(i) = Potential(i) + Ea * Ubrack                    ! LJ Potential	
   end if
  end if
 end if

end do

!$OMP END PARALLEL DO

end subroutine LJ_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LJ_force_wall(sphere,box,zdim,fwallz,Ew,R3,R9,r9R9,r3R3,Potential,k)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
real(kind=dp) , dimension(:) , intent(in) :: Box                ! Length of simulation box
real(kind=dp) , dimension(:,:) , intent(in) :: Zdim             ! Pos. of monomer in z direction
real(kind=dp) , intent(inout) :: fwallz                         ! Force of monomer with wall, z
integer , intent(in) :: sphere                                  ! Incoming sphere identity
real(kind=dp) , intent(in) :: Ew                                ! LJ pot. of wall with monomer
real(kind=dp) , intent(in) :: R3                                ! R^3
real(kind=dp) , intent(in) :: R9                                ! R^9
real(kind=dp) , intent(in) :: r9R9                              ! 9*R^9
real(kind=dp) , intent(in) :: r3R3                              ! 3*R^3
real(kind=dp) , dimension(:) , intent(inout) :: Potential       ! Potential of each particle

! Definitions for internal variables
integer(kind=dp) :: i                                           ! Loop variable
real(kind=dp) :: dR                                             ! wall-monomer distance
real(kind=dp) :: IdRsq                                          ! Inverse of the square of dRsq
real(kind=dp) :: IdR3                                           ! Inverse of (dR)**6
real(kind=dp) :: IdR9                                           ! Inverse of (dR)**8
real(kind=dp) :: IdR5                                           ! Inverse of (dR)**12
real(kind=dp) :: IdR11                                          ! Inverse of (dR)**14
real(kind=dp) :: Fbrack                                         ! LJ force calc., part with radii
integer :: k                                                    ! Identity of dimer motor

! Calculate force on upper wall of z-axis from monomer
dR = box(3) - zdim(k,sphere)

IdRsq = 1. / (dR**2)
IdR3 = IdRsq ** (3./2.)
IdR9 = IdR3 ** 3
IdR5 = IdR3 * IdRsq
IdR11 = IdR9 * IdRsq

Fbrack = 4*Ew*(r9R9*IdR11-r3R3*IdR5)
fwallz = Fbrack * dR
Potential(sphere) = 4*Ew*(R9*IdR9-R3*IdR3)

! Calculate force between lower wall of z-axis with monomer
dR = 0. - zdim(k,sphere)

IdRsq = 1. / (dR**2)
IdR3 = IdRsq ** (3./2.)
IdR9 = IdR3 ** 3
IdR5 = IdR3 * IdRsq
IdR11 = IdR9 * IdRsq

Fbrack = 4*Ew*(r9R9*IdR11-r3R3*IdR5)
fwallz = fwallz + Fbrack * dR
Potential(sphere) = Potential(sphere) + 4*Ew*(R9*IdR9-R3*IdR3)

end subroutine LJ_force_wall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine harmonic_wall(zdim,zdimeq,k_harm,f_harm,pot_harm)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
real(kind=dp) , intent(in) :: zdim             ! Pos. of monomer in z direction
real(kind=dp) , intent(in) :: zdimeq                            ! Equil. pos. of mon., z dir.
real(kind=dp) , intent(in) :: k_harm                            ! Force constant
real(kind=dp) , intent(out) :: f_harm                           ! Force on monomer
real(kind=dp) , intent(out) :: pot_harm                         ! Potential of each monomer

! Internal variables
real(kind=dp) :: zdiff                                          ! Diff. from equil. pos., z

! Difference from equilibrium position
zdiff = zdim - zdimeq

! Calculate force on monomer from harmonic potential
f_harm = -1. * k_harm * zdiff 
pot_harm = 0.5 * k_harm * zdiff**2

end subroutine harmonic_wall

end module
