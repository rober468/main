MODULE Global
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                             !
! This module contains all the shared constants used through the main program and subroutines !
!                                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Precision of variables
integer, parameter :: sp = selected_real_kind(6,37)	! Single precision
integer, parameter :: dp = selected_real_kind(15,307)	! Double precision

! Random number seed
integer(kind=selected_int_kind(9)) :: seed

! Simulation Parameters
real(kind=dp), parameter :: a = 3.0			! Quartic potential parameter (quartic)
real(kind=dp), parameter :: b = 2.0			! Quartic potential parameter (quadratic)
real(kind=dp), parameter :: omega_max = 3.0		! Frequency cutoff
real(kind=dp), parameter :: xi = 1.0			! Kondo coupling parameter
real(kind=dp), parameter :: beta = 6.			! Inverse temperature
real(kind=dp), parameter :: N = 1.			! Number of oscillators

! Other parameters
integer(kind=dp), parameter :: runs = 10000		! Number of runs for average
integer(kind=dp), parameter :: time = 5000		! Time for MD simulations
real(kind=dp), parameter :: tstep = 0.005		! MD time step

! Other Constants
real(kind=dp), parameter :: pi = 3.14159265358979323826264

! Global Variables
integer(kind=dp) :: i,j,k,l
integer(kind=dp) :: cstep
real(kind=dp) :: meanoscpos
real(kind=dp) :: varoscpos
real(kind=dp) :: meanoscmom
real(kind=dp) :: varoscmom
real(kind=dp) :: meanrxnmom
real(kind=dp) :: varrxnmom
real(kind=dp), dimension(:,:), allocatable :: R_i
real(kind=dp), dimension(:,:), allocatable :: P_i
real(kind=dp), dimension(:,:), allocatable :: F_i
real(kind=dp), dimension(:), allocatable :: R_0
real(kind=dp), dimension(:), allocatable :: P_0
real(kind=dp), dimension(:), allocatable :: P_00
real(kind=dp), dimension(:), allocatable :: F_0
real(kind=dp), dimension(:), allocatable :: k_f
real(kind=dp), dimension(:), allocatable :: theta_0
real(kind=dp) :: M
real(kind=dp) :: av
real(kind=dp) :: R_A
real(kind=dp) :: d2UdR02
real(kind=dp) :: W_R_A
real(kind=dp), dimension(:), allocatable :: st_err

END MODULE Global
