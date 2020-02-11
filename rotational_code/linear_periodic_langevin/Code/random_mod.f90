module random_mod
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gauss_rand(mean,var,Vrandom)
use mt95
implicit none

! This subroutine chooses a normal random variable (mean=x, st=0) from a randomly distributed variable on [0,1) using Box Muller method.

! precision
integer, parameter :: sp = selected_real_kind(6,37)     ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)   ! Double precision

! declaring intent of variables between program and subroutine
real(kind=dp) , intent(in) :: mean                      ! Mean of Gaussian distribution
real(kind=dp) , intent(in) ::var                        ! Variance of Gaussian distribution
real(kind=dp) , intent(out) :: Vrandom                  ! Random velocity from Gaussian distr.

! other variables
real(kind=dp) ::xx                                      ! Variable for random number
real(kind=dp) ::yy                                      ! Variable for random number
real(kind=dp) ::theta                                   ! Value for polar coord. rand. variable
real(kind=dp) ::R                                       ! Value for polar coord. rand. variable
real(kind=dp) ::z                                       ! Random normal number
real(kind=dp) , parameter :: pi=4*DATAN(1.d0)           ! pi=3.14159265358979323846264...

call genrand_real1(xx)                                  ! Randomize x for xE[0,1)
call genrand_real1(yy)                                  ! Randomize y for yE[0,1)

theta=2*pi*xx                                           ! Independent random to theta
R=SQRT(-2*DLOG(yy))                                     ! Independent random to R

z=R*DCOS(theta)                                         ! Random normal number

Vrandom=SQRT(var)*z+mean                                ! Random velocity from Gaussian distr.

end subroutine gauss_rand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine circumference_2d_drand(rand)
use mt95
implicit none

! This subroutine calculates a double precision random angle on the circumference of a circle
! from [0,2pi).

! precision
integer, parameter :: dp = selected_real_kind(15,307)   ! Double precision

! external variables
real(kind=dp) :: rand

! internal variables
real(kind=dp) , parameter :: pi = 4.*atan(1.d0)

call genrand_real1(rand)

rand = 2.*pi*rand

end subroutine circumference_2d_drand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module random_mod
