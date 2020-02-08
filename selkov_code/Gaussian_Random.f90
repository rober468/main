SUBROUTINE Gaussian_Random (mean,var,Vrandom)
USE mt95
IMPLICIT NONE

!This subroutine chooses a normal random variable (mean=x, st=0) from a randomly distributed variable on [0,1) using Box Muller method.

!Define kind value for precision
INTEGER, PARAMETER::SGL=4
INTEGER, PARAMETER::DBL=8

!Declaring intent of variables between program and subroutine
REAL(KIND=DBL),INTENT(IN)::mean			!Mean of Gaussian distribution
REAL(KIND=DBL),INTENT(IN)::var			!Variance of Gaussian distribution
REAL(KIND=DBL),INTENT(OUT)::Vrandom		!Random velocity from Gaussian distr.

!Other variables
REAL(KIND=DBL)::xx				!Variable for random number
REAL(KIND=DBL)::yy				!Variable for random number
REAL(KIND=DBL)::theta				!Value for polar coord. rand. variable
REAL(KIND=DBL)::R				!Value for polar coord. rand. variable
REAL(KIND=DBL)::z				!Random normal number
REAL(KIND=DBL),PARAMETER::pi=4*DATAN(1.d0)	!pi=3.14159265358979323846264...

CALL genrand_real1(xx)				!Randomize x for xE[0,1)
CALL genrand_real1(yy)				!Randomize y for yE[0,1)

theta=2*pi*xx					!Independent random to theta
R=SQRT(-2*DLOG(yy))				!Independent random to R

z=R*DCOS(theta)					!Random normal number

Vrandom=SQRT(var)*z+mean				!Random velocity from Gaussian distr.

END SUBROUTINE Gaussian_Random
