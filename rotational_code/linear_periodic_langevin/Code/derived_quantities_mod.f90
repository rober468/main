module derived_quantities_mod
implicit none



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diffusion_coefficient_MPC(kBT,tau,alpha,dens,a,m,D)
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! external variables
real(kind=dp), intent(in) :: kBT        ! temperature
real(kind=dp), intent(in) :: tau        ! MPC time
real(kind=dp), intent(in) :: dens       ! number density of solvent
real(kind=dp), intent(in) :: alpha      ! rotation angle
real(kind=dp), intent(in) :: a          ! length of MPC cell
real(kind=dp), intent(in) :: m          ! mass of solvent
real(kind=dp), intent(out) :: D         ! diffusion coefficient

! internal variables
real(kind=dp) :: vol                    ! volume of MPC cell
real(kind=dp) :: avnum                  ! average number of particles per cell

! derived quantities
vol = a*a*a
avnum = dens*vol

! diffusion coefficient
D = (kBT*tau/(2.*m)) * (3*avnum/((avnum-1.+dexp(-avnum))*(1.-dcos(alpha)))-1.)

end subroutine diffusion_coefficient_MPC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine derjaguin_length_N(Rn,Ean,Ebn,derj_len)
use dqng_mod
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! external variables
real(kind=dp), intent(in) :: Rn          ! radius of sphere
real(kind=dp), intent(in) :: Ean         ! interaction energy of A with sphere
real(kind=dp), intent(in) :: Ebn         ! interaction energy of B with sphere
real(kind=dp), intent(out) :: derj_len   ! Derjaguin length

! internal variables
real(kind=dp) :: sigma_Rn
real(kind=dp) :: abserr
integer :: nevals,ier

! derived quantity
sigma_Rn = Rn*2.**(1./6.)

call dqng(derj_integrand_N,0.d0,sigma_Rn,1.d-10,1.d-8,derj_len,abserr,nevals,ier)

end subroutine derjaguin_length_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kappa_decay(k1,k_1,D,kappa)
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! external variables
real(kind=dp), intent(in) :: k1
real(kind=dp), intent(in) :: k_1
real(kind=dp), intent(in) :: D
real(kind=dp), intent(out) :: kappa

! internal variables

kappa = sqrt( (k1+k_1) / D )

end subroutine kappa_decay

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine conc_exp_prefactor(p_plus,p_minus,Rc,kBT,m,k_D,D,kappa,c_A_SS,c_B_SS,Z)
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! external variables
real(kind=dp), intent(in) :: p_plus
real(kind=dp), intent(in) :: p_minus
real(kind=dp), intent(in) :: Rc
real(kind=dp), intent(in) :: kBT
real(kind=dp), intent(in) :: m
real(kind=dp), intent(in) :: k_D
real(kind=dp), intent(in) :: D
real(kind=dp), intent(in) :: kappa
real(kind=dp), intent(in) :: c_A_SS
real(kind=dp), intent(in) :: c_B_SS
real(kind=dp), intent(out) :: Z

! internal variables
real(kind=dp) :: k_int_mf
real(kind=dp) :: k_int_mr
real(kind=dp) :: k_Mf
real(kind=dp) :: k_Mr
real(kind=dp), parameter :: pi = 4.*atan(1.d0)

! calculate intrinsic rates
k_int_mf = p_plus*(Rc**2)*sqrt(8.*pi*kBT/m)
k_int_mr = p_minus*(Rc**2)*sqrt(8.*pi*kBT/m)

! calculate prefactor
Z = ((k_int_mf*c_A_SS-k_int_mr*c_B_SS)*k_D) / ( 4.*pi*D* (k_D*(1.+kappa*Rc)+k_int_mf+k_int_mr) )

end subroutine conc_exp_prefactor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function derj_integrand_N(r)
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! external variables
real(kind=dp) :: kBT                    ! temperature
real(kind=dp) :: Rn                     ! radius of catalytic sphere
real(kind=dp) :: Ea                     ! interaction energy of A withs C sphere
real(kind=dp) :: Eb                     ! interaction energy of B withs C sphere

! internal variables
integer :: i
real(kind=dp) :: sigma_Rn               ! potential cutoff
real(kind=dp) :: r                      ! indepedent variable
real(kind=dp) :: derj_integrand_N

! read external variables
open(unit=10,file="input.txt",status="old",action="read")
do i =1,4
 read(10,*)
end do
read(10,*)kBT
read(10,*)
read(10,*)Rn
do i =1,6
 read(10,*)
end do
read(10,*)Ea
read(10,*)Eb
close(10)

sigma_Rn = Rn*2.**(1./6.)

derj_integrand_N = r * &
(dexp(-(1./kBT)*(4.*Eb*((Rn/r)**12-(Rn/r)**6+0.25))) - &
dexp(-(1./kBT)*(4.*Ea*((Rn/r)**12-(Rn/r)**6+0.25))))

end function derj_integrand_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module derived_quantities_mod
