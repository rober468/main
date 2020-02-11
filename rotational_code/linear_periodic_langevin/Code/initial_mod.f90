module initial_mod
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize(seed,numrotors,r_CN,len_rotor,half_sim_box_len,sim_box_len,Rn,Rc,kappa,kBT,&
                      derj_len_N,Z,mean,var,sep,cstep,rotor_CM,angle,F,w)
use mt95
use random_mod
use forces_mod
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! external variables
integer, intent(in) :: seed
integer, intent(in) :: numrotors
real(kind=dp), intent(in) :: r_CN
real(kind=dp), intent(in) :: len_rotor
real(kind=dp), intent(in) :: half_sim_box_len
real(kind=dp), intent(in) :: sim_box_len
real(kind=dp), intent(in) :: Rn
real(kind=dp), intent(in) :: Rc
real(kind=dp), intent(in) :: kappa
real(kind=dp), intent(in) :: kBT
real(kind=dp), intent(in) :: derj_len_N
real(kind=dp), intent(in) :: Z
real(kind=dp), intent(in) :: mean
real(kind=dp), intent(in) :: var
real(kind=dp), intent(in) :: sep
integer, intent(out) :: cstep
real(kind=dp), dimension(:,:), intent(inout) :: rotor_CM
real(kind=dp), dimension(:), intent(out) :: angle
real(kind=dp), dimension(:) , intent(out) :: F
real(kind=dp), dimension(:) , intent(out) :: w

! internal variables
integer :: i                                  ! loop variables

! initial timestep
cstep = 0

! initialize random number generator
call genrand_init( put=seed )

! random initial angle
do i = 1,numrotors,1
 call circumference_2d_drand( angle(i) )
end do

! calculate centre-of-mass of each rotor
do i = 1,numrotors
 rotor_CM(i,1) = 0.5*(sep+len_rotor) + (i-1)*(sep+len_rotor)
 rotor_CM(i,2) = 0.
end do

! calculate force at time zero
call force(numrotors,rotor_CM,r_CN,angle,half_sim_box_len,sim_box_len,Rn,Rc,kappa,kBT,&
           derj_len_N,Z,F)

! calculate standard normal noise at time zero
do i =1,numrotors,1
 call gauss_rand(mean,var,w(i))
end do

end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module initial_mod
