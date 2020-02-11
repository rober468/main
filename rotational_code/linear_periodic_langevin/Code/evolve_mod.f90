module evolve_mod
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine evolve(timestep,numrotors,r_CN,D_R,kBT,rotor_CM,half_sim_box_len,sim_box_len,Rn,Rc,kappa,&
                  derj_len_N,Z,mean,var,F,w,angle)
use forces_mod
use random_mod
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)             ! single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! double precision

! external variables
real(kind=dp), intent(in) :: timestep
integer, intent(in) :: numrotors
real(kind=dp), intent(in) :: r_CN
real(kind=dp), intent(in) :: D_R
real(kind=dp), intent(in) :: kBT
real(kind=dp), intent(in) :: half_sim_box_len
real(kind=dp), intent(in) :: sim_box_len
real(kind=dp), intent(in) :: Rn
real(kind=dp), intent(in) :: Rc
real(kind=dp), intent(in) :: kappa
real(kind=dp), intent(in) :: derj_len_N
real(kind=dp), intent(in) :: Z
real(kind=dp), intent(in) :: mean
real(kind=dp), intent(in) :: var
real(kind=dp), dimension(:,:), intent(in) :: rotor_CM
real(kind=dp), dimension(:), intent(inout) :: F
real(kind=dp), dimension(:), intent(inout) :: w
real(kind=dp), dimension(:), intent(inout) :: angle

! internal variables
integer :: i                                                    ! loop variable
real(kind=dp) :: F_coeff                                        ! Coefficient of F
real(kind=dp) :: w_coeff                                        ! Coefficient of w
real(kind=dp), parameter :: pi=4.*atan(1.d0)


! derived quanitites
F_coeff = (r_CN*D_R*timestep)/kBT
w_coeff = sqrt(2.*D_R*timestep)

! compute angle
do i = 1,numrotors
 angle(i) = modulo( angle(i) + F_coeff*F(i) + w_coeff * w(i) , 2.*pi )
end do

! calculate force and random standard normal variable
call force(numrotors,rotor_CM,r_CN,angle,half_sim_box_len,sim_box_len,Rn,Rc,kappa,kBT,&
           derj_len_N,Z,F)

do i = 1,numrotors
 call gauss_rand(mean,var,w(i))
end do

end subroutine evolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module evolve_mod
