module forces_mod
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine force(numrotors,rotor_CM,r_CN,angle,half_sim_box_len,sim_box_len,Rn,Rc,kappa,kBT,&
                 derj_len_N,Z,F)
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! external variables
integer, intent(in) :: numrotors
real(kind=dp), dimension(:,:), intent(in) :: rotor_CM
real(kind=dp), intent(in) :: r_CN
real(kind=dp), dimension(:), intent(inout) :: angle
real(kind=dp), intent(in) :: half_sim_box_len
real(kind=dp), intent(in) :: sim_box_len
real(kind=dp), intent(in) :: Rn
real(kind=dp), intent(in) :: Rc
real(kind=dp), intent(in) :: kappa
real(kind=dp), intent(in) :: kBT
real(kind=dp), intent(in) :: derj_len_N
real(kind=dp), intent(in) :: Z
real(kind=dp), dimension(:) , intent(out) :: F

! internal variables
integer :: i,j
real(kind=dp), dimension(:,:,:), allocatable :: sphere_pos
real(kind=dp), dimension(:,:,:), allocatable :: d_ij_vec
real(kind=dp), dimension(:,:), allocatable :: d_ij
real(kind=dp), dimension(:,:), allocatable :: sin_theta_ij
real(kind=dp), dimension(:,:,:), allocatable :: s_ij_Rn
real(kind=dp), parameter :: pi = 4.*atan(1.d0)


! allocate arrays
allocate(sphere_pos(numrotors,2,2))
allocate(d_ij_vec(numrotors,numrotors,2))
allocate(d_ij(numrotors,numrotors))
allocate(sin_theta_ij(numrotors,numrotors))
allocate(s_ij_Rn(numrotors,numrotors,2))

! calculate sphere positions
do i = 1,numrotors,1
 ! N sphere
 sphere_pos(i,1,1) = rotor_CM(i,1) + r_CN * cos(angle(i))
 sphere_pos(i,1,2) = rotor_CM(i,2) + r_CN * sin(angle(i))
 ! C sphere
 sphere_pos(i,2,1) = rotor_CM(i,1) - r_CN * cos(angle(i))
 sphere_pos(i,2,2) = rotor_CM(i,2) - r_CN * sin(angle(i))
end do

! 1. calculate displacement between vectors d_ij_vec
! 2. introduce periodic condition
! 3. calculate distance d_ij
! 4. calculate (d_ii_vec x d_ij_vec)/(d_ii*d_ij) = sin_theta_ij \hat{n}
! 5. calculate limits for integral in f_ij
do i = 1,numrotors,1
 do j = 1,numrotors,1
  d_ij_vec(i,j,1) = sphere_pos(i,1,1) - sphere_pos(j,2,1)
  d_ij_vec(i,j,2) = sphere_pos(i,1,2) - sphere_pos(j,2,2)
  if ( d_ij_vec(i,j,1) <= -1.*half_sim_box_len ) then
   d_ij_vec(i,j,1) = d_ij_vec(i,j,1) + sim_box_len
  else if ( d_ij_vec(i,j,1) > half_sim_box_len ) then
   d_ij_vec(i,j,1) = d_ij_vec(i,j,1) - sim_box_len
  end if
  d_ij(i,j) = sqrt( d_ij_vec(i,j,1)**2 + d_ij_vec(i,j,2)**2 )
  !write(*,*)d_ij(i,j)
 end do
end do

! 1. sine of angle between z_ii and z_ij with correct sign for normal direction
! 2. s_{ij}^{+/-}(R_N) for limits of integration
do i = 1,numrotors,1
 do j = 1,numrotors,1
  sin_theta_ij(i,j) = ( &
d_ij_vec(i,i,1)*d_ij_vec(i,j,2)-d_ij_vec(i,i,2)*d_ij_vec(i,j,1) ) &
                      / ( d_ij(i,i)*d_ij(i,j) )
  s_ij_Rn(i,j,1) = sqrt( Rn**2 + d_ij(i,j)**2 + 2.*Rn*d_ij(i,j) )  ! plus
  s_ij_Rn(i,j,2) = sqrt( Rn**2 + d_ij(i,j)**2 - 2.*Rn*d_ij(i,j) )  ! minus
 end do
end do

! force on the ith N sphere due to the jth C sphere
F = 0.
do i = 1,numrotors,1
 do j = 1,numrotors,1
  if ( i /= j ) then
   F(i) = F(i) + (sin_theta_ij(i,j)/d_ij(i,j)**2 ) * ( &
           ((Rn**2+d_ij(i,j)**2-s_ij_Rn(i,j,1)**2)/kappa - 2.*s_ij_Rn(i,j,1)/(kappa)**2 - 2./(kappa)**3) &
            * exp(-1.*kappa*(s_ij_Rn(i,j,1)-Rc)) - &
           ((Rn**2+d_ij(i,j)**2-s_ij_Rn(i,j,2)**2)/kappa - 2.*s_ij_Rn(i,j,2)/(kappa)**2 - 2./(kappa)**3) &
            * exp(-1.*kappa*(s_ij_Rn(i,j,2)-Rc)) &
           )
  end if
 end do
 F(i) = ((3.*pi*kBT*derj_len_N*Z)/(Rn**2)) * F(i)
end do

! deallocate arrays
deallocate(sphere_pos)
deallocate(d_ij_vec)
deallocate(d_ij)
deallocate(sin_theta_ij)
deallocate(s_ij_Rn)

end subroutine force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module forces_mod
