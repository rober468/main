module initial_mod
implicit none

contains

subroutine initial(Box,numdims,Dcn,Na0,Nb0,Ns0,NT,cutoff,Rc,Rn,Rw,varA,varB,varS,mean,flag,Xdim,Ydim,Zdim, &
                   VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,Xsol,Ysol,Zsol,VXsol,VYsol,VZsol, &
                   FpairX,FpairY,FpairZ,motorCM,mt,mti)
use mtmod
use gauss_rand_mod
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)     ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)   ! Double precision

! Definitions for in/out variables
real(kind=dp) , dimension(:) , intent(in) :: Box        ! Length of simulation box;x=1,y=2,z=3
integer , intent(in) :: numdims                         ! Number of dimers
real(kind=dp) , intent(in) :: Dcn                       ! Internuclear distance
integer(kind=dp) , intent(in) :: Na0                    ! Initial number of particle A
integer(kind=dp) , intent(in) :: Nb0                    ! Initial number of particle B
integer(kind=dp) , intent(in) :: Ns0                    ! Initial number of particle S
integer(kind=sp) , intent(in) :: NT                     ! Total number of particles
real(kind=dp) , intent(in) :: cutoff                    ! Cut. rad. C sphere LJ pot., squared
real(kind=dp) , intent(in) :: Rc                        ! Radius C sphere
real(kind=dp) , intent(in) :: Rn                        ! Radius N sphere
real(kind=dp) , intent(in) :: Rw
real(kind=dp) , intent(in) :: varA                      ! Variance of A vel. Gaussian distr.
real(kind=dp) , intent(in) :: varB                      ! Variance of A vel. Gauss. distr.
real(kind=dp) , intent(in) :: varS                      ! Variance of S vel. Gauss. distr.
real(kind=dp) , intent(in) :: mean                      ! Mean of vel. Gaussian distr.
integer(kind=sp) , dimension(:), intent(out) :: flag    ! Identity of solvent
real(kind=dp) , dimension(:,:) , intent(out) :: Xdim    ! Dimer spheres positions, x
real(kind=dp) , dimension(:,:) , intent(out) :: Ydim    ! Dimer spheres positions, y
real(kind=dp) , dimension(:,:) , intent(out) :: Zdim    ! Dimer spheres positions, z
real(kind=dp) , dimension(:,:) , intent(out) :: VXdim   ! Dimer spheres velocities, x
real(kind=dp) , dimension(:,:) , intent(out) :: VYdim   ! Dimer spheres velocities, y
real(kind=dp) , dimension(:,:) , intent(out) :: VZdim   ! Dimer spheres velocities, z
real(kind=dp) , dimension(:,:) , intent(out) :: FXdim   ! Dimer spheres forces, x
real(kind=dp) , dimension(:,:) , intent(out) :: FYdim   ! Dimer spheres forces, y
real(kind=dp) , dimension(:,:) , intent(out) :: FZdim   ! Dimer spheres forces, z
real(kind=dp) , dimension(:) , intent(out) :: Xsol      ! Solvent positions, x
real(kind=dp) , dimension(:) , intent(out) :: Ysol      ! Solvent positions, y
real(kind=dp) , dimension(:) , intent(out) :: Zsol      ! Solvent positions, z
real(kind=dp) , dimension(:) , intent(out) :: VXsol     ! Solvent positions, x
real(kind=dp) , dimension(:) , intent(out) :: VYsol     ! Solvent positions, y
real(kind=dp) , dimension(:) , intent(out) :: VZsol     ! Solvent positions, z
real(kind=dp) , dimension(:,:,:) , intent(out) :: FpairX! Solvent positions, x
real(kind=dp) , dimension(:,:,:) , intent(out) :: FpairY! Solvent positions, y
real(kind=dp) , dimension(:,:,:) , intent(out) :: FpairZ! Solvent positions, z
real(kind=dp) , dimension(:,:) , intent(inout) :: motorCM
integer , dimension(:) , intent(inout) :: mt
integer , intent(inout) :: mti

! Definitions for internal variables
real(kind=dp) :: x,y,z,theta,phi                        ! Variables for random numbers
integer(kind=dp) :: i,j,k                               ! Loop variables
real(kind=dp) , parameter :: pi = 4.d0*datan(1.d0)      ! pi
real(kind=dp) :: dVXsol                                 ! Variable used for iSumdVXsol
real(kind=dp) :: dVYsol                                 ! Variable used for iSumdVYsol
real(kind=dp) :: dVZsol                                 ! Variable used for iSumdVZsol
real(KIND=dp) , dimension(numdims,3) :: dR              ! Diff. in dim. and solv. positions
real(kind=dp) :: sumVXsol                               ! Sum of initial velocities, x
real(kind=dp) :: sumVYsol                               ! Sum of initial velocities, y
real(kind=dp) :: sumVZsol                               ! Sum of initial velocities, z
real(kind=dp) :: cutoffC2                               ! Cut. C sphere LJ pot. squared
real(kind=dp) :: cutoffN2                               ! Cut. C sphere LJ pot. squared
real(kind=dp) :: wall_cutoffC2                          ! Cut. N sphere LJ pot. squared
real(kind=dp) :: wall_cutoffN2                          ! Cut. N sphere LJ pot. squared

!!!!!! Derived quantities
cutoffC2 = (Rc * cutoff)**2                             ! Cut. C sphere LJ pot. squared
cutoffN2 = (Rn * cutoff)**2                             ! Cut. N sphere LJ pot. squared
wall_cutoffC2 = ((Rw+Rc)*3.d0**(1./6.))**2
wall_cutoffN2 = ((Rw+Rn)*3.d0**(1./6.))**2

!!!!!! Initial position of dimer with random position and initial orientation
open(unit=10,file="motor_positions.dat",status="replace",position="append")
do i = 1,numdims
 400 x = grnd(mt,mti)
     y = grnd(mt,mti)
     z = grnd(mt,mti)
     phi = grnd(mt,mti)
     theta = grnd(mt,mti)
     phi = 2.d0*pi*phi
     theta = 2.d0*theta-1.d0
     Xdim(i,1) = box(1)*x + 0.5*Dcn*(dsqrt(1.d0-theta*theta)*dcos(phi))
     Ydim(i,1) = box(2)*y + 0.5*Dcn*(dsqrt(1.d0-theta*theta)*dsin(phi))
     Zdim(i,1) = box(3)*z + 0.5*Dcn*theta
     Xdim(i,2) = box(1)*x - 0.5*Dcn*(dsqrt(1.d0-theta*theta)*dcos(phi))
     Ydim(i,2) = box(2)*y - 0.5*Dcn*(dsqrt(1.d0-theta*theta)*dsin(phi))
     Zdim(i,2) = box(3)*z - 0.5*Dcn*theta
     motorCM(i,1) = 0.5 * (Xdim(i,1)+Xdim(i,2))
     motorCM(i,2) = 0.5 * (Ydim(i,1)+Ydim(i,2))
     motorCM(i,3) = 0.5 * (Zdim(i,1)+Zdim(i,2))
     Xdim(i,1) = Xdim(i,1) - Box(1)*floor(Xdim(i,1)/Box(1))
     Ydim(i,1) = Ydim(i,1) - Box(2)*floor(Ydim(i,1)/Box(2))
     Xdim(i,2) = Xdim(i,2) - Box(1)*floor(Xdim(i,2)/Box(1))
     Ydim(i,2) = Ydim(i,2) - Box(2)*floor(Ydim(i,2)/Box(2))
     if ( i == 1 ) then
      if ( Zdim(i,1) < sqrt(wall_cutoffC2)+Rc*cutoff .or. Zdim(i,1) > box(3)-sqrt(wall_cutoffC2)-Rc*cutoff &
      .or. Zdim(i,2) < sqrt(wall_cutoffN2)+Rn*cutoff .or. Zdim(i,2) > box(3)-sqrt(wall_cutoffN2)-Rn*cutoff ) then
       go to 400
      end if
     else if ( i /= 1 ) then
      do j = 1,i-1
       if ((Xdim(i,1)-Xdim(j,1))**2+(Ydim(i,1)-Ydim(j,1))**2+(Zdim(i,1)-Zdim(j,1))**2 <(2.*Rc*cutoff)**2 &
       .or. (Xdim(i,1)-Xdim(j,2))**2+(Ydim(i,1)-Ydim(j,2))**2+(Zdim(i,1)-Zdim(j,2))**2 < ((Rn+Rc)*cutoff)**2 &
       .or. (Xdim(i,2)-Xdim(j,1))**2+(Ydim(i,2)-Ydim(j,1))**2+(Zdim(i,2)-Zdim(j,1))**2 < ((Rc+Rn)*cutoff)**2 &
       .or. (Xdim(i,2)-Xdim(j,2))**2+(Ydim(i,2)-Ydim(j,2))**2+(Zdim(i,2)-Zdim(j,2))**2 < (2.*Rn*cutoff)**2 &
       .or. Zdim(i,1) < sqrt(wall_cutoffC2)+Rc*cutoff .or. Zdim(i,1) > box(3)-sqrt(wall_cutoffC2)-Rc*cutoff &
       .or. Zdim(i,2) < sqrt(wall_cutoffN2)+Rn*cutoff .or. Zdim(i,2) > box(3)-sqrt(wall_cutoffN2)-Rn*cutoff ) then
        goto 400
       end if
      end do
     end if
     write(10,*)motorCM(i,1),motorCM(i,2),motorCM(i,3)
end do
close(10)

!!!!!! Initial position of solvent A
do i = 1,Na0
 flag(i) = 0                                            ! Solvent particle is initially A

 100 x = grnd(mt,mti)
     y = grnd(mt,mti)
     z = grnd(mt,mti)

     Xsol(i) = Box(1) * x               !Random position of solvent in box, x
     Ysol(i) = Box(2) * y               !Random position of solvent in box, y
     Zsol(i) = Box(3) * z               !Random position of solvent in box, z

 ! Make sure A solvent isn't overlapping the C dimer sphere; if alright, go to N
 do k = 1,numdims
  dR(k,1) = Xsol(i) - Xdim(k,1)
  dR(k,2) = Ysol(i) - Ydim(k,1)
  dR(k,3) = Zsol(i) - Zdim(k,1)
 
  do j=1,2
   if ( dR(k,j) > 0.5d0*Box(j) ) then
    dR(k,j) = dR(k,j) - Box(j)
   else if ( dR(k,j) < -0.5d0*Box(j) ) then
    dR(k,j) = dR(k,j) + Box(j)
   end if
  end do

  if ( (dR(k,1)*dR(k,1) + dR(k,2)*dR(k,2) + dR(k,3)*dR(k,3)) < cutoffC2) then!If outside of C sphere cutoff distance, then
   goto 100                                                              !move on to check with N; if not, restart
  else

   ! Make sure A solvent isn't overlapping N dimer sphere; if alright,then position good
   dR(k,1) = Xsol(i) - Xdim(k,2)
   dR(k,2) = Ysol(i) - Ydim(k,2)
   dR(k,3) = Zsol(i) - Zdim(k,2)

   do j=1,2
    if ( dR(k,j) > 0.5d0*Box(j) ) then
     dR(k,j) = dR(k,j) - Box(j)
    else if ( dR(k,j) < -0.5d0*Box(j) ) then
     dR(k,j) = dR(k,j) + Box(j)
    end if
   end do

   if ( (dR(k,1)*dR(k,1) + dR(k,2)*dR(k,2) + dR(k,3)*dR(k,3)) < cutoffN2 ) then!If outside of C sphere cutoff distance, then
    goto 100                                                             !move on to check with N; if not, restart
   end if
  end if
 end do
end do

!!!!!!Initial position of B
do i = Na0+1,Na0+Nb0
 flag(i) = 1                    !Solvent particle is initially B

 200 x = grnd(mt,mti)
     y = grnd(mt,mti)
     z = grnd(mt,mti)

     Xsol(i) = Box(1)*x         !Random position of solvent in box, x
     Ysol(i) = Box(2)*y         !Random position of solvent in box, y
     Zsol(i) = Box(3)*z         !Random position of solvent in box, z

 !Make sure B solvent isn't overlapping the C dimer sphere; if alright, go to N
 do k = 1,numdims
  dR(k,1) = Xsol(i) - Xdim(k,1)
  dR(k,2) = Ysol(i) - Ydim(k,1)
  dR(k,3) = Zsol(i) - Zdim(k,1)
 
  do j=1,2
   if ( dR(k,j) > 0.5d0*Box(j) ) then
    dR(k,j) = dR(k,j) - Box(j)
   else if ( dR(k,j) < -0.5d0*Box(j) ) then
    dR(k,j) = dR(k,j) + Box(j)
   end if
  end do

  if ( (dR(k,1)*dR(k,1) + dR(k,2)*dR(k,2) + dR(k,3)*dR(k,3)) < cutoffC2) then!If outside of C sphere cutoff distance, then
   goto 200                                                              !move on to check with N; if not, restart
  else

   ! Make sure B solvent isn't overlapping N dimer sphere; if alright,then position good
   dR(k,1) = Xsol(i) - Xdim(k,2)
   dR(k,2) = Ysol(i) - Ydim(k,2)
   dR(k,3) = Zsol(i) - Zdim(k,2)

   do j=1,2
    if ( dR(k,j) > 0.5d0*Box(j) ) then
     dR(k,j) = dR(k,j) - Box(j)
    else if ( dR(k,j) < -0.5d0*Box(j) ) then
     dR(k,j) = dR(k,j) + Box(j)
    end if
   end do

   if ( (dR(k,1)*dR(k,1) + dR(k,2)*dR(k,2) + dR(k,3)*dR(k,3)) < cutoffN2 ) then!If outside of C sphere cutoff distance, then
    goto 200                                                             !move on to check with N; if not, restart
   end if
  end if
 end do
end do

!!!!!!Initial position of S
do i = Na0+Nb0+1,NT
 flag(i) = 2                    !Solvent particle is initially B

 300 x = grnd(mt,mti)
     y = grnd(mt,mti)
     z = grnd(mt,mti)

     Xsol(i) = Box(1)*x         !Random position of solvent in box, x
     Ysol(i) = Box(2)*y         !Random position of solvent in box, y
     Zsol(i) = Box(3)*z         !Random position of solvent in box, z

 !Make sure B solvent isn't overlapping the C dimer sphere; if alright, go to N
 do k = 1,numdims
  dR(k,1) = Xsol(i) - Xdim(k,1)
  dR(k,2) = Ysol(i) - Ydim(k,1)
  dR(k,3) = Zsol(i) - Zdim(k,1)

  do j=1,2
   if ( dR(k,j) > 0.5d0*Box(j) ) then
    dR(k,j) = dR(k,j) - Box(j)
   else if ( dR(k,j) < -0.5d0*Box(j) ) then
    dR(k,j) = dR(k,j) + Box(j)
   end if
  end do

  if ( (dR(k,1)*dR(k,1) + dR(k,2)*dR(k,2) + dR(k,3)*dR(k,3)) < cutoffC2) then!If outside of C sphere cutoff distance, then
   goto 300                                                              !move on to check with N; if not, restart
  else

   ! Make sure B solvent isn't overlapping N dimer sphere; if alright,then
   ! position good
   dR(k,1) = Xsol(i) - Xdim(k,2)
   dR(k,2) = Ysol(i) - Ydim(k,2)
   dR(k,3) = Zsol(i) - Zdim(k,2)

   do j=1,2
    if ( dR(k,j) > 0.5d0*Box(j) ) then
     dR(k,j) = dR(k,j) - Box(j)
    else if ( dR(k,j) < -0.5d0*Box(j) ) then
     dR(k,j) = dR(k,j) + Box(j)
    end if
   end do

   if ( (dR(k,1)*dR(k,1) + dR(k,2)*dR(k,2) + dR(k,3)*dR(k,3)) < cutoffN2 ) then!If outside of C sphere cutoff distance, then
    goto 300                                                             !move on to check with N; if not, restart
   end if
  end if
 end do
end do


!!!!!!Initial velocity of dimer
VXdim = 0.
VYdim = 0.
VZdim = 0.

sumVXsol = 0.
sumVYsol = 0.
sumVZsol = 0.

!!!!!!Initial velocity of solvent A
!For x direction
do i = 1,Na0                                    ! Rand. A velocity from Gauss. distr.,x
 call gauss_rand(mean,varA,VXsol(i),mt,mti)
 sumVXsol = sumVXsol + VXsol(i)
end do
!For y direction
do i = 1,Na0                                    ! Rand. A velocity from Gauss. distr.,y
 call gauss_rand(mean,varA,VYsol(i),mt,mti)
 sumVYsol = sumVYsol + VYsol(i)
end do
!For z direction
do i = 1,Na0                                    ! Rand. A velocity from Gauss. distr.,z
 call gauss_rand(mean,varA,VZsol(i),mt,mti)
 sumVZsol = sumVZsol + VZsol(i)
end do

!!!!!!Initial velocity of B solvent
!For x direction
do i = Na0+1,Na0+Nb0                            ! Rand. B velocity from Gauss. distr.,x
 call gauss_rand(mean,varB,VXsol(i),mt,mti)
 sumVXsol = sumVXsol + VXsol(i)
end do
!For y direction
do i = Na0+1,Na0+Nb0                            ! Rand. B velocity from Gauss. distr.,y
 call gauss_rand(mean,varB,VYsol(i),mt,mti)
 sumVYsol = sumVYsol + VYsol(i)
end do
!For z direction
do i = Na0+1,Na0+Nb0                            ! Rand. B velocity from Gauss. distr.,z
 call gauss_rand(mean,varB,VZsol(i),mt,mti)
 sumVZsol = sumVZsol + VZsol(i)
end do

!!!!!!Initial velocity of S solvent
!For x direction
do i = Na0+Nb0+1,NT                             ! Rand. S velocity from Gauss. distr.,x
 call gauss_rand(mean,varS,VXsol(i),mt,mti)
 sumVXsol = sumVXsol + VXsol(i)
end do
!For y direction
do i = Na0+Nb0+1,NT                             ! Rand. S velocity from Gauss. distr.,y
 call gauss_rand(mean,varS,VYsol(i),mt,mti)
 sumVYsol = sumVYsol + VYsol(i)
end do
!For z direction
do i = Na0+Nb0+1,NT                             ! Rand. S velocity from Gauss. distr.,z
 call gauss_rand(mean,varS,VZsol(i),mt,mti)
 sumVZsol = sumVZsol + VZsol(i)
end do

!!!!!!Make total momentum zero
! Average velocity in each direction
sumVXsol = sumVXsol/dble(NT)
sumVYsol = sumVYsol/dble(NT)
sumVZsol = sumVZsol/dble(NT)

! Shift and adjust
do i = 1,NT
 VXsol(i) = VXsol(i) - sumVXsol
 VYsol(i) = VYsol(i) - sumVYsol
 VZsol(i) = VZsol(i) - sumVZsol
end do

!!!!!!Initial force on solvent and dimer
FpairX=0.                                       ! Initial force on solvent is zero
FpairY=0.
FpairZ=0.
FXdim=0.                                        ! Initial force on dimer is zero
FYdim=0.
FZdim=0.

end subroutine initial

end module initial_mod
