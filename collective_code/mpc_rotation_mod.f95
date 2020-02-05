module mpc_rotation_mod
implicit none

contains

subroutine mpc_rotation(NT,rotang,numcell,a0,Box,Xsol,Ysol,Zsol,VXsol,VYsol,VZsol,VXCM,VYCM,VZCM, &
                        tCellShift,NewCellN,UVX,UVY,UVZ,mt,mti)
use mtmod
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
real(kind=dp) , intent(in) :: rotang                            ! Rotation angle
integer(kind=dp) , intent(in) :: numcell                        ! Number of subcells
real(kind=dp) , intent(in) :: a0                                ! MPC cell size
real(kind=dp) , dimension(:) , intent(in) :: Box                ! Length of simulation box
real(kind=dp) , dimension(:) , intent(in) :: Xsol               ! Solvent positions, x
real(kind=dp) , dimension(:) , intent(in) :: Ysol               ! Solvent positions, y
real(kind=dp) , dimension(:) , intent(in) :: Zsol               ! Solvent positions, z
real(kind=dp) , dimension(:) , intent(inout) :: VXsol           ! Solvent positions, x
real(kind=dp) , dimension(:) , intent(inout) :: VYsol           ! Solvent positions, y
real(kind=dp) , dimension(:) , intent(inout) :: VZsol           ! Solvent positions, z
real(kind=dp) , dimension(:) , intent(inout) :: VXCM            ! Center of mass velocity, x
real(kind=dp) , dimension(:) , intent(inout) :: VYCM            ! Center of mass velocity, y
real(kind=dp) , dimension(:) , intent(inout) :: VZCM            ! Center of mass velocity, z
real(kind=dp) , dimension(:) , intent(inout) :: tCellShift      ! Tot.solv.in cellN after shift
integer(kind=dp) , dimension(:) , intent(inout) :: NewCellN     ! ith particle's new cell N
real(kind=dp) , dimension(:) , intent(inout) :: UVX             ! Unit vector, x
real(kind=dp) , dimension(:) , intent(inout) :: UVY             ! Unit vector, y
real(kind=dp) , dimension(:) , intent(inout) :: UVZ             ! Unit vector, z
integer , dimension(:,:) , intent(inout) :: mt
integer , dimension(:) , intent(inout) :: mti

! Definitions for internal variables
real(kind=dp) , parameter :: pi = 4.d0*datan(1.d0)              ! pi
integer(kind=dp) :: i,j                                         ! Loop variable
real(kind=dp) :: random                                         ! random number
real(kind=dp) :: sinrot                                         ! Sine of rotation angle
real(kind=dp) :: cosrot                                         ! Cosine of rotation angle
real(kind=dp) :: Xshift                                         ! Shift in x direction
real(kind=dp) :: Yshift                                         ! Shift in y direction
real(kind=dp) :: Zshift                                         ! Shift in z direction
real(kind=dp) :: Xsolshift                                      ! Shifted position, x
real(kind=dp) :: Ysolshift                                      ! Shifted position, y
real(kind=dp) :: Zsolshift                                      ! Shifted position, z
integer(kind=dp) :: CellN                                       ! Identify subcell number
real(kind=dp) :: x,y                                            ! Variables for random numbers
real(kind=dp) :: phi                                            ! Azimuthal angle
real(kind=dp) :: theta                                          ! Polar angle
real(kind=dp) :: dVXsolVXCM                                     ! Diff. betw. CM vel. and ith particle vel., x
real(kind=dp) :: dVYsolVYCM                                     ! Diff. betw. CM vel. and ith particle vel., y
real(kind=dp) :: dVZsolVZCM                                     ! Diff. betw. CM vel. and ith particle vel., z
real(kind=dp) :: dVXsolVXCMrot                                  ! Diff. betw. CM vel. and ith particle vel., x, rotated
real(kind=dp) :: dVYsolVYCMrot                                  ! Diff. betw. CM vel. and ith particle vel., y, rotated
real(kind=dp) :: dVZsolVZCMrot                                  ! Diff. betw. CM vel. and ith particle vel., z, rotated
integer , dimension(8) :: values
real(kind=dp) , dimension(:,:) , allocatable :: tcellshift_part
real(kind=dp) , dimension(:,:) , allocatable :: vxcm_part
real(kind=dp) , dimension(:,:) , allocatable :: vycm_part
real(kind=dp) , dimension(:,:) , allocatable :: vzcm_part
integer :: num_threads,thread_num
integer :: omp_get_num_threads,omp_get_thread_num

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(num_threads)
num_threads = omp_get_num_threads()
!$OMP END PARALLEL

allocate(tcellshift_part(num_threads,numcell))
allocate(vxcm_part(num_threads,numcell))
allocate(vycm_part(num_threads,numcell))
allocate(vzcm_part(num_threads,numcell))

!This subroutine performs the multiparticle collisions for the solvent particles.
!*********************************************************************************!
!******************************MPC************************************************!
!*********************************************************************************!

sinrot = dsin(rotang)
cosrot = dcos(rotang)

!**************************General_Collision**************************************!

!!!!!!Shift particles between [-a0/2,a0/2] and calculate CM pre-collision velocity
!Set CM pre-collision velocities and total particles in cell N to 0
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(vxcm_part,vycm_part,vzcm_part,tcellshift_part) &
!$OMP PRIVATE(i)
do i = 1,numcell
 VXCM_part(:,i) = 0.
 VYCM_part(:,i) = 0.
 VZCM_part(:,i) = 0.
 tCellshift_part(:,i) = 0
end do
!$OMP END PARALLEL DO

!Calculate random shift between [-a0/2,a0/2]
Xshift = grnd(mt(:,1),mti(1))
Xshift = a0*(Xshift-0.5d0)                                      !Rand. # betw. [-a0/2,a0/2], x
Yshift = grnd(mt(:,1),mti(1))
Yshift = a0*(Yshift-0.5d0)                                      !Rand. # betw. [-a0/2,a0/2], y
Zshift = grnd(mt(:,1),mti(1))
Zshift = a0*(Zshift-0.5d0)                                      !Rand. # betw. [-a0/2,a0/2], z

!Apply shift and PBC to ith solvent particle
if ( zshift < 0 ) then

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i,Xsolshift,Ysolshift,Zsolshift,CellN) &
!$OMP SHARED(NT,Xsol,Ysol,Zsol,Xshift,Yshift,Zshift,NewCellN,Box) &
!$OMP SCHEDULE(static)
do i = 1,NT
 Xsolshift = Xsol(i) + Xshift
 Xsolshift = Xsolshift - Box(1)*floor(Xsolshift/Box(1))         !Periodic BC in x direction 
 Ysolshift = Ysol(i) + Yshift
 Ysolshift = Ysolshift - Box(2)*floor(Ysolshift/Box(2))         !Periodic BC in y direction
 Zsolshift = Zsol(i) + Zshift

 !Find new cell N for solvent and find total solv. in new cell N
 CellN = dint( floor(Xsolshift)+Box(1)*floor(Ysolshift)+Box(1)*Box(2)*(floor(Zsolshift+1))+1 )
 NewCellN(i) = CellN
end do
!$OMP END PARALLEL DO

else if ( zshift >= 0 ) then

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i,Xsolshift,Ysolshift,Zsolshift,CellN) &
!$OMP SHARED(NT,Xsol,Ysol,Zsol,Xshift,Yshift,Zshift,NewCellN,Box) &
!$OMP SCHEDULE(static)
do i = 1,NT
 Xsolshift = Xsol(i) + Xshift
 Xsolshift = Xsolshift - Box(1)*floor(Xsolshift/Box(1))         !Periodic BC in x direction 
 Ysolshift = Ysol(i) + Yshift
 Ysolshift = Ysolshift - Box(2)*floor(Ysolshift/Box(2))         !Periodic BC in y direction
 Zsolshift = Zsol(i) + Zshift

 !Find new cell N for solvent and find total solv. in new cell N
 CellN = dint( floor(Xsolshift)+Box(1)*floor(Ysolshift)+Box(1)*Box(2)*(floor(Zsolshift))+1 )
 NewCellN(i) = CellN
end do
!$OMP END PARALLEL DO

end if

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(NT,newcelln,vxsol,vysol,vzsol,tcellshift_part,vxcm_part,vycm_part,vzcm_part) &
!$OMP PRIVATE(i,thread_num)
do i = 1,NT
 thread_num = omp_get_thread_num()+1
 tCellshift_part(thread_num,NewCellN(i)) = tCellshift_part(thread_num,NewCellN(i)) + 1
 !Sum the velocities in each new cell N; used next for center of mass velocity calc.
 VXCM_part(thread_num,NewCellN(i)) = VXCM_part(thread_num,NewCellN(i)) + VXsol(i)
 VYCM_part(thread_num,NewCellN(i)) = VYCM_part(thread_num,NewCellN(i)) + VYsol(i)
 VZCM_part(thread_num,NewCellN(i)) = VZCM_part(thread_num,NewCellN(i)) + VZsol(i)
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numcell,tcellshift,vxcm,vycm,vzcm,uvx,uvy,uvz) &
!$OMP SHARED(vxcm_part,vycm_part,vzcm_part,tcellshift_part,num_threads,mt,mti) &
!$OMP PRIVATE(i,j,phi,theta,thread_num,random)
do i = 1,numcell
 tcellshift(i) = tcellshift_part(1,i)
 vxcm(i) = vxcm_part(1,i)
 vycm(i) = vycm_part(1,i)
 vzcm(i) = vzcm_part(1,i)
 do j = 2,num_threads
  tcellshift(i) = tcellshift(i) + tcellshift_part(j,i)
  vxcm(i) = vxcm(i) + vxcm_part(j,i)
  vycm(i) = vycm(i) + vycm_part(j,i)
  vzcm(i) = vzcm(i) + vzcm_part(j,i)
 end do
 if ( tCellshift(i) /= 0 ) then
  VXCM(i) = VXCM(i)/tCellshift(i)
  VYCM(i) = VYCM(i)/tCellshift(i)
  VZCM(i) = VZCM(i)/tCellshift(i)
 end if

 !!!!!!Find the new velocities for ith particle (from Kapral, 2008, "long MPC article")
 !Calculate random unit vector in Cart.coord. from spherical polar coord. (r=1)
 thread_num = omp_get_thread_num()+1
 random = grnd(mt(:,thread_num),mti(thread_num))
 phi = 2.d0*pi*random                           !Random angle between 0 and 2*pi
 random = grnd(mt(:,thread_num),mti(thread_num))
 theta = 2.d0*random-1.d0
 UVX(i) = dsqrt(1.d0-theta*theta)*dcos(phi)     !x component of random unit vector
 UVY(i) = dsqrt(1.d0-theta*theta)*dsin(phi)     !y component of random unit vector
 UVZ(i) = theta                                 !z component of random unit vector
end do
!$OMP END PARALLEL DO
!Find the rotated comp. of the difference between CM velocity and particle i velocity

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE (i,dVXsolVXCM,dVYsolVYCM,dVZsolVZCM) &
!$OMP PRIVATE (dVXsolVXCMrot,dVYsolVYCMrot,dVZsolVZCMrot) &
!$OMP SHARED (NT,VXCM,VYCM,VZCM,UVX,UVY,UVZ,VXsol,VYsol,VZsol) &
!$OMP SHARED (NewCellN,sinrot,cosrot) &
!$OMP SCHEDULE(static)

do i = 1,NT

 ! Difference between CM velocity and particle i velocity, unrotated
 dVXsolVXCM=VXsol(i)-VXCM(NewCellN(i))
 dVYsolVYCM=VYsol(i)-VYCM(NewCellN(i))
 dVZsolVZCM=VZsol(i)-VZCM(NewCellN(i))

 ! Rotation of difference vector (Rodrigues' Rotation Formula)
 dVXsolVXCMrot = dVXsolVXCM * ( cosrot + UVX(NewCellN(i))**2*(1.d0-cosrot) ) + &
                 dVYsolVYCM * ( UVX(NewCellN(i))*UVY(NewCellN(i))*(1.d0-cosrot) - UVZ(NewCellN(i))*sinrot ) + &
                 dVZsolVZCM * ( UVX(NewCellN(i))*UVZ(NewCellN(i))*(1.d0-cosrot) + UVY(NewCellN(i))*sinrot )

 dVYsolVYCMrot = dVXsolVXCM * ( UVX(NewCellN(i))*UVY(NewCellN(i))*(1.d0-cosrot) + UVZ(NewCellN(i))*sinrot ) + &
                 dVYsolVYCM * ( cosrot + UVY(NewCellN(i))**2*(1.d0-cosrot) ) + &
                 dVZsolVZCM * ( UVY(NewCellN(i))*UVZ(NewCellN(i))*(1.d0-cosrot) - UVX(NewCellN(i))*sinrot )

 dVZsolVZCMrot = dVXsolVXCM * ( UVX(NewCellN(i))*UVZ(NewCellN(i))*(1.d0-cosrot) - UVY(NewCellN(i))*sinrot ) + &
                 dVYsolVYCM * ( UVY(NewCellN(i))*UVZ(NewCellN(i))*(1.d0-cosrot) + UVX(NewCellN(i))*sinrot ) + &
                 dVZsolVZCM * ( cosrot + UVZ(NewCellN(i))**2*(1.d0-cosrot) )

 !New velocities after collision
 VXsol(i)=VXCM(NewCellN(i))+dVXsolVXCMrot
 VYsol(i)=VYCM(NewCellN(i))+dVYsolVYCMrot
 VZsol(i)=VZCM(NewCellN(i))+dVZsolVZCMrot

end do

!$OMP END PARALLEL DO

deallocate(tcellshift_part)
deallocate(vxcm_part)
deallocate(vycm_part)
deallocate(vzcm_part)

end subroutine mpc_rotation

end module mpc_rotation_mod
