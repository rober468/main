module evolve_mod
implicit none

contains

subroutine evolve(cstep,time,NT,MDt,MPCtstep,numdims,Ma,Mb,Mc,Mn,Dcn,Box,cutoff,Rc,Rn,EaC,EbC,EaN,EbN,&
                zdimeq,k_harm,RPab,RPba,sphere,flag,type_B,xsol,Ysol,Zsol,VXsol,VYsol,VZsol,Xdim,&
                  Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,motorCM,TPotential,Potential,FpairX,&
                  FpairY,FpairZ,random,skinrad,max_solvent_dist,neighbour_num,neighbour,time_count,&
                  time_penetrate,outside_particle_num,outside_particle,inside_particle_num,inside_particle,&
                  vol,sort_step,values)
use sort_hilbert_mod
use neighbour_list_mod
use forces_mod
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer , dimension(8) :: values 
integer(kind=sp) :: cstep                                       ! Current time step
integer(kind=dp) :: time                                        ! Time loop variable
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
real(kind=dp) , intent(in) :: MDt                               ! MD time step
integer(kind=sp) :: MPCtstep                                    ! MPC step happens every nth step
integer , intent(in) :: numdims
real(kind=dp) , intent(in) :: Ma                                ! Mass of one A particle
real(kind=dp) , intent(in) :: Mb                                ! Mass of one B particle
real(kind=dp) , intent(in) :: Mc                                ! Mass of C sphere
real(kind=dp) , intent(in) :: Mn                                ! Mass of N sphere
real(kind=dp) , intent(in) :: Dcn                               ! Internuclear bond length
real(kind=dp) , dimension(:), intent(in) :: Box                 ! Length of simulation box
real(kind=dp) , intent(in) :: cutoff                            ! Cut. rad. C sphere LJ pot., squared
real(kind=dp) , intent(in) :: Rc                                ! Radius C sphere
real(kind=dp) , intent(in) :: Rn                                ! Radius N sphere
real(kind=dp) , intent(in) :: EaC                               ! Reduced LJ Potential of A with C
real(kind=dp) , intent(in) :: EbC                               ! Reduced LJ Potential of B with C
real(kind=dp) , intent(in) :: EaN                               ! Reduced LJ Potential of A with N
real(kind=dp) , intent(in) :: EbN                               ! Reduced LJ Potential of B with N
real(kind=dp) , intent(in) :: zdimeq                            ! Equil. pos. of mon., z dir.
real(kind=dp) , intent(in) :: k_harm                            ! Force constant
real(kind=dp) , intent(in) :: RPab                              ! Probability of reaction from F to G
real(kind=dp) , intent(in) :: RPba                              ! Probability of reaction from G to F
integer , dimension(:) , intent(in) :: sphere                   ! Sphere identity; C=0;N=1
integer(kind=sp) , dimension(:), intent(inout) :: flag          ! Identity of solvent
integer(kind=sp) , dimension(:), intent(inout) :: type_B
real(kind=dp) , dimension(:) , intent(inout) :: Xsol            ! Force on solvent in x direction
real(kind=dp) , dimension(:) , intent(inout) :: Ysol            ! Force on solvent in y direction
real(kind=dp) , dimension(:) , intent(inout) :: Zsol            ! Force on solvent in z direction
real(kind=dp) , dimension(:,:) , intent(inout) :: Xdim          ! Force on solvent in x direction
real(kind=dp) , dimension(:,:) , intent(inout) :: Ydim          ! Force on solvent in y direction
real(kind=dp) , dimension(:,:) , intent(inout) :: Zdim          ! Force on solvent in z direction
real(kind=dp) , dimension(:) , intent(inout) :: VXsol           ! Force on solvent in x direction
real(kind=dp) , dimension(:) , intent(inout) :: VYsol           ! Force on solvent in y direction
real(kind=dp) , dimension(:) , intent(inout) :: VZsol           ! Force on solvent in z direction
real(kind=dp) , dimension(:,:) , intent(inout) :: VXdim         ! Force on solvent in x direction
real(kind=dp) , dimension(:,:) , intent(inout) :: VYdim         ! Force on solvent in y direction
real(kind=dp) , dimension(:,:) , intent(inout) :: VZdim         ! Force on solvent in z direction
real(kind=dp) , dimension(:,:) , intent(inout) :: FXdim         ! Force on solvent in x direction
real(kind=dp) , dimension(:,:) , intent(inout) :: FYdim         ! Force on solvent in y direction
real(kind=dp) , dimension(:,:) , intent(inout) :: FZdim         ! Force on solvent in z direction
real(kind=dp) , dimension(numdims,3) , intent(in) :: motorCM
real(kind=dp) , intent(out) :: TPotential                       ! Total potential
real(kind=dp) , dimension(:,:,:) , intent(inout) :: Potential   ! Potential per particle
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairX      ! LJ force in x direction
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairY      ! LJ force in y direction
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairZ      ! LJ force in z direction
real(kind=dp) , dimension(:) , intent(inout) :: random          ! Vector of random numbers
integer(kind=dp) , dimension(:,:) , intent(inout) :: neighbour_num! Number of solvent particles in list
integer(kind=dp) , dimension(:,:,:) , intent(inout) :: neighbour! Identify solvent number in list
real(kind=dp) , intent(in) :: skinrad                           ! Skin radius for neighbour list
real(kind=dp) , dimension(:) , intent(inout) :: max_solvent_dist! Max distance of solvent
integer , intent(inout) :: time_count
integer , intent(inout) :: time_penetrate
integer(kind=dp) , intent(inout) :: outside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: outside_particle
integer(kind=dp) :: vol                                         ! Volume of simulation box
integer :: sort_step
integer(kind=dp) , intent(inout) :: inside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: inside_particle

! Definitions for internal variables
integer(kind=dp) :: i,j                                         ! Loop variable
integer :: k                                                    ! Loop variable
real(kind=dp) :: hIMaMDt                                        ! Value used a lot for calculations
real(kind=dp) :: hIMbMDt                                        ! Value used a lot for calculations
real(kind=dp) :: hIMcMDt                                        ! Value used a lot for calculations
real(kind=dp) :: hIMnMDt                                        ! Value used a lot for calculations
real(kind=dp) :: time_aft_coll                                  ! Time after collision
real(kind=dp) :: IIMcIMn                                        ! Value used for first Lagr. mult.	
real(kind=dp) , dimension(:) , allocatable :: prXdim            ! Prev. force on solvent in x direction
real(kind=dp) , dimension(:) , allocatable :: prYdim            ! Prev. force on solvent in y direction
real(kind=dp) , dimension(:) , allocatable :: prZdim            ! Prev. force on solvent in z direction
real(kind=dp) , dimension(:) , allocatable :: diffdimP          ! Diff. in dimer position at prev. timestep
real(kind=dp) , dimension(:) , allocatable :: diffdim           ! Diff. in dimer position at curr. timestep
integer(kind=dp) :: Con1Acc                                     ! Number of iter. for first constraint
integer(kind=dp) , parameter :: Con1Max = 200000000             ! Max numb. of iter. for first constr.
real(kind=dp) , parameter :: Con1Conv = 1.d-6                   ! Tolerance from bond length
real(kind=dp) :: Con1diff                                       ! Difference of new and old bond length
real(kind=dp) , parameter :: Con1Convx = 1.d-6                  ! Tolerance from bond length
real(kind=dp) :: Con1diffx                                      ! Difference of new and old bond length
real(kind=dp) , parameter :: Con1Convy = 1.d-6                  ! Tolerance from bond length
real(kind=dp) :: Con1diffy                                      ! Difference of new and old bond length
real(kind=dp) , parameter :: Con1Convz = 1.d-6                  ! Tolerance from bond length
real(kind=dp) :: Con1diffz                                      ! Difference of new and old bond length
real(kind=dp) :: dVXdim                                         ! Difference in dimer velocities, x
real(kind=dp) :: dVYdim                                         ! Difference in dimer velocities, y
real(kind=dp) :: dVZdim                                         ! Difference in dimer velocities, z
real(kind=dp) , dimension(numdims,2) :: QXdim                   ! The q(i) in Andersen paper, x
real(kind=dp) , dimension(numdims,2) :: QYdim                   ! The q(i) in Andersen paper, y
real(kind=dp) , dimension(numdims,2) :: QZdim                   ! The q(i) in Andersen paper, z
real(kind=dp) , dimension(3) :: dimCM
integer(kind=dp) :: Con2Acc                                     ! Number of iter. for second constraint
integer(kind=dp),PARAMETER :: Con2Max=200000000                 ! Max numb. of iter. for first constr.
real(kind=dp),PARAMETER :: Con2Conv=1.d-8                       ! Tolerance from dot prod. constr.
real(kind=dp) :: Con2diff                                       ! Difference of new and old velocities
real(kind=dp),PARAMETER :: Con2Convx=1.d-8                      ! Tolerance from dot prod. constr.
real(kind=dp) :: Con2diffx                                      ! Difference of new and old velocities
real(kind=dp),PARAMETER :: Con2Convy=1.d-8                      ! Tolerance from dot prod. constr.
real(kind=dp) :: Con2diffy                                      ! Difference of new and old velocities
real(kind=dp),PARAMETER :: Con2Convz=1.d-8                      ! Tolerance from dot prod. constr.
real(kind=dp) :: Con2diffz                                      ! Difference of new and old velocities
real(kind=dp) :: G                                          
real(kind=dp) :: H
real(kind=dp) :: Hx
real(kind=dp) :: Hy
real(kind=dp) :: Hz
real(kind=dp) :: O                                          
real(kind=dp) :: P                                          
real(kind=dp) :: Px                                          
real(kind=dp) :: Py                                         
real(kind=dp) :: Pz
real(kind=dp) :: Rn6                                            ! Rn^6
real(kind=dp) :: Rn12                                           ! Rn^12
real(kind=dp) :: r12Rn12                                        ! 12*Rn^12
real(kind=dp) :: r6Rn6                                          ! 6*Rn^6
real(kind=dp) :: Rc6                                            ! Rc^6
real(kind=dp) :: Rc12                                           ! Rc^12
real(kind=dp) :: r12Rc12                                        ! 12*Rc^12
real(kind=dp) :: r6Rc6                                          ! 6*Rc^6
real(kind=dp) :: cutoffC                                        ! Cut. C sphere LJ pot.
real(kind=dp) :: cutoffN                                        ! Cut. N sphere LJ pot.
real(kind=dp) :: cutoffC2                                       ! Cut. C sphere LJ pot. squared
real(kind=dp) :: cutoffN2                                       ! Cut. N sphere LJ pot. squared
real(kind=dp) , dimension(:,:) , allocatable :: fwallz          ! Force of wall on monomer, z
real(kind=dp) , dimension(:,:) , allocatable :: wall_pot        ! Potential of wall with monomer
real(kind=dp) , dimension(:,:) , allocatable :: TPotential_part
real(kind=dp) :: maxv

allocate(prXdim(2),prYdim(2),prZdim(2),diffdimP(3),diffdim(3))
allocate(fwallz(numdims,2),wall_pot(numdims,2),TPotential_part(numdims,2))

!!!!!! Define variables within subroutine
hIMaMDt = (0.5d0*MDt)/Ma                                        ! Value used a lot for calculations
hIMbMDt = (0.5d0*MDt)/Mb                                        ! Value used a lot for calculations
hIMcMDt = (0.5d0*MDt)/Mc                                        ! Value used a lot for calculations
hIMnMDt = (0.5d0*MDt)/Mn                                        ! Value used a lot for calculations
IIMcIMn = 1.d0/((1.d0/Mc)+(1.d0/Mn))                            ! Value used for first Lagr. mult.
Rn6 = Rn**6
Rn12 = Rn**12
r6Rn6 = 6.d0*Rn6
r12Rn12 = 12.d0*Rn12
Rc6 = Rc**6
Rc12 = Rc**12
r6Rc6 = 6.d0*Rc6
r12Rc12 = 12.d0*Rc12
cutoffC = Rc * cutoff
cutoffN = Rn * cutoff
cutoffC2 = cutoffC**2                                           ! Cut. C sphere LJ pot. squared
cutoffN2 = cutoffN**2                                           ! Cut. N sphere LJ pot. squared

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Calculations before force calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Motion of the solvent molecules
do k = 1,numdims
 do i = 1,2
 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(hIMaMDt,hIMbMDt) &
 !$OMP SHARED(VXsol,VYsol,VZsol,FpairX,FpairY,FpairZ,flag) &
 !$OMP SHARED(neighbour,neighbour_num,i,k) &
 !$OMP PRIVATE(j) SCHEDULE(static)

!do k = 1,numdims
! do i = 1,2
  do j = 1,neighbour_num(k,i)
   if ( flag(neighbour(k,i,j)) == 0 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMaMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMaMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMaMDt * FpairZ(k,i,j)
   else if ( flag(neighbour(k,i,j)) == 1 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMbMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMbMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMbMDt * FpairZ(k,i,j)
   end if
  end do
! end do
!end do
!$OMP END PARALLEL DO
end do
end do

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(cstep,NT,MDt,Box,Xsol,Ysol,Zsol,flag,max_solvent_dist) &
!$OMP SHARED(VXsol,VYsol,VZsol,inside_particle_num,inside_particle) & 
!$OMP PRIVATE(i,time_aft_coll) SCHEDULE(static)

 ! Evolve inside particles to current time step
 do i = 1,inside_particle_num
  Xsol(inside_particle(i)) = Xsol(inside_particle(i)) + VXsol(inside_particle(i))*MDt
  if ( Xsol(inside_particle(i)) >= Box(1) ) then
   Xsol(inside_particle(i)) = Xsol(inside_particle(i)) - Box(1)
  else if ( Xsol(inside_particle(i)) < 0. ) then
   Xsol(inside_particle(i)) = Xsol(inside_particle(i)) + Box(1)
  end if
  Ysol(inside_particle(i)) = Ysol(inside_particle(i)) + VYsol(inside_particle(i))*MDt
  if ( Ysol(inside_particle(i)) >= Box(2) ) then
   Ysol(inside_particle(i)) = Ysol(inside_particle(i)) - Box(2)
  else if ( Ysol(inside_particle(i)) < 0. ) then
   Ysol(inside_particle(i)) = Ysol(inside_particle(i)) + Box(2)
  end if
  Zsol(inside_particle(i)) = Zsol(inside_particle(i)) + VZsol(inside_particle(i))*MDt

  ! Bounce back from wall
  if ( Zsol(inside_particle(i)) < 0. ) then
   time_aft_coll = Zsol(inside_particle(i)) / VZsol(inside_particle(i))
   VXsol(inside_particle(i)) = -1. * VXsol(inside_particle(i))
   VYsol(inside_particle(i)) = -1. * VYsol(inside_particle(i))
   VZsol(inside_particle(i)) = -1. * VZsol(inside_particle(i))
   Xsol(inside_particle(i)) = Xsol(inside_particle(i)) + 2. * VXsol(inside_particle(i)) * time_aft_coll
   Ysol(inside_particle(i)) = Ysol(inside_particle(i)) + 2. * VYsol(inside_particle(i)) * time_aft_coll
   Zsol(inside_particle(i)) = Zsol(inside_particle(i)) + 2. * VZsol(inside_particle(i)) * time_aft_coll
   if ( Xsol(inside_particle(i)) >= Box(1) ) then
    Xsol(inside_particle(i)) = Xsol(inside_particle(i)) - Box(1)
   else if ( Xsol(inside_particle(i)) < 0. ) then
    Xsol(inside_particle(i)) = Xsol(inside_particle(i)) + Box(1)
   end if
   if ( Ysol(inside_particle(i)) >= Box(2) ) then
    Ysol(inside_particle(i)) = Ysol(inside_particle(i)) - Box(2)
   else if ( Ysol(inside_particle(i)) < 0. ) then
    Ysol(inside_particle(i)) = Ysol(inside_particle(i)) + Box(2)
   end if
  else if ( Zsol(inside_particle(i)) > box(3) ) then
   time_aft_coll = (Zsol(inside_particle(i)) - box(3)) / VZsol(inside_particle(i))
   VXsol(inside_particle(i)) = -1. * VXsol(inside_particle(i))
   VYsol(inside_particle(i)) = -1. * VYsol(inside_particle(i))
   VZsol(inside_particle(i)) = -1. * VZsol(inside_particle(i))
   Xsol(inside_particle(i)) = Xsol(inside_particle(i)) + 2. * VXsol(inside_particle(i)) * time_aft_coll
   Ysol(inside_particle(i)) = Ysol(inside_particle(i)) + 2. * VYsol(inside_particle(i)) * time_aft_coll
   Zsol(inside_particle(i)) = Zsol(inside_particle(i)) + 2. * VZsol(inside_particle(i)) * time_aft_coll
   if ( Xsol(inside_particle(i)) >= Box(1) ) then
    Xsol(inside_particle(i)) = Xsol(inside_particle(i)) - Box(1)
   else if ( Xsol(inside_particle(i)) < 0. ) then
    Xsol(inside_particle(i)) = Xsol(inside_particle(i)) + Box(1)
   end if
   if ( Ysol(inside_particle(i)) >= Box(2) ) then
    Ysol(inside_particle(i)) = Ysol(inside_particle(i)) - Box(2)
   else if ( Ysol(inside_particle(i)) < 0. ) then
    Ysol(inside_particle(i)) = Ysol(inside_particle(i)) + Box(2)
   end if
  end if
  max_solvent_dist(inside_particle(i)) = max_solvent_dist(inside_particle(i)) + &
  (VXsol(inside_particle(i))**2+VYsol(inside_particle(i))**2+VZsol(inside_particle(i))**2)*MDt**2
 end do
!$OMP END PARALLEL DO

!!!!!! Motion of the dimer before force calculation using RATTLE

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(Xdim,Ydim,Zdim,QXdim,QYdim,QZdim,VXdim,VYdim,VZdim,hIMcMDt,hIMnMDt) &
!$OMP SHARED(FXdim,FYdim,FZdim,MDt,box,Mc,Mn,Dcn,motorCM,IIMcIMn,numdims) &
!$OMP PRIVATE(k,i,prXdim,prYdim,prZdim,diffdimP,Con1Acc,diffdim,dimCM,Con1diff) &
!$OMP PRIVATE(Con1diffx,Con1diffy,Con1diffz,G,Hx,Hy,Hz)

do k = 1,numdims
!!!!!! Set old dimer pos. and vel. to previous time step
 do i = 1,2
  prXdim(i) = Xdim(k,i)
  prYdim(i) = Ydim(k,i)
  prZdim(i) = Zdim(k,i)
 end do
 ! Difference in position at prev. step
 diffdimP(1) = Xdim(k,1) - Xdim(k,2)
 diffdimP(2) = Ydim(k,1) - Ydim(k,2)
 diffdimP(3) = Zdim(k,1) - Zdim(k,2)

 ! Periodic BC for diff. in pos. at pstep
 do i = 1,2
  if ( diffdimP(i) > 0.5d0*Box(i) ) then
   diffdimP(i) = diffdimP(i) - Box(i)
  else if ( diffdimP(i) < -0.5d0*Box(i) ) then
   diffdimP(i) = diffdimP(i) + Box(i)
  end if
 end do

 ! The value for q without added Lagrange multiplier term
 QXdim(k,1) = VXdim(k,1) + hIMcMDt*FXdim(k,1)
 QYdim(k,1) = VYdim(k,1) + hIMcMDt*FYdim(k,1)
 QZdim(k,1) = VZdim(k,1) + hIMcMDt*FZdim(k,1)
 QXdim(k,2) = VXdim(k,2) + hIMnMDt*FXdim(k,2)
 QYdim(k,2) = VYdim(k,2) + hIMnMDt*FYdim(k,2)
 QZdim(k,2) = VZdim(k,2) + hIMnMDt*FZdim(k,2)

 ! Advance position by one time step
 do i=1,2
  Xdim(k,i) = Xdim(k,i) + MDt*QXdim(k,i)
  Ydim(k,i) = Ydim(k,i) + MDt*QYdim(k,i)
  Zdim(k,i) = Zdim(k,i) + MDt*QZdim(k,i)
  ! Apply periodic BC (if particle exits sim. box, bring it back into box)
  Xdim(k,i) = Xdim(k,i) - Box(1)*floor(Xdim(k,i)/Box(1))
  Ydim(k,i) = Ydim(k,i) - Box(2)*floor(Ydim(k,i)/Box(2))
 end do

 Con1Acc=0
 do
  ! Difference in dimer positions
  diffdim(1) = Xdim(k,1)-Xdim(k,2)
  diffdim(2) = Ydim(k,1)-Ydim(k,2)
  diffdim(3) = Zdim(k,1)-Zdim(k,2)
  ! Center of mass
  dimCM(1) = 0.5*( Xdim(k,1) + Xdim(k,2) )
  dimCM(2) = 0.5*( Ydim(k,1) + Ydim(k,2) )
  dimCM(3) = 0.5*( Zdim(k,1) + Zdim(k,2) )

 !Apply periodic conditions
 do i = 1,2
  if ( diffdim(i) > 0.5d0*Box(i) ) then
   diffdim(i) = diffdim(i) - Box(i)
  else if ( diffdim(i) < -0.5d0*Box(i) ) then
   diffdim(i) = diffdim(i) + Box(i)
  end if
  if ( dimCM(i) >= Box(i) ) then
   dimCM(i) = dimCM(i) - Box(i)
  else if ( dimCM(i) < 0. ) then
   dimCM(i) = dimCM(i) + Box(i)
  end if
 end do

 ! Find distance between the two dimers and compare to bond length
 Con1diff = (diffdim(1)**2 + diffdim(2)**2 + diffdim(3)**2) - Dcn**2
 ! Find distance between COM and stationary starting point
 Con1diffx = dimCM(1) - motorCM(k,1)
 Con1diffy = dimCM(2) - motorCM(k,2)
 Con1diffz = dimCM(3) - motorCM(k,3)

 !If distance is within constraint, exit and move on; if not, repeat with Lagrange mul
 if ( (abs(Con1diff) < Con1Conv .AND. abs(Con1diffx) < Con1Convx .and. abs(Con1diffy) < Con1Convy .and. &
       abs(Con1diffz) < Con1Convz) .OR. Con1Acc == Con1Max ) exit

  !Langrange multiplier formula for G
  G = - 0.25 * Con1diff*IIMcIMn / (diffdim(1)*diffdimP(1)+diffdim(2)*diffdimP(2)+diffdim(3)*diffdimP(3))
  Hx = - 4.0 * IIMcIMn * Con1diffx
  Hy = - 4.0 * IIMcIMn * Con1diffy
  Hz = - 4.0 * IIMcIMn * Con1diffz

  ! Position update
  Xdim(k,1) = prXdim(1) + QXdim(k,1)*MDt + G*diffdimP(1)/Mc + 0.5*Hx/Mc
  Ydim(k,1) = prYdim(1) + QYdim(k,1)*MDt + G*diffdimP(2)/Mc + 0.5*Hy/Mc
  Zdim(k,1) = prZdim(1) + QZdim(k,1)*MDt + G*diffdimP(3)/Mc + 0.5*Hz/Mc
  Xdim(k,2) = prXdim(2) + QXdim(k,2)*MDt - G*diffdimP(1)/Mn + 0.5*Hx/Mn
  Ydim(k,2) = prYdim(2) + QYdim(k,2)*MDt - G*diffdimP(2)/Mn + 0.5*Hy/Mn
  Zdim(k,2) = prZdim(2) + QZdim(k,2)*MDt - G*diffdimP(3)/Mn + 0.5*Hz/Mn
  do i = 1,2
   Xdim(k,i) = Xdim(k,i) - Box(1)*floor(Xdim(k,i)/Box(1))
   Ydim(k,i) = Ydim(k,i) - Box(2)*floor(Ydim(k,i)/Box(2))
   Zdim(k,i) = Zdim(k,i) - Box(3)*floor(Zdim(k,i)/Box(3))
  end do
  ! Velocity update
  QXdim(k,1) = QXdim(k,1) + G*diffdimP(1)/(Mc*MDt) + 0.5*Hx/(Mc*MDt)
  QYdim(k,1) = QYdim(k,1) + G*diffdimP(2)/(Mc*MDt) + 0.5*Hy/(Mc*MDt)
  QZdim(k,1) = QZdim(k,1) + G*diffdimP(3)/(Mc*MDt) + 0.5*Hz/(Mc*MDt)
  QXdim(k,2) = QXdim(k,2) - G*diffdimP(1)/(Mn*MDt) + 0.5*Hx/(Mn*MDt)
  QYdim(k,2) = QYdim(k,2) - G*diffdimP(2)/(Mn*MDt) + 0.5*Hy/(Mn*MDt)
  QZdim(k,2) = QZdim(k,2) - G*diffdimP(3)/(Mn*MDt) + 0.5*Hz/(Mn*MDt)
  
  Con1Acc = Con1Acc + 1
 end do

end do

!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Force Calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!Set forces and potential to 0 at time t+h
do k = 1,numdims
 do i = 1,2
  FXdim(k,i) = 0.
  FYdim(k,i) = 0.
  FZdim(k,i) = 0.
  fwallz(k,i) = 0.
  wall_pot(k,i) = 0.
  TPotential_part(k,i) = 0.
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(FpairX,FpairY,FpairZ,Potential,i,neighbour,neighbour_num,k) &
!$OMP PRIVATE(j) SCHEDULE(static)

 do j = 1,neighbour_num(k,i)
  FpairX(k,i,j) = 0.
  FpairY(k,i,j) = 0.
  FpairZ(k,i,j) = 0.
  Potential(k,i,j) = 0.
 end do

 !$OMP END PARALLEL DO
 end do
end do
TPotential = 0.

maxv = maxval(max_solvent_dist)
! Check for neighbour list construction
if ( maxv > (0.5*skinrad)**2 .or. time_count == time_penetrate .or. mod(cstep,MPCtstep) == 0 ) then

 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(cstep,NT,MDt,Box,Xsol,Ysol,Zsol) &
 !$OMP SHARED(VXsol,VYsol,VZsol,time_count,outside_particle_num,outside_particle) &
 !$OMP PRIVATE(i,time_aft_coll) SCHEDULE(static)

 ! Evolve outside particles to current time step
 do i = 1,outside_particle_num
  Xsol(outside_particle(i)) = Xsol(outside_particle(i)) + VXsol(outside_particle(i))*MDt*time_count                      
  if ( Xsol(outside_particle(i)) >= Box(1) ) then                          
   Xsol(outside_particle(i)) = Xsol(outside_particle(i)) - Box(1)
  else if ( Xsol(outside_particle(i)) < 0. ) then
   Xsol(outside_particle(i)) = Xsol(outside_particle(i)) + Box(1)
  end if
  Ysol(outside_particle(i)) = Ysol(outside_particle(i)) + VYsol(outside_particle(i))*MDt*time_count                       
  if ( Ysol(outside_particle(i)) >= Box(2) ) then                          
   Ysol(outside_particle(i)) = Ysol(outside_particle(i)) - Box(2)
  else if ( Ysol(outside_particle(i)) < 0. ) then
   Ysol(outside_particle(i)) = Ysol(outside_particle(i)) + Box(2)
  end if
  Zsol(outside_particle(i)) = Zsol(outside_particle(i)) + VZsol(outside_particle(i))*MDt*time_count                       

  ! Bounce back from wall
  if ( Zsol(outside_particle(i)) < 0. ) then
   time_aft_coll = Zsol(outside_particle(i)) / VZsol(outside_particle(i))
   VXsol(outside_particle(i)) = -1. * VXsol(outside_particle(i))
   VYsol(outside_particle(i)) = -1. * VYsol(outside_particle(i))
   VZsol(outside_particle(i)) = -1. * VZsol(outside_particle(i))
   Xsol(outside_particle(i)) = Xsol(outside_particle(i)) + 2. * VXsol(outside_particle(i)) * time_aft_coll
   Ysol(outside_particle(i)) = Ysol(outside_particle(i)) + 2. * VYsol(outside_particle(i)) * time_aft_coll
   Zsol(outside_particle(i)) = Zsol(outside_particle(i)) + 2. * VZsol(outside_particle(i)) * time_aft_coll
   if ( Xsol(outside_particle(i)) >= Box(1) ) then                          
    Xsol(outside_particle(i)) = Xsol(outside_particle(i)) - Box(1)
   else if ( Xsol(outside_particle(i)) < 0. ) then
    Xsol(outside_particle(i)) = Xsol(outside_particle(i)) + Box(1)
   end if
   if ( Ysol(outside_particle(i)) >= Box(2) ) then                          
    Ysol(outside_particle(i)) = Ysol(outside_particle(i)) - Box(2)
   else if ( Ysol(outside_particle(i)) < 0. ) then
    Ysol(outside_particle(i)) = Ysol(outside_particle(i)) + Box(2)
   end if
  else if ( Zsol(outside_particle(i)) > box(3) ) then
   time_aft_coll = (Zsol(outside_particle(i)) - box(3)) / VZsol(outside_particle(i))
   VXsol(outside_particle(i)) = -1. * VXsol(outside_particle(i))
   VYsol(outside_particle(i)) = -1. * VYsol(outside_particle(i))
   VZsol(outside_particle(i)) = -1. * VZsol(outside_particle(i))
   Xsol(outside_particle(i)) = Xsol(outside_particle(i)) + 2. * VXsol(outside_particle(i)) * time_aft_coll
   Ysol(outside_particle(i)) = Ysol(outside_particle(i)) + 2. * VYsol(outside_particle(i)) * time_aft_coll
   Zsol(outside_particle(i)) = Zsol(outside_particle(i)) + 2. * VZsol(outside_particle(i)) * time_aft_coll
   if ( Xsol(outside_particle(i)) >= Box(1) ) then                          
    Xsol(outside_particle(i)) = Xsol(outside_particle(i)) - Box(1)
   else if ( Xsol(outside_particle(i)) < 0. ) then
    Xsol(outside_particle(i)) = Xsol(outside_particle(i)) + Box(1)
   end if
   if ( Ysol(outside_particle(i)) >= Box(2) ) then                          
    Ysol(outside_particle(i)) = Ysol(outside_particle(i)) - Box(2)
   else if ( Ysol(outside_particle(i)) < 0. ) then
    Ysol(outside_particle(i)) = Ysol(outside_particle(i)) + Box(2)
   end if
  end if
 end do
 !$OMP END PARALLEL DO
! Sorting algorithm
if ( mod(cstep,sort_step) == 0 ) then
 call sort_hilbert(vol,NT,box,xsol,ysol,zsol,flag,type_B,vxsol,vysol,vzsol)
end if

 ! Equate total number of neighbours and neighbour list to zero
 neighbour_num = 0
 neighbour = 0
 outside_particle_num = 0
 outside_particle = 0
 inside_particle_num = 0
 inside_particle = 0
 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(numdims,NT,box,sphere,cutoffC,cutoffN,skinrad,xsol,ysol,zsol) &
 !$OMP SHARED(xdim,ydim,zdim,neighbour_num,neighbour) &
 !$OMP PRIVATE(k)
 do k = 1,numdims
  ! Catalytic sphere
  call neighbour_list(NT,Box,cutoffC,skinrad,Xsol,Ysol,Zsol,Xdim(k,1),Ydim(k,1),Zdim(k,1),neighbour_num(k,1),neighbour(k,1,:))
  ! Noncatalytic sphere
  call neighbour_list(NT,Box,cutoffN,skinrad,Xsol,Ysol,Zsol,Xdim(k,2),Ydim(k,2),Zdim(k,2),neighbour_num(k,2),neighbour(k,2,:))
 end do
 !$OMP END PARALLEL DO
 call inout_particle_list(NT,box,cutoffC,cutoffN,skinrad,xsol,ysol,zsol,xdim,ydim,zdim,&
outside_particle_num,outside_particle,inside_particle_num,inside_particle,numdims)
 max_solvent_dist = 0
 time_penetrate = floor((0.5*skinrad / (sqrt(maxval(VXsol*VXsol+VYsol*VYsol+VZsol*VZsol)))) / MDt)
 time_count = 0
end if

!$OMP PARALLEL DO DEFAULT(NONE)&
!$OMP SHARED(sphere,nt,box,flag,type_B,xsol,ysol,zsol,xdim,ydim,zdim,fpairx,fpairy,fpairz) &
!$OMP SHARED(potential,rpab,rpba,cutoffc2,cutoffn2,eac,ebc,rc6,rc12,r12rc12,r6rc6) &
!$OMP SHARED(ean,ebn,rn6,rn12,r12rn12,r6rn6,random,neighbour_num,neighbour,zdimeq) &
!$OMP SHARED(k_harm,fwallz,wall_pot) &
!$OMP PRIVATE(k)

do k = 1,numdims
 ! Solvent with catalytic sphere
 call LJ_force(k,sphere(1),NT,Box,flag,type_B,Xsol,Ysol,Zsol,Xdim(k,1),Ydim(k,1),Zdim(k,1),FpairX(k,1,:),FpairY(k,1,:), &
               FpairZ(k,1,:),Potential(k,1,:),RPab,RPba,cutoffC2,EaC,EbC,Rc6,Rc12,r12Rc12,r6Rc6,random, &
               neighbour_num(k,1),neighbour(k,1,:))

 ! Solvent with noncatalytic sphere
 call LJ_force(k,sphere(2),NT,Box,flag,type_B,Xsol,Ysol,Zsol,Xdim(k,2),Ydim(k,2),Zdim(k,2),FpairX(k,2,:),FpairY(k,2,:), &
               FpairZ(k,2,:),Potential(k,2,:),RPab,RPba,cutoffN2,EaN,EbN,Rn6,Rn12,r12Rn12,r6Rn6,random, &
               neighbour_num(k,2),neighbour(k,2,:))

 ! Harmonic potential for z direction constraint, C sphere
 call harmonic_wall(zdim(k,1),zdimeq,k_harm,fwallz(k,1),wall_pot(k,1))

 ! Harmonic potential for z direction constraint, N sphere
 call harmonic_wall(zdim(k,2),zdimeq,k_harm,fwallz(k,2),wall_pot(k,2))
end do

!$OMP END PARALLEL DO

! Calculate total potential and force on dimer
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numdims,neighbour_num,fxdim,fydim,fzdim,fpairx,fpairy,fpairz) &
!$OMP SHARED(fwallz,wall_pot,Potential,Tpotential_part) &
!$OMP PRIVATE(k,i,j)
do k = 1,numdims
 do i = 1,2
  do j = 1,neighbour_num(k,i)
   Tpotential_part(k,i) = Tpotential_part(k,i) + Potential(k,i,j)
   FXdim(k,i) = FXdim(k,i) - FpairX(k,i,j)                          ! Addition of LJ force on dim., x
   FYdim(k,i) = FYdim(k,i) - FpairY(k,i,j)                          ! Addition of LJ force on dim., y
   FZdim(k,i) = FZdim(k,i) - FpairZ(k,i,j)                          ! Addition of LJ force on dim., z
  end do
  FZdim(k,i) = FZdim(k,i) + fwallz(k,i)
  TPotential_part(k,i) = TPotential_part(k,i) + wall_pot(k,i)
 end do
end do
!$OMP END PARALLEL DO
do k = 1,numdims
 do i = 1,2
  TPotential = TPotential + TPotential_part(k,i)
 end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Calculations after force calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!Velocity of solvent molecules
!Find velocities according to velocity Verlet algorithm

do k = 1,numdims
 do i = 1,2
 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(hIMaMDt,hIMbMDT) &
 !$OMP SHARED(VXsol,VYsol,VZsol,FpairX,FpairY,FpairZ,flag) &
 !$OMP SHARED(neighbour,neighbour_num,i,k) &
 !$OMP PRIVATE(j) SCHEDULE(static)

!do k = 1,numdims
! do i = 1,2
  do j = 1,neighbour_num(k,i)
   if ( flag(neighbour(k,i,j)) == 0 ) then                 ! Velocity half-step update
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMaMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMaMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMaMDt * FpairZ(k,i,j)
   else if ( flag(neighbour(k,i,j)) == 1 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMbMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMbMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMbMDt * FpairZ(k,i,j)
   end if
  end do
!end do
!end do
 !$OMP END PARALLEL DO
end do
end do
!!!!!!Velocity of Dimer(following Andersen, 1983, Appendix C)
!Find value of velocity without Langrange multiplier

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numdims,vxdim,vydim,vzdim,qxdim,qydim,qzdim,fxdim,fydim,fzdim) &
!$OMP SHARED(hIMcMDt,hIMnMDt,box,IIMcIMn,Dcn,Mc,Mn,xdim,ydim,zdim) &
!$OMP PRIVATE(k,Con2Acc,dvxdim,dvydim,dvzdim,diffdim,Con2diff,Con2diffx,Con2diffy) &
!$OMP PRIVATE(Con2diffz,O,Px,Py,Pz)

do k = 1,numdims
 VXdim(k,1) = QXdim(k,1) + hIMcMDt*FXdim(k,1)
 VYdim(k,1) = QYdim(k,1) + hIMcMDt*FYdim(k,1)
 VZdim(k,1) = QZdim(k,1) + hIMcMDt*FZdim(k,1)
 VXdim(k,2) = QXdim(k,2) + hIMnMDt*FXdim(k,2)
 VYdim(k,2) = QYdim(k,2) + hIMnMDt*FYdim(k,2)
 VZdim(k,2) = QZdim(k,2) + hIMnMDt*FZdim(k,2)

 Con2Acc = 0                                             !Begin with zeroeth time

 do
  !Calculate difference in dimer velocities
  dVXdim = VXdim(k,1) - VXdim(k,2)
  dVYdim = VYdim(k,1) - VYdim(k,2)
  dVZdim = VZdim(k,1) - VZdim(k,2)
  !Difference in dimer positions at t+h
  diffdim(1) = Xdim(k,1) - Xdim(k,2)
  diffdim(2) = Ydim(k,1) - Ydim(k,2)
  diffdim(3) = Zdim(k,1) - Zdim(k,2)
  !Apply periodic conditions
  do i = 1,2
   if ( diffdim(i) > 0.5d0*Box(i) ) then
    diffdim(i) = diffdim(i) - Box(i)
   else if ( diffdim(i) < -0.5d0*Box(i) ) then
    diffdim(i) = diffdim(i) + Box(i)
   end if
  end do

  ! Check if dot product=0 is within acceptable range
  Con2diff = dVXdim*diffdim(1) + dVYdim*diffdim(2) + dVZdim*diffdim(3)
  Con2diffx = 0.5*(VXdim(k,1)+VXdim(k,2))
  Con2diffy = 0.5*(VYdim(k,1)+VYdim(k,2))
  Con2diffz = 0.5*(VZdim(k,1)+VZdim(k,2))

  if ( (dabs(Con2diff) < Con2Conv .and. dabs(Con2diffx) < Con2Convx .and. abs(Con2diffy) < Con2Convy &
.and. abs(Con2diffz) < Con2Convz) .or. Con2Acc == Con2Max ) exit

  O = - 0.25 * (IIMcIMn*Con2diff) / (Dcn**2)                     !Calculate Lagrange mult.
  Px = - 4.0 * Con2diffx * IIMcIMn
  Py = - 4.0 * Con2diffy * IIMcIMn
  Pz = - 4.0 * Con2diffz * IIMcIMn

  VXdim(k,1) = VXdim(k,1) + O*diffdim(1)/Mc + 0.5*Px/Mc
  VYdim(k,1) = VYdim(k,1) + O*diffdim(2)/Mc + 0.5*Py/Mc
  VZdim(k,1) = VZdim(k,1) + O*diffdim(3)/Mc + 0.5*Pz/Mc
  VXdim(k,2) = VXdim(k,2) - O*diffdim(1)/Mn + 0.5*Px/Mn
  VYdim(k,2) = VYdim(k,2) - O*diffdim(2)/Mn + 0.5*Py/Mn
  VZdim(k,2) = VZdim(k,2) - O*diffdim(3)/Mn + 0.5*Pz/Mn

  Con2Acc = Con2Acc+1
 end do

end do

!$OMP END PARALLEL DO

deallocate(prXdim,prYdim,prZdim,diffdimP,diffdim)
deallocate(fwallz,wall_pot,Tpotential_part)

end subroutine evolve

end module
