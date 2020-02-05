module evolve_mod
implicit none

contains

subroutine evolve(cstep,time,NT,MDt,MPCtstep,numdims,Ma,Mb,Ms,Mc,Mn,Dcn,Box,cutoff,Rc,Rn,EaC,EbC,EsC,EaN,EbN,EsN,&
                  Rw,Ew,Rm_delta,Em,RPab,RPba,sphere,flag,xsol,Ysol,Zsol,VXsol,VYsol,VZsol,Xdim,&
                  Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,motorCM,TPotential,Potential,FpairX,&
                  FpairY,FpairZ,skinrad,max_solvent_dist,neighbour_num,neighbour,time_count,&
                  time_penetrate,outside_particle_num,outside_particle,inside_particle_num,inside_particle,&
                  vol,sort_step,react_box,react_dir)
use sort_hilbert_mod
use neighbour_list_mod
use reactive_grid_mod
use forces_mod
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=sp) :: cstep                                       ! Current time step
integer(kind=dp) :: time                                        ! Time loop variable
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
real(kind=dp) , intent(in) :: MDt                               ! MD time step
integer(kind=sp) :: MPCtstep                                    ! MPC step happens every nth step
integer , intent(in) :: numdims
real(kind=dp) , intent(in) :: Ma                                ! Mass of one A particle
real(kind=dp) , intent(in) :: Mb                                ! Mass of one B particle
real(kind=dp) , intent(in) :: Ms                                ! Mass of one S particle
real(kind=dp) , intent(in) :: Mc                                ! Mass of C sphere
real(kind=dp) , intent(in) :: Mn                                ! Mass of N sphere
real(kind=dp) , intent(in) :: Dcn                               ! Internuclear bond length
real(kind=dp) , dimension(:), intent(in) :: Box                 ! Length of simulation box
real(kind=dp) , intent(in) :: cutoff                            ! Cut. rad. C sphere LJ pot., squared
real(kind=dp) , intent(in) :: Rc                                ! Radius C sphere
real(kind=dp) , intent(in) :: Rn                                ! Radius N sphere
real(kind=dp) , intent(in) :: EaC                               ! Reduced LJ Potential of A with C
real(kind=dp) , intent(in) :: EbC                               ! Reduced LJ Potential of B with C
real(kind=dp) , intent(in) :: EsC                               ! Reduced LJ Potential of S with C
real(kind=dp) , intent(in) :: EaN                               ! Reduced LJ Potential of A with N
real(kind=dp) , intent(in) :: EbN                               ! Reduced LJ Potential of B with N
real(kind=dp) , intent(in) :: EsN                               ! Reduced LJ Potential of S with N
real(kind=dp) , intent(in) :: Rw                                ! Cutoff distance for wall
real(kind=dp) , intent(in) :: Ew                                ! Reduced LJ Potential for wall
real(kind=dp) , intent(in) :: Rm_delta
real(kind=dp) , intent(in) :: Em                                ! Reduced LJ Potential for motor-motor
real(kind=dp) , intent(in) :: RPab                              ! Probability of reaction from F to G
real(kind=dp) , intent(in) :: RPba                              ! Probability of reaction from G to F
integer , dimension(:) , intent(in) :: sphere                   ! Sphere identity; C=0;N=1
integer(kind=sp) , dimension(:), intent(inout) :: flag          ! Identity of solvent
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
real(kind=dp) , dimension(:,:) , intent(in) :: motorCM
real(kind=dp) , intent(out) :: TPotential                       ! Total potential
real(kind=dp) , dimension(:,:,:) , intent(inout) :: Potential   ! Potential per particle
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairX      ! LJ force in x direction
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairY      ! LJ force in y direction
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairZ      ! LJ force in z direction
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
integer , dimension(:) , intent(in) :: react_box
integer , intent(in) :: react_dir

! Definitions for internal variables
integer(kind=dp) :: i,j,k,l                                     ! Loop variable
real(kind=dp) :: hIMaMDt                                        ! Value used a lot for calculations
real(kind=dp) :: hIMbMDt                                        ! Value used a lot for calculations
real(kind=dp) :: hIMsMDt                                        ! Value used a lot for calculations
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
real(kind=dp) :: G                                          
real(kind=dp) :: O                                          
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
real(kind=dp) :: Rcc6
real(kind=dp) :: Rcc12
real(kind=dp) :: Rcn6
real(kind=dp) :: Rcn12
real(kind=dp) :: Rnn6
real(kind=dp) :: Rnn12
real(kind=dp) :: r6Rcc6
real(kind=dp) :: r12Rcc12
real(kind=dp) :: r6Rcn6
real(kind=dp) :: r12Rcn12
real(kind=dp) :: r6Rnn6
real(kind=dp) :: r12Rnn12
real(kind=dp) :: mot_cutoff_CC2
real(kind=dp) :: mot_cutoff_CN2
real(kind=dp) :: mot_cutoff_NN2
real(kind=dp) :: Rwc3
real(kind=dp) :: Rwc9
real(kind=dp) :: r3Rwc3
real(kind=dp) :: r9Rwc9
real(kind=dp) :: Rwn3
real(kind=dp) :: Rwn9
real(kind=dp) :: r3Rwn3
real(kind=dp) :: r9Rwn9
real(kind=dp) :: wall_cutoffC2
real(kind=dp) :: wall_cutoffN2
real(kind=dp) , dimension(:,:,:,:) , allocatable :: FmotX       ! Mot-mot LJ force, x 
real(kind=dp) , dimension(:,:,:,:) , allocatable :: FmotY       ! Mot-mot LJ force, y 
real(kind=dp) , dimension(:,:,:,:) , allocatable :: FmotZ       ! Mot-mot LJ force, z 
real(kind=dp) , dimension(:,:,:,:) , allocatable :: Fmot_pot    ! Mot-mot LJ potential 
real(kind=dp) , dimension(:,:) , allocatable :: fwallz          ! Force of wall on monomer, z
real(kind=dp) , dimension(:,:) , allocatable :: wall_pot        ! Potential of wall with monomer
real(kind=dp) , dimension(:,:) , allocatable :: TPotential_part
real(kind=dp) :: maxv
real(kind=dp) :: maxt
real(kind=dp) :: vel

allocate(prXdim(2),prYdim(2),prZdim(2),diffdimP(3),diffdim(3))
allocate(fwallz(numdims,2),wall_pot(numdims,2),TPotential_part(numdims,2))
allocate(FmotX(numdims,2,numdims,2),FmotY(numdims,2,numdims,2),FmotZ(numdims,2,numdims,2),&
Fmot_pot(numdims,2,numdims,2))

!!!!!! Define variables within subroutine
hIMaMDt = (0.5d0*MDt)/Ma                                        ! Value used a lot for calculations
hIMbMDt = (0.5d0*MDt)/Mb                                        ! Value used a lot for calculations
hIMsMDt = (0.5d0*MDt)/Ms                                        ! Value used a lot for calculations
hIMcMDt = (0.5d0*MDt)/Mc                                        ! Value used a lot for calculations
hIMnMDt = (0.5d0*MDt)/Mn                                        ! Value used a lot for calculations
IIMcIMn = 1.d0/((1.d0/Mc)+(1.d0/Mn))                            ! Value used for first Lagr. mult.
! Motor-solvent interaction parameters
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
! Motor-motor
Rcc6 = (Rc+Rc+Rm_delta)**6
Rcc12 = (Rc+Rc+Rm_delta)**12
Rcn6 = (Rc+Rn+Rm_delta)**6
Rcn12 = (Rc+Rn+Rm_delta)**12
Rnn6 = (Rn+Rn+Rm_delta)**6
Rnn12 = (Rn+Rn+Rm_delta)**12
r6Rcc6 = 6.d0*Rcc6
r12Rcc12 = 12.d0*Rcc12
r6Rcn6 = 6.d0*Rcn6
r12Rcn12 = 12.d0*Rcn12
r6Rnn6 = 6.d0*Rnn6
r12Rnn12 = 12.d0*Rnn12
mot_cutoff_CC2 = ((Rc+Rc+Rm_delta)*cutoff)**2
mot_cutoff_CN2 = ((Rc+Rn+Rm_delta)*cutoff)**2
mot_cutoff_NN2 = ((Rn+Rn+Rm_delta)*cutoff)**2
! Motor-wall interaction parameters
Rwc3 = (Rw+Rc)**3
Rwc9 = (Rw+Rc)**9
r3Rwc3 = 3.d0*Rwc3
r9Rwc9 = 9.d0*Rwc9
Rwn3 = (Rw+Rn)**3
Rwn9 = (Rw+Rn)**9
r3Rwn3 = 3.d0*Rwn3
r9Rwn9 = 9.d0*Rwn9
wall_cutoffC2 = ((Rw+Rc)*3.d0**(1./6.))**2
wall_cutoffN2 = ((Rw+Rn)*3.d0**(1./6.))**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Calculations before force calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Motion of the solvent molecules
do k = 1,numdims
 do i = 1,2
 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(hIMaMDt,hIMbMDt,hIMsMDt) &
 !$OMP SHARED(VXsol,VYsol,VZsol,FpairX,FpairY,FpairZ,flag) &
 !$OMP SHARED(neighbour,neighbour_num,i,k) &
 !$OMP PRIVATE(j) SCHEDULE(static)

  do j = 1,neighbour_num(k,i)
   if ( flag(neighbour(k,i,j)) == 0 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMaMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMaMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMaMDt * FpairZ(k,i,j)
   else if ( flag(neighbour(k,i,j)) == 1 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMbMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMbMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMbMDt * FpairZ(k,i,j)
   else if ( flag(neighbour(k,i,j)) == 2 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMsMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMsMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMsMDt * FpairZ(k,i,j)
   end if
  end do
!$OMP END PARALLEL DO
end do
end do

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(cstep,NT,MDt,Box,Xsol,Ysol,Zsol,flag,max_solvent_dist,react_box,react_dir) &
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
   call reactive_grid(react_box,react_dir,box,xsol(inside_particle(i)),ysol(inside_particle(i)),&
            vxsol(inside_particle(i)),vysol(inside_particle(i)),flag(inside_particle(i)),time_aft_coll)
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
!$OMP PRIVATE(k,i,prXdim,prYdim,prZdim,diffdimP,Con1Acc,diffdim,dimCM,Con1diff,G)

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
 end do

 ! Find distance between the two dimers and compare to bond length
 Con1diff = (diffdim(1)**2 + diffdim(2)**2 + diffdim(3)**2) - Dcn**2

 !If distance is within constraint, exit and move on; if not, repeat with Lagrange mul
 if ( (abs(Con1diff)) < Con1Conv .OR. Con1Acc == Con1Max ) exit

  !Langrange multiplier formula for G
  G = - 0.25 * Con1diff*IIMcIMn / (diffdim(1)*diffdimP(1)+diffdim(2)*diffdimP(2)+diffdim(3)*diffdimP(3))

  ! Position update
  Xdim(k,1) = prXdim(1) + QXdim(k,1)*MDt + 2.0*G*diffdimP(1)/Mc
  Ydim(k,1) = prYdim(1) + QYdim(k,1)*MDt + 2.0*G*diffdimP(2)/Mc
  Zdim(k,1) = prZdim(1) + QZdim(k,1)*MDt + 2.0*G*diffdimP(3)/Mc
  Xdim(k,2) = prXdim(2) + QXdim(k,2)*MDt - 2.0*G*diffdimP(1)/Mn
  Ydim(k,2) = prYdim(2) + QYdim(k,2)*MDt - 2.0*G*diffdimP(2)/Mn
  Zdim(k,2) = prZdim(2) + QZdim(k,2)*MDt - 2.0*G*diffdimP(3)/Mn
  do i = 1,2
   Xdim(k,i) = Xdim(k,i) - Box(1)*floor(Xdim(k,i)/Box(1))
   Ydim(k,i) = Ydim(k,i) - Box(2)*floor(Ydim(k,i)/Box(2))
  end do
  ! Velocity update
  QXdim(k,1) = QXdim(k,1) + 2.0*G*diffdimP(1)/(Mc*MDt)
  QYdim(k,1) = QYdim(k,1) + 2.0*G*diffdimP(2)/(Mc*MDt)
  QZdim(k,1) = QZdim(k,1) + 2.0*G*diffdimP(3)/(Mc*MDt)
  QXdim(k,2) = QXdim(k,2) - 2.0*G*diffdimP(1)/(Mn*MDt)
  QYdim(k,2) = QYdim(k,2) - 2.0*G*diffdimP(2)/(Mn*MDt)
  QZdim(k,2) = QZdim(k,2) - 2.0*G*diffdimP(3)/(Mn*MDt)
  
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
  do j = 1,numdims
   do l = 1,2
    FmotX(k,i,j,l) = 0.
    FmotY(k,i,j,l) = 0.
    FmotZ(k,i,j,l) = 0.
    Fmot_pot(k,i,j,l) = 0.
   end do
  end do

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

maxv = max_solvent_dist(1)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP REDUCTION(MAX:maxv) SHARED(max_solvent_dist,NT) PRIVATE(i)
do i = 1,NT
 maxv = max(maxv,max_solvent_dist(i))
end do
!$OMP END PARALLEL DO
! Check for neighbour list construction
if ( maxv > (0.5*skinrad)**2 .or. time_count == time_penetrate .or. mod(cstep,MPCtstep) == 0 ) then

 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(cstep,NT,MDt,Box,Xsol,Ysol,flag,Zsol,react_box,react_dir) &
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
   call reactive_grid(react_box,react_dir,box,xsol(outside_particle(i)),ysol(outside_particle(i)),&
            vxsol(outside_particle(i)),vysol(outside_particle(i)),flag(outside_particle(i)),time_aft_coll)
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
 call sort_hilbert(vol,NT,box,xsol,ysol,zsol,flag,vxsol,vysol,vzsol)
end if

 ! Equate total number of neighbours and neighbour list to zero
 neighbour_num = 0
 outside_particle_num = 0
 inside_particle_num = 0
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numdims,NT,box,sphere,cutoffC,cutoffN,skinrad,xsol,ysol,zsol) &
!$OMP SHARED(xdim,ydim,zdim,neighbour_num,neighbour) &
!$OMP PRIVATE(k)
 do k = 1,numdims
  ! Catalytic sphere
  call neighbour_list(NT,Box,cutoffC,skinrad,Xsol,Ysol,Zsol,Xdim(k,1),Ydim(k,1),Zdim(k,1),neighbour_num(k,1),&
neighbour(k,1,:))
  ! Noncatalytic sphere
  call neighbour_list(NT,Box,cutoffN,skinrad,Xsol,Ysol,Zsol,Xdim(k,2),Ydim(k,2),Zdim(k,2),neighbour_num(k,2),&
neighbour(k,2,:))
 end do
!$OMP END PARALLEL DO
 call inout_particle_list(NT,box,cutoffC,cutoffN,skinrad,xsol,ysol,zsol,xdim,ydim,zdim,&
outside_particle_num,outside_particle,inside_particle_num,inside_particle,numdims)
 max_solvent_dist = 0
 maxt = vxsol(1)**2+vysol(1)**2+vzsol(1)**2
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP REDUCTION(MAX:maxt) SHARED(NT,vxsol,vysol,vzsol) PRIVATE(vel)
do i = 1,NT
 vel = vxsol(i)**2+vysol(i)**2+vzsol(i)**2
 maxt = max(maxt,vel)
end do
!$OMP END PARALLEL DO
 time_penetrate = floor((0.5*skinrad / (sqrt(maxt))) / MDt)
 time_count = 0
end if

do k = 1,numdims
 ! Solvent with catalytic sphere
 call LJ_mot_sol(sphere(1),NT,Box,flag,Xsol,Ysol,Zsol,Xdim(k,1),Ydim(k,1),Zdim(k,1),FpairX(k,1,:),FpairY(k,1,:), &
               FpairZ(k,1,:),Potential(k,1,:),RPab,RPba,cutoffC2,EaC,EbC,EsC,Rc6,Rc12,r12Rc12,r6Rc6, &
               neighbour_num(k,1),neighbour(k,1,:))

 ! Solvent with noncatalytic sphere
 call LJ_mot_sol(sphere(2),NT,Box,flag,Xsol,Ysol,Zsol,Xdim(k,2),Ydim(k,2),Zdim(k,2),FpairX(k,2,:),FpairY(k,2,:), &
               FpairZ(k,2,:),Potential(k,2,:),RPab,RPba,cutoffN2,EaN,EbN,EsN,Rn6,Rn12,r12Rn12,r6Rn6, &
               neighbour_num(k,2),neighbour(k,2,:))

 ! Motor-motor potential interactions
 do j = 1,numdims
  if ( j /= k ) then
   ! C sphere with other C spheres
   call LJ_mot_mot(box,Xdim(k,1),Ydim(k,1),Zdim(k,1),Xdim(j,1),Ydim(j,1),Zdim(j,1),FmotX(k,1,j,1),&
                   FmotY(k,1,j,1),FmotZ(k,1,j,1),Fmot_pot(k,1,j,1),mot_cutoff_CC2,Em,Rcc6,Rcc12,r6Rcc6,r12Rcc12)
   ! C sphere with other N spheres
   call LJ_mot_mot(box,Xdim(k,1),Ydim(k,1),Zdim(k,1),Xdim(j,2),Ydim(j,2),Zdim(j,2),FmotX(k,1,j,2),&
                   FmotY(k,1,j,2),FmotZ(k,1,j,2),Fmot_pot(k,1,j,2),mot_cutoff_CN2,Em,Rcn6,Rcn12,r6Rcn6,r12Rcn12)
   ! N sphere with other C spheres
   call LJ_mot_mot(box,Xdim(k,2),Ydim(k,2),Zdim(k,2),Xdim(j,1),Ydim(j,1),Zdim(j,1),FmotX(k,2,j,1),&
                   FmotY(k,2,j,1),FmotZ(k,2,j,1),Fmot_pot(k,2,j,1),mot_cutoff_CN2,Em,Rcn6,Rcn12,r6Rcn6,r12Rcn12)
   ! N sphere with other N spheres
   call LJ_mot_mot(box,Xdim(k,2),Ydim(k,2),Zdim(k,2),Xdim(j,2),Ydim(j,2),Zdim(j,2),FmotX(k,2,j,2),&
                   FmotY(k,2,j,2),FmotZ(k,2,j,2),Fmot_pot(k,2,j,2),mot_cutoff_NN2,Em,Rnn6,Rnn12,r6Rnn6,r12Rnn12)
  end if 
 end do

 ! Catalytic sphere interaction with wall in z-direction
 call LJ_mot_wall(box,zdim(k,1),fwallz(k,1),Ew,Rwc3,Rwc9,r3Rwc3,r9Rwc9,wall_cutoffC2,wall_pot(k,1))

 ! Noncatalytic sphere interaction with wall in z-direction
 call LJ_mot_wall(box,zdim(k,2),fwallz(k,2),Ew,Rwn3,Rwn9,r3Rwn3,r9Rwn9,wall_cutoffN2,wall_pot(k,2))

end do

! Calculate total potential and force on dimer
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numdims,neighbour_num,fxdim,fydim,fzdim,fpairx,fpairy,fpairz) &
!$OMP SHARED(fwallz,wall_pot,Potential,Tpotential_part,fmotx,fmoty,fmotz,fmot_pot) &
!$OMP PRIVATE(k,i,j,l)
do k = 1,numdims
 do i = 1,2
  do j = 1,neighbour_num(k,i)
   Tpotential_part(k,i) = Tpotential_part(k,i) + Potential(k,i,j)
   FXdim(k,i) = FXdim(k,i) - FpairX(k,i,j)                          ! Addition of LJ force on dim., x
   FYdim(k,i) = FYdim(k,i) - FpairY(k,i,j)                          ! Addition of LJ force on dim., y
   FZdim(k,i) = FZdim(k,i) - FpairZ(k,i,j)                          ! Addition of LJ force on dim., z
  end do
  do j = 1,numdims
   do l = 1,2
    FXdim(k,i) = FXdim(k,i) + FmotX(k,i,j,l)
    FYdim(k,i) = FYdim(k,i) + FmotY(k,i,j,l)
    FZdim(k,i) = FZdim(k,i) + FmotZ(k,i,j,l)
    if ( j < k ) then
     Tpotential_part(k,i) = TPotential_part(k,i) + Fmot_pot(k,i,j,l)
    end if
   end do
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
 !$OMP SHARED(hIMaMDt,hIMbMDt,hIMsMDt) &
 !$OMP SHARED(VXsol,VYsol,VZsol,FpairX,FpairY,FpairZ,flag) &
 !$OMP SHARED(neighbour,neighbour_num,i,k) &
 !$OMP PRIVATE(j) SCHEDULE(static)

  do j = 1,neighbour_num(k,i)
   if ( flag(neighbour(k,i,j)) == 0 ) then                 ! Velocity half-step update
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMaMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMaMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMaMDt * FpairZ(k,i,j)
   else if ( flag(neighbour(k,i,j)) == 1 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMbMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMbMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMbMDt * FpairZ(k,i,j)
   else if ( flag(neighbour(k,i,j)) == 2 ) then
    VXsol(neighbour(k,i,j)) = VXsol(neighbour(k,i,j)) + hIMsMDt * FpairX(k,i,j)
    VYsol(neighbour(k,i,j)) = VYsol(neighbour(k,i,j)) + hIMsMDt * FpairY(k,i,j)
    VZsol(neighbour(k,i,j)) = VZsol(neighbour(k,i,j)) + hIMsMDt * FpairZ(k,i,j)
   end if
  end do
 !$OMP END PARALLEL DO
end do
end do
!!!!!!Velocity of Dimer(following Andersen, 1983, Appendix C)
!Find value of velocity without Langrange multiplier

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(numdims,vxdim,vydim,vzdim,qxdim,qydim,qzdim,fxdim,fydim,fzdim) &
!$OMP SHARED(hIMcMDt,hIMnMDt,box,IIMcIMn,Dcn,Mc,Mn,xdim,ydim,zdim) &
!$OMP PRIVATE(k,Con2Acc,dvxdim,dvydim,dvzdim,diffdim,Con2diff,O)

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
  Con2diff = 2.0 * (dVXdim*diffdim(1) + dVYdim*diffdim(2) + dVZdim*diffdim(3) )

  if ( (dabs(Con2diff)) < Con2Conv .or. Con2Acc == Con2Max ) exit

  O = - 0.25 * (IIMcIMn*Con2diff) / (Dcn**2)                     !Calculate Lagrange mult.

  VXdim(k,1) = VXdim(k,1) + 2.0*O*diffdim(1)/Mc
  VYdim(k,1) = VYdim(k,1) + 2.0*O*diffdim(2)/Mc
  VZdim(k,1) = VZdim(k,1) + 2.0*O*diffdim(3)/Mc
  VXdim(k,2) = VXdim(k,2) - 2.0*O*diffdim(1)/Mn
  VYdim(k,2) = VYdim(k,2) - 2.0*O*diffdim(2)/Mn
  VZdim(k,2) = VZdim(k,2) - 2.0*O*diffdim(3)/Mn

  Con2Acc = Con2Acc+1
 end do

end do

!$OMP END PARALLEL DO

deallocate(prXdim,prYdim,prZdim,diffdimP,diffdim)
deallocate(fwallz,wall_pot,Tpotential_part)
deallocate(FmotX,FmotY,FmotZ,Fmot_pot)

end subroutine evolve

end module
