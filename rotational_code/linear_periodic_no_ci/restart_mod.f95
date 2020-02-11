module restart_mod

contains

        subroutine restart(cstep,iTotA,iTotB,NT,MDt,Box,vol,flag,type_B,sphere,numdims,cutoff,Rc,Rn,EaC,EbC,&
EaN,EbN,zdimeq,k_harm,RPab,RPba,Xdim,Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,Xsol,&
Ysol,Zsol,VXsol,VYsol,VZsol,TPotential,Potential,FpairX,FpairY,FpairZ,motorCM,random,neighbour_num,&
neighbour,max_neighbour,skinrad,max_solvent_dist,time_count,time_penetrate,outside_particle_num,&
outside_particle,inside_particle_num,inside_particle,sort_step,Freq1,Freq4,data_id,sub_group_id,&
sys_dspace_id,sys_crp_list,sys_dset_id,conc_dspace_id,conc_crp_list,conc_dset_id,solv_dspace_id,&
solv_crp_list,solv_dset_id,dimer_dspace_id,dimer_crp_list,dimer_dset_id)
use hdf5
use neighbour_list_mod
use forces_mod
use sort_hilbert_mod
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=sp) , intent(inout) :: cstep                       ! Current timestep
integer(kind=dp) , intent(in) :: iTotA                          ! Total number of A
integer(kind=dp) , intent(in) :: iTotB                          ! Total number of B
integer(kind=sp) , intent(out) :: NT                            ! Total number of A,B particles
integer(kind=dp) :: Vol
real(kind=dp) :: MDt                                            ! MD time step
real(kind=dp) , dimension(:), intent(in) :: Box                 ! Length of simulation box
integer(kind=sp) , intent(in) :: Freq1                          ! Frequency of solvent data written
integer(kind=sp) , intent(in) :: Freq4                          ! Frequency of dimer data written
integer(kind=sp) , dimension(:), intent(out) :: flag            ! Identity of solvent
integer(kind=sp) , dimension(:), intent(out) :: type_B
integer , dimension(:) , intent(in) :: sphere                   ! Sphere identity; C=0;N=1
integer , intent(in) :: numdims                                 ! Number of dimers
real(kind=dp) , intent(in) :: cutoff                            ! Cut. rad. C sphere LJ pot., squared
real(kind=dp) , intent(in) :: Rc                                ! Radius C sphere
real(kind=dp) , intent(in) :: Rn                                ! Radius N sphere
real(kind=dp) , intent(in) :: EaC                               ! Reduced LJ Potential of A with C
real(kind=dp) , intent(in) :: EbC                               ! Reduced LJ Potential of B with C
real(kind=dp) , intent(in) :: EaN                               ! Reduced LJ Potential of A with N
real(kind=dp) , intent(in) :: EbN                               ! Reduced LJ Potential of B with N
real(kind=dp) , intent(in) :: zdimeq                            ! Equil. pos. of mon., z dir.
real(kind=dp) , intent(in) :: k_harm                            ! Force constant
real(kind=dp) , intent(in) :: RPab                              ! Probability of reaction from A to B
real(kind=dp) , intent(in) :: RPba                              ! Probability of reaction from B to A
real(kind=dp) , dimension(:,:) , intent(out) :: Xdim            ! Dimer spheres positions, x
real(kind=dp) , dimension(:,:) , intent(out) :: Ydim            ! Dimer spheres positions, y
real(kind=dp) , dimension(:,:) , intent(out) :: Zdim            ! Dimer spheres positions, z
real(kind=dp) , dimension(:,:) , intent(out) :: VXdim           ! Dimer spheres velocities, x
real(kind=dp) , dimension(:,:) , intent(out) :: VYdim           ! Dimer spheres velocities, y
real(kind=dp) , dimension(:,:) , intent(out) :: VZdim           ! Dimer spheres velocities, z
real(kind=dp) , dimension(:,:) , intent(out) :: FXdim           ! Dimer spheres positions, x
real(kind=dp) , dimension(:,:) , intent(out) :: FYdim           ! Dimer spheres positions, y
real(kind=dp) , dimension(:,:) , intent(out) :: FZdim           ! Dimer spheres positions, z
real(kind=dp) , dimension(:) , intent(out) :: Xsol              ! Solvent positions, x
real(kind=dp) , dimension(:) , intent(out) :: Ysol              ! Solvent positions, y
real(kind=dp) , dimension(:) , intent(out) :: Zsol              ! Solvent positions, z
real(kind=dp) , dimension(:) , intent(out) :: VXsol             ! Solvent positions, x
real(kind=dp) , dimension(:) , intent(out) :: VYsol             ! Solvent positions, y
real(kind=dp) , dimension(:) , intent(out) :: VZsol             ! Solvent positions, z
real(kind=dp) , intent(out) :: TPotential                       ! Total potential
real(kind=dp) , dimension(:,:,:) , intent(inout) :: Potential   ! Potential per particle
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairX      ! LJ force in x direction
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairY      ! LJ force in y direction
real(kind=dp) , dimension(:,:,:) , intent(inout) :: FpairZ      ! LJ force in z direction
real(kind=dp) , dimension(numdims,3) , intent(inout) :: motorCM
real(kind=dp) , dimension(:) , intent(inout) :: random          ! Vector of random numbers
integer(kind=dp) , dimension(:,:) , intent(inout) :: neighbour_num! Number of solvent particles in list
integer(kind=dp) , dimension(:,:,:) , intent(inout) :: neighbour! Identify solvent number in list
integer(kind=dp) , intent(in) :: max_neighbour
real(kind=dp) , intent(in) :: skinrad                           ! Skin radius for neighbour list
real(kind=dp) , dimension(:) , intent(inout) :: max_solvent_dist! Max distance of solvent
integer , intent(inout) :: time_count
integer , intent(inout) :: time_penetrate
integer(kind=dp) , intent(inout) :: outside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: outside_particle
integer :: sort_step
integer(kind=dp) , intent(inout) :: inside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: inside_particle

! Definitions for internal variables
integer(kind=dp) :: rcstep
integer(kind=dp) :: i,j
integer :: k
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
real(kind=dp) :: cutoffN2                                       ! Cut. N sphere LJ pot. squared9
real(kind=dp) , dimension(numdims,2) :: fwallz                  ! Force of wall on monomer, z
real(kind=dp) , dimension(numdims,2) :: wall_pot                ! Potential of wall with monomer

! HDF5 file writing
! Overall File
character(len=7) , parameter :: datafile = "data.h5"
integer(hid_t) , intent(inout) :: data_id
integer(kind=sp) :: err
integer(kind=sp) , dimension(3) , parameter :: rank = (/1,2,3/)
integer(hsize_t) , dimension(1) :: maxdims_1
integer(hsize_t) , dimension(2) :: maxdims_2
integer(hsize_t) , dimension(3) :: maxdims_3
integer(hsize_t) , dimension(1) :: sys_chunk_1
integer(hsize_t) , dimension(1) :: conc_chunk_1
integer(hsize_t) , dimension(2) :: conc1d_chunk_2
integer(hsize_t) , dimension(3) :: conc2d_chunk_3
integer(hsize_t) , dimension(2) :: solv_chunk_2
integer(hsize_t) , dimension(3) :: dimer_chunk_1
integer(hsize_t) , dimension(3) :: dimer_chunk_2
! Group initiation variables
integer(hid_t) , dimension(4) , intent(inout) :: sub_group_id
character(len=13) , dimension(4) , parameter :: sub_group_names = (/"Solvent      ",&
"Concentration","System       ","Dimer        "/)
! System group
integer(hid_t) , dimension(6) , intent(inout) :: sys_dspace_id
integer(hsize_t) , dimension(1) , parameter :: system_init_dims = (/1/)
integer(hid_t) , dimension(6) , intent(inout) :: sys_crp_list
character(len=19) , dimension(6) , parameter :: sys_dset_names = (/"pxsol              ",&
"pysol              ","pzsol              ","total_energy       ","solvent_temperature",&
"dimer_temperature  "/)
integer(hid_t) , dimension(6) , intent(inout) :: sys_dset_id
! Concentration group
integer(hid_t) , dimension(8) , intent(inout):: conc_dspace_id
integer(hsize_t) , dimension(1) , parameter :: conc_init_dims = (/1/)
integer(hid_t) , dimension(8) , intent(inout) :: conc_crp_list
character(len=18) , dimension(8) , parameter :: conc_dset_names = (/"total_A           ",&
"total_B           ","concentration_A   ","concentration_B   ",&
"1D_concentration_A","1D_concentration_B",&
"2D_concentration_A","2D_concentration_B"/)
integer(hid_t) , dimension(8) , intent(inout) :: conc_dset_id
! Solvent group
integer(hid_t) , dimension(8) , intent(inout) :: solv_dspace_id
integer(hid_t) , dimension(8) , intent(inout) :: solv_crp_list
character(len=10) , dimension(8) , parameter :: solv_dset_names = (/"position_x",&
        "position_y","position_z","flag      ","type_B    ","velocity_x","velocity_y","velocity_z"/)
integer(hid_t) , dimension(8) , intent(inout) :: solv_dset_id
! Dimer group
integer(hid_t) , dimension(7) , intent(inout) :: dimer_dspace_id
integer(hid_t) , dimension(7) , intent(inout) :: dimer_crp_list
character(len=10) , dimension(7) , parameter :: dimer_dset_names = (/"position_x",&
"position_y","position_z","velocity_x","velocity_y","velocity_z","nucvec    "/)
integer(hid_t) , dimension(7) , intent(inout) :: dimer_dset_id
integer(hsize_t) , dimension(2) :: solvent_extend_offset
integer(hsize_t) , dimension(2) :: solvent_extend_count
integer(hsize_t) , dimension(3) :: dimer_dimer_extend_offset
integer(hsize_t) , dimension(3) :: dimer_dimer_extend_count
integer(hsize_t) , dimension(1) :: dims_1
integer(hsize_t) , dimension(2) :: dims_2
integer(hsize_t) , dimension(3) :: dims_3
integer(kind=sp) :: layout
integer(hsize_t) , dimension(1) :: solvent_read_dims
integer(hsize_t) , dimension(2) :: dimer_read_dims
integer(hid_t) :: solvent_memspace
integer(hid_t) :: dimer_memspace

!!!!!! Define variables within subroutine
rcstep = cstep / Freq1 + 1
NT = iTotA+iTotB
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
sys_chunk_1 = (/75000/)
conc_chunk_1 = (/75000/)
conc1d_chunk_2 = (/120,800/)
conc2d_chunk_3 = (/120,120,7/)
solv_chunk_2 = (/100000,1/)
dimer_chunk_1 = (/numdims,3,floor(100000./(numdims*3))/)
dimer_chunk_2 = (/numdims,2,floor(100000./(numdims*2))/)

call h5open_f(err)

! Open file
call h5fopen_f(datafile,h5f_acc_rdwr_f,data_id,err)

! Open groups
do i = 1,4,1
 call h5gopen_f(data_id,sub_group_names(i),sub_group_id(i),err)
end do

! System group
do i = 1,6,1
 call h5dopen_f(sub_group_id(3),sys_dset_names(i),sys_dset_id(i),err)
 call h5dget_space_f(sys_dset_id(i),sys_dspace_id(i),err)
 call h5sget_simple_extent_dims_f(sys_dspace_id(i),dims_1,maxdims_1,err)
 call h5dget_create_plist_f(sys_dset_id(i),sys_crp_list(i),err)
 call h5pget_layout_f(sys_crp_list(i),layout,err)
 if ( h5d_chunked_f == layout ) then
  call h5pget_chunk_f(sys_crp_list(i),rank(1),sys_chunk_1,err)
 end if
 call h5screate_simple_f(rank(1),dims_1,sys_dspace_id(i),err,maxdims_1)
end do

! Concentration group
do i = 1,8,1
 call h5dopen_f(sub_group_id(2),conc_dset_names(i),conc_dset_id(i),err)
 call h5dget_space_f(conc_dset_id(i),conc_dspace_id(i),err)
 if ( i >= 1 .and. i <= 4 ) then
  call h5sget_simple_extent_dims_f(conc_dspace_id(i),dims_1,maxdims_1,err)
  call h5dget_create_plist_f(conc_dset_id(i),conc_crp_list(i),err)
  call h5pget_layout_f(conc_crp_list(i),layout,err)
  if ( h5d_chunked_f == layout ) then
   call h5pget_chunk_f(conc_crp_list(i),rank(1),conc_chunk_1,err)
  end if
  call h5screate_simple_f(rank(1),dims_1,conc_dspace_id(i),err,maxdims_1)
 else if ( i >= 5 .and. i <= 6 ) then
  call h5sget_simple_extent_dims_f(conc_dspace_id(i),dims_2,maxdims_2,err)
  call h5dget_create_plist_f(conc_dset_id(i),conc_crp_list(i),err)
  call h5pget_layout_f(conc_crp_list(i),layout,err)
  if ( h5d_chunked_f == layout ) then
   call h5pget_chunk_f(conc_crp_list(i),rank(2),conc1d_chunk_2,err)
  end if
  call h5screate_simple_f(rank(2),dims_2,conc_dspace_id(i),err,maxdims_2)
 else if ( i >= 7 .and. i <= 8 ) then
  call h5sget_simple_extent_dims_f(conc_dspace_id(i),dims_3,maxdims_3,err)
  call h5dget_create_plist_f(conc_dset_id(i),conc_crp_list(i),err)
  call h5pget_layout_f(conc_crp_list(i),layout,err)
  if ( h5d_chunked_f == layout ) then
   call h5pget_chunk_f(conc_crp_list(i),rank(3),conc2d_chunk_3,err)
  end if
  call h5screate_simple_f(rank(3),dims_3,conc_dspace_id(i),err,maxdims_3)
 end if
end do

! Solvent group
do i = 1,8,1
 call h5dopen_f(sub_group_id(1),solv_dset_names(i),solv_dset_id(i),err)
 call h5dget_space_f(solv_dset_id(i),solv_dspace_id(i),err)
 call h5sget_simple_extent_dims_f(solv_dspace_id(i),dims_2,maxdims_2,err)
 call h5dget_create_plist_f(solv_dset_id(i),solv_crp_list(i),err)
 call h5pget_layout_f(solv_crp_list(i),layout,err)
 if ( h5d_chunked_f == layout ) then
  call h5pget_chunk_f(solv_crp_list(i),rank(2),solv_chunk_2,err)
 end if
 call h5screate_simple_f(rank(2),dims_2,solv_dspace_id(i),err,maxdims_2)
end do

! Dimer group
do i = 1,7,1
 call h5dopen_f(sub_group_id(4),dimer_dset_names(i),dimer_dset_id(i),err)
 call h5dget_space_f(dimer_dset_id(i),dimer_dspace_id(i),err)
 if ( i /= 7 ) then
  call h5sget_simple_extent_dims_f(dimer_dspace_id(i),dims_3,maxdims_3,err)
  call h5dget_create_plist_f(dimer_dset_id(i),dimer_crp_list(i),err)
  call h5pget_layout_f(dimer_crp_list(i),layout,err)
  if ( h5d_chunked_f == layout ) then
   call h5pget_chunk_f(dimer_crp_list(i),rank(3),dimer_chunk_2,err)
  end if
  call h5screate_simple_f(rank(3),dims_3,dimer_dspace_id(i),err,maxdims_3)
 else if ( i == 7 ) then
  call h5sget_simple_extent_dims_f(dimer_dspace_id(i),dims_3,maxdims_3,err)
  call h5dget_create_plist_f(dimer_dset_id(i),dimer_crp_list(i),err)
  call h5pget_layout_f(dimer_crp_list(i),layout,err)
  if ( h5d_chunked_f == layout ) then
   call h5pget_chunk_f(dimer_crp_list(i),rank(3),dimer_chunk_1,err)
  end if
  call h5screate_simple_f(rank(3),dims_3,dimer_dspace_id(i),err,maxdims_3)
 end if
end do

! Read data for restart
! Solvent
solvent_read_dims = (/NT/)
solvent_extend_offset = (/0,int(cstep/Freq1)/)
solvent_extend_count = (/NT,1/)
do i = 1,8,1
 call h5screate_simple_f(rank(1),solvent_read_dims,solvent_memspace,err)
 call h5sselect_hyperslab_f(solv_dspace_id(i),h5s_select_set_f,solvent_extend_offset,solvent_extend_count,err)
 if ( i == 1 ) then             ! xsol
  call h5dread_f(solv_dset_id(i),h5t_native_double,xsol,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 else if ( i == 2 ) then
  call h5dread_f(solv_dset_id(i),h5t_native_double,ysol,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 else if ( i == 3 ) then
  call h5dread_f(solv_dset_id(i),h5t_native_double,zsol,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 else if ( i == 4 ) then
  call h5dread_f(solv_dset_id(i),h5t_native_integer,flag,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 else if ( i == 5 ) then
         call h5dread_f(solv_dset_id(i),h5t_native_integer,type_B,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 else if ( i == 6 ) then
  call h5dread_f(solv_dset_id(i),h5t_native_double,vxsol,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 else if ( i == 7 ) then
  call h5dread_f(solv_dset_id(i),h5t_native_double,vysol,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 else if ( i == 8 ) then
  call h5dread_f(solv_dset_id(i),h5t_native_double,vzsol,solvent_read_dims,err,solvent_memspace,solv_dspace_id(i))
 end if
 call h5sclose_f(solvent_memspace,err)
end do 
! Dimer
dimer_read_dims = (/numdims,2/)
dimer_dimer_extend_offset = (/0,0,int(cstep/Freq4)/)
dimer_dimer_extend_count = (/numdims,2,1/)
do i = 1,6,1
 call h5screate_simple_f(rank(2),dimer_read_dims,dimer_memspace,err)
 call h5sselect_hyperslab_f(dimer_dspace_id(i),h5s_select_set_f,dimer_dimer_extend_offset,dimer_dimer_extend_count,err)
 if ( i == 1 ) then
  call h5dread_f(dimer_dset_id(i),h5t_native_double,xdim,dimer_read_dims,err,dimer_memspace,dimer_dspace_id(i))
 else if ( i == 2 ) then
  call h5dread_f(dimer_dset_id(i),h5t_native_double,ydim,dimer_read_dims,err,dimer_memspace,dimer_dspace_id(i))
 else if ( i == 3 ) then
  call h5dread_f(dimer_dset_id(i),h5t_native_double,zdim,dimer_read_dims,err,dimer_memspace,dimer_dspace_id(i))
 else if ( i == 4 ) then
  call h5dread_f(dimer_dset_id(i),h5t_native_double,vxdim,dimer_read_dims,err,dimer_memspace,dimer_dspace_id(i))
 else if ( i == 5 ) then
  call h5dread_f(dimer_dset_id(i),h5t_native_double,vydim,dimer_read_dims,err,dimer_memspace,dimer_dspace_id(i))
 else if ( i == 6 ) then
  call h5dread_f(dimer_dset_id(i),h5t_native_double,vzdim,dimer_read_dims,err,dimer_memspace,dimer_dspace_id(i))
 end if
 call h5sclose_f(dimer_memspace,err)
end do

open(unit=10,file="motor_positions.dat",status="old")
do i = 1 , numdims
 read(10,*)motorCM(i,1),motorCM(i,2),motorCM(i,3)
end do
close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Force Calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sorting algorithm
if ( mod(cstep,sort_step) == 0 ) then
 call sort_hilbert(vol,NT,box,xsol,ysol,zsol,flag,type_B,vxsol,vysol,vzsol)
end if

!!!!!!Set forces and potential to 0 at time t+h
FXdim = 0.
FYdim = 0.
FZdim = 0.
fwallz = 0.
wall_pot = 0.
FpairX = 0.
FpairY = 0.
FpairZ = 0.
Potential = 0.
TPotential = 0.

! Equate total number of neighbours and neighbour list to zero
neighbour_num = 0
neighbour = 0
outside_particle_num = 0
outside_particle = 0
inside_particle_num = 0
inside_particle = 0
do k = 1,numdims
 ! Catalytic sphere
 call neighbour_list(NT,Box,cutoffC,skinrad,Xsol,Ysol,Zsol,Xdim(k,1),Ydim(k,1),Zdim(k,1),neighbour_num(k,1),neighbour(k,1,:))
 ! Noncatalytic sphere
 call neighbour_list(NT,Box,cutoffN,skinrad,Xsol,Ysol,Zsol,Xdim(k,2),Ydim(k,2),Zdim(k,2),neighbour_num(k,2),neighbour(k,2,:))
end do
call inout_particle_list(NT,box,cutoffC,cutoffN,skinrad,xsol,ysol,zsol,xdim,ydim,zdim,&
outside_particle_num,outside_particle,inside_particle_num,inside_particle,numdims)
max_solvent_dist = 0
time_penetrate = floor((0.5*skinrad / (sqrt(maxval(VXsol*VXsol+VYsol*VYsol+VZsol*VZsol)))) / MDt)
time_count = 0
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
! Calculate total potential and force on dimer
do k = 1,numdims
 do i = 1,2
  do j = 1,neighbour_num(k,i)

   TPotential = TPotential + Potential(k,i,j)
   FXdim(k,i) = FXdim(k,i) - FpairX(k,i,j)                          ! Addition of LJ force on dim., x
   FYdim(k,i) = FYdim(k,i) - FpairY(k,i,j)                          ! Addition of LJ force on dim., y
   FZdim(k,i) = FZdim(k,i) - FpairZ(k,i,j)                          ! Addition of LJ force on dim., z
  end do
  TPotential = TPotential + wall_pot(k,i)
  FZdim(k,i) = FZdim(k,i) + fwallz(k,i)
 end do
end do
end subroutine restart

end module restart_mod
