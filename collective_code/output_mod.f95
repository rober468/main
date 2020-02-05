module output_mod
use hdf5
implicit none

! Precision
integer, parameter :: sgl = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dbl = selected_real_kind(15,307)           ! Double precision

! Non-HDF5 variables
integer :: i,j,k

! HDF5 file writing
! Overall File
character(len=7) , parameter :: datafile = "data.h5"
integer(hid_t) :: data_id
integer(kind=sgl) :: err
integer(kind=sgl) , dimension(3) , parameter :: rank = (/1,2,3/)
integer(hsize_t) , dimension(1) :: maxdims_1
integer(hsize_t) , dimension(2) :: maxdims_2
integer(hsize_t) , dimension(3) :: maxdims_3
integer(hsize_t) , dimension(1) , parameter :: sys_chunk_1 = (/75000/)
integer(hsize_t) , dimension(1) , parameter :: conc_chunk_1 = (/75000/)
integer(hsize_t) , dimension(2) , parameter :: conc1d_chunk_2 = (/120,800/)
integer(hsize_t) , dimension(3) , parameter :: conc2d_chunk_3 = (/120,120,7/)
integer(hsize_t) , dimension(2) , parameter :: solv_chunk_2 = (/100000,1/)
integer(hid_t) , dimension(6) :: data_attr_id
integer(hsize_t) , dimension(2) , parameter :: data_attr_box_dims = (/3,1/)
integer(hsize_t) , dimension(1) , parameter :: data_attr_1_dims = 1
integer(hid_t) , dimension(6) :: data_attr_dspace_id
character(len=21) , dimension(6) , parameter :: data_attr_names =(/"simulation_box_length",&
"solvent_temperature  ","MD_time              ","MPC_time             ","total_#_particles    ",&
"number_of_dimers     "/)
! Group initiation variables
integer(hid_t) , dimension(4) :: sub_group_id
character(len=13) , dimension(4) , parameter :: sub_group_names = (/"Solvent      ",&
"Concentration","System       ","Dimer        "/)
! System group
integer(hid_t) :: sys_attr_id
integer(hsize_t) , dimension(1) , parameter :: sys_attr_1_dims = 1
integer(hid_t) :: sys_attr_dspace_id
character(len=12) , parameter :: sys_attr_freq = "io_frequency"
integer(hid_t) , dimension(6) :: sys_dspace_id
integer(hsize_t) , dimension(1) , parameter :: system_init_dims = (/1/)
integer(hid_t) , dimension(6) :: sys_crp_list
character(len=19) , dimension(6) , parameter :: sys_dset_names = (/"pxsol              ",&
"pysol              ","pzsol              ","total_energy       ","solvent_temperature",&
"dimer_temperature  "/)
integer(hid_t) , dimension(6) :: sys_dset_id
! Concentration group
integer(hid_t) :: conc_attr_id
integer(hsize_t) , dimension(1) , parameter :: conc_attr_1_dims = 1
integer(hid_t) :: conc_attr_dspace_id
character(len=12) , parameter :: conc_attr_freq = "io_frequency"
integer(hid_t) , dimension(12) :: conc_dspace_id
integer(hsize_t) , dimension(1) , parameter :: conc_init_dims = (/1/)
integer(hsize_t) , dimension(2) :: conc1d_init_dims
integer(hsize_t) , dimension(3) :: conc2d_init_dims_xy
integer(hsize_t) , dimension(3) :: conc2d_init_dims_xz
integer(hsize_t) , dimension(3) :: conc2d_init_dims_yz
integer(hid_t) , dimension(12) :: conc_crp_list
character(len=21) , dimension(12) , parameter :: conc_dset_names = (/"total_A              ",&
"total_B              ","concentration_A      ","concentration_B      ",&
"1D_concentration_A   ","1D_concentration_B   ",&
"2D_concentration_A_xy","2D_concentration_B_xy","2D_concentration_A_xz",&
"2D_concentration_B_xz","2D_concentration_A_yz","2D_concentration_B_yz"/)
integer(hid_t) , dimension(12) :: conc_dset_id
! Solvent group
integer(hid_t) :: solv_attr_id
integer(hsize_t) , dimension(1) , parameter :: solv_attr_1_dims = 1
integer(hid_t) :: solv_attr_dspace_id
character(len=12) , parameter :: solv_attr_freq = "io_frequency"
integer(hid_t) , dimension(7) :: solv_dspace_id
integer(hsize_t) , dimension(2) :: solvent_init_dims
integer(hid_t) , dimension(7) :: solv_crp_list
character(len=10) , dimension(7) , parameter :: solv_dset_names = (/"position_x",&
"position_y","position_z","flag      ","velocity_x","velocity_y","velocity_z"/)
integer(hid_t) , dimension(7) :: solv_dset_id
! Dimer group
integer(hid_t) :: dimer_attr_id
integer(hsize_t) , dimension(1) , parameter :: dimer_attr_1_dims = 1
integer(hid_t) :: dimer_attr_dspace_id
character(len=12) , parameter :: dimer_attr_freq = "io_frequency"
integer(hid_t) , dimension(7) :: dimer_dspace_id
integer(hsize_t) , dimension(3) :: dimer_init_dims_dimer
integer(hsize_t) , dimension(3) :: dimer_init_dims_nucvec
integer(hid_t) , dimension(7) :: dimer_crp_list
character(len=10) , dimension(7) , parameter :: dimer_dset_names = (/"position_x",&
"position_y","position_z","velocity_x","velocity_y","velocity_z","nucvec    "/)
integer(hid_t) , dimension(7) :: dimer_dset_id


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_hdf5_initiate(cstep,Freq1,Freq2,Freq3,Freq4,box,T,MDt,MPCt,NT,numdims,x_width,y_width,z_width)
implicit none

! Non-HDF5 variables
integer(kind=sgl) , intent(in) :: cstep
integer(kind=sgl) , intent(in) :: Freq1
integer(kind=sgl) , intent(in) :: Freq2
integer(kind=sgl) , intent(in) :: Freq3
integer(kind=sgl) , intent(in) :: Freq4
real(kind=dbl) , dimension(3) , intent(in) :: box
real(kind=dbl) , intent(in) :: T
real(kind=dbl) , intent(in) :: MDt
real(kind=dbl) , intent(in) :: MPCt
integer(kind=sgl) , intent(in) :: NT
integer , intent(in) :: numdims
real(kind=dbl) :: x_width                                        ! Width of x slice, concentration fields
real(kind=dbl) :: y_width                                        ! Width of y slice, concentration fields
real(kind=dbl) :: z_width                                        ! Width of z slice, concentration fields
integer(hsize_t) , dimension(3) :: dimer_chunk_1
integer(hsize_t) , dimension(3) :: dimer_chunk_2

! Open HDF5
call h5open_f(err)
maxdims_1 = (/h5s_unlimited_f/)
maxdims_2 = (/h5s_unlimited_f,h5s_unlimited_f/)
maxdims_3 = (/h5s_unlimited_f,h5s_unlimited_f,h5s_unlimited_f/)
solvent_init_dims = (/NT,1/)
conc1d_init_dims = (/int(box(1)*2.),1/)
conc2d_init_dims_xy = (/int(box(1)*x_width),int(box(2)*y_width),1/)
conc2d_init_dims_xz = (/int(box(1)*x_width),int(box(3)*z_width),1/)
conc2d_init_dims_yz = (/int(box(2)*y_width),int(box(3)*z_width),1/)
dimer_init_dims_dimer = (/numdims,2,1/)
dimer_init_dims_nucvec = (/numdims,3,1/)
dimer_chunk_1 = (/numdims,3,floor(100000./(3*numdims))/)
dimer_chunk_2 = (/numdims,2,floor(100000./(numdims*2))/)

! Create file
call h5fcreate_f(datafile,h5f_acc_trunc_f,data_id,err)
! General simulation attributes
! Simulation box length
call h5screate_simple_f(rank(2),data_attr_box_dims,data_attr_dspace_id(1),err)
call h5acreate_f(data_id,data_attr_names(1),h5t_native_double,data_attr_dspace_id(1),data_attr_id(1),err)
call h5awrite_f(data_attr_id(1),h5t_native_double,box,data_attr_box_dims,err)
call h5aclose_f(data_attr_id(1),err)
! Others
do i = 2,6,1
 call h5screate_simple_f(rank(1),data_attr_1_dims,data_attr_dspace_id(i),err)
 call h5acreate_f(data_id,data_attr_names(i),h5t_native_double,data_attr_dspace_id(i),data_attr_id(i),err)
 if ( i == 2 ) then             ! Temperature
  call h5awrite_f(data_attr_id(i),h5t_native_double,T,data_attr_1_dims,err)
 else if ( i == 3 ) then        ! MD time
  call h5awrite_f(data_attr_id(i),h5t_native_double,MDt,data_attr_1_dims,err)
 else if ( i == 4 ) then        ! MPC time
  call h5awrite_f(data_attr_id(i),h5t_native_double,MPCt,data_attr_1_dims,err)
 else if ( i == 5 ) then        ! Total number of particles
  call h5awrite_f(data_attr_id(i),h5t_native_integer,NT,data_attr_1_dims,err)
 else if ( i == 6 ) then        ! Total number of particles
  call h5awrite_f(data_attr_id(i),h5t_native_integer,numdims,data_attr_1_dims,err)
 end if
call h5aclose_f(data_attr_id(i),err)
end do
! Create subgroups
if ( numdims > 0 ) then
 do i = 1,4,1
  call h5gcreate_f(data_id,sub_group_names(i),sub_group_id(i),err)
 end do
else
 do i = 1,3,1
  call h5gcreate_f(data_id,sub_group_names(i),sub_group_id(i),err)
 end do
end if

! Create groups, dataspaces,  datasets

! System group
! Group attributes
call h5screate_simple_f(rank(1),sys_attr_1_dims,sys_attr_dspace_id,err)
call h5acreate_f(sub_group_id(3),sys_attr_freq,h5t_native_integer,sys_attr_dspace_id,sys_attr_id,err)
call h5awrite_f(sys_attr_id,h5t_native_integer,Freq3,sys_attr_1_dims,err)
call h5aclose_f(sys_attr_id,err)
do i = 1,6,1
 ! Dataspace creation
 call h5screate_simple_f(rank(1),system_init_dims,sys_dspace_id(i),err,maxdims_1)
 ! Dataset creation properties
 call h5pcreate_f(h5p_dataset_create_f,sys_crp_list(i),err)
 call h5pset_chunk_f(sys_crp_list(i),rank(1),sys_chunk_1,err)
 ! Dataset creation
 call h5dcreate_f(sub_group_id(3),sys_dset_names(i),h5t_native_double,sys_dspace_id(i),sys_dset_id(i),err,sys_crp_list(i))
end do

! Concentration group
! Group attributes
call h5screate_simple_f(rank(1),conc_attr_1_dims,conc_attr_dspace_id,err)
call h5acreate_f(sub_group_id(2),conc_attr_freq,h5t_native_integer,conc_attr_dspace_id,conc_attr_id,err)
call h5awrite_f(conc_attr_id,h5t_native_integer,Freq2,conc_attr_1_dims,err)
call h5aclose_f(conc_attr_id,err)
do i = 1,12,1
 ! Dataspace creation, dataset creation property, chunking, dataset creation
 if ( i <= 4 ) then
 call h5screate_simple_f(rank(1),conc_init_dims,conc_dspace_id(i),err,maxdims_1)
 call h5pcreate_f(h5p_dataset_create_f,conc_crp_list(i),err)
 call h5pset_chunk_f(conc_crp_list(i),rank(1),conc_chunk_1,err)
 call h5dcreate_f(sub_group_id(2),conc_dset_names(i),h5t_native_double,conc_dspace_id(i),conc_dset_id(i),err,conc_crp_list(i))
 else if ( i >= 5 .and. i  <= 6 ) then
 call h5screate_simple_f(rank(2),conc1d_init_dims,conc_dspace_id(i),err,maxdims_2)
 call h5pcreate_f(h5p_dataset_create_f,conc_crp_list(i),err)
 call h5pset_chunk_f(conc_crp_list(i),rank(2),conc1d_chunk_2,err)
 call h5dcreate_f(sub_group_id(2),conc_dset_names(i),h5t_native_double,conc_dspace_id(i),conc_dset_id(i),err,conc_crp_list(i))
 else if ( i >= 7 .and. i  <= 8 ) then
 call h5screate_simple_f(rank(3),conc2d_init_dims_xy,conc_dspace_id(i),err,maxdims_3)
 call h5pcreate_f(h5p_dataset_create_f,conc_crp_list(i),err)
 call h5pset_chunk_f(conc_crp_list(i),rank(3),conc2d_chunk_3,err)
 call h5dcreate_f(sub_group_id(2),conc_dset_names(i),h5t_native_double,conc_dspace_id(i),conc_dset_id(i),err,conc_crp_list(i))
 else if ( i >= 9 .and. i  <= 10 ) then
 call h5screate_simple_f(rank(3),conc2d_init_dims_xz,conc_dspace_id(i),err,maxdims_3)
 call h5pcreate_f(h5p_dataset_create_f,conc_crp_list(i),err)
 call h5pset_chunk_f(conc_crp_list(i),rank(3),conc2d_chunk_3,err)
 call h5dcreate_f(sub_group_id(2),conc_dset_names(i),h5t_native_double,conc_dspace_id(i),conc_dset_id(i),err,conc_crp_list(i))
 else if ( i >= 11 .and. i  <= 12 ) then
 call h5screate_simple_f(rank(3),conc2d_init_dims_yz,conc_dspace_id(i),err,maxdims_3)
 call h5pcreate_f(h5p_dataset_create_f,conc_crp_list(i),err)
 call h5pset_chunk_f(conc_crp_list(i),rank(3),conc2d_chunk_3,err)
 call h5dcreate_f(sub_group_id(2),conc_dset_names(i),h5t_native_double,conc_dspace_id(i),conc_dset_id(i),err,conc_crp_list(i))
 end if
end do

! Solvent group
! Group attributes
call h5screate_simple_f(rank(1),solv_attr_1_dims,solv_attr_dspace_id,err)
call h5acreate_f(sub_group_id(1),solv_attr_freq,h5t_native_integer,solv_attr_dspace_id,solv_attr_id,err)
call h5awrite_f(solv_attr_id,h5t_native_integer,Freq1,solv_attr_1_dims,err)
call h5aclose_f(solv_attr_id,err)
do i = 1,7,1
 ! Dataspace creation, dataset creation property, chunking, dataset creation
 call h5screate_simple_f(rank(2),solvent_init_dims,solv_dspace_id(i),err,maxdims_2)
 call h5pcreate_f(h5p_dataset_create_f,solv_crp_list(i),err)
 call h5pset_chunk_f(solv_crp_list(i),rank(2),solv_chunk_2,err)
 if ( i /= 4 ) then
 call h5dcreate_f(sub_group_id(1),solv_dset_names(i),h5t_native_double,solv_dspace_id(i),solv_dset_id(i),err,solv_crp_list(i))
 else if ( i == 4 ) then
 call h5dcreate_f(sub_group_id(1),solv_dset_names(i),h5t_native_integer,solv_dspace_id(i),solv_dset_id(i),err,solv_crp_list(i))
 end if
end do

if ( numdims > 0 ) then

! Dimer group
! Group attributes
call h5screate_simple_f(rank(1),dimer_attr_1_dims,dimer_attr_dspace_id,err)
call h5acreate_f(sub_group_id(4),dimer_attr_freq,h5t_native_integer,dimer_attr_dspace_id,dimer_attr_id,err)
call h5awrite_f(dimer_attr_id,h5t_native_integer,Freq4,dimer_attr_1_dims,err)
call h5aclose_f(dimer_attr_id,err)
do i = 1,7,1
 if ( i /= 7 ) then
  call h5screate_simple_f(rank(3),dimer_init_dims_dimer,dimer_dspace_id(i),err,maxdims_3)
  call h5pcreate_f(h5p_dataset_create_f,dimer_crp_list(i),err)
  call h5pset_chunk_f(dimer_crp_list(i),rank(3),dimer_chunk_2,err)
  call h5dcreate_f(sub_group_id(4),dimer_dset_names(i),h5t_native_double,dimer_dspace_id(i),dimer_dset_id(i),err,dimer_crp_list(i))
 else if ( i == 7 ) then
  call h5screate_simple_f(rank(3),dimer_init_dims_nucvec,dimer_dspace_id(i),err,maxdims_3)
  call h5pcreate_f(h5p_dataset_create_f,dimer_crp_list(i),err)
  call h5pset_chunk_f(dimer_crp_list(i),rank(3),dimer_chunk_1,err)
  call h5dcreate_f(sub_group_id(4),dimer_dset_names(i),h5t_native_double,dimer_dspace_id(i),dimer_dset_id(i),err,dimer_crp_list(i))
 end if
end do

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_hdf5_close(numdims)
use hdf5
implicit none

! Precision
integer, parameter :: sgl = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dbl = selected_real_kind(15,307)           ! Double precision

! External variables
integer , intent(in) :: numdims

! System group
do i = 1,6,1
 call h5dclose_f(sys_dset_id(i),err)
 call h5sclose_f(sys_dspace_id(i),err)
end do

! Concentration group
do i = 1,12,1
 call h5dclose_f(conc_dset_id(i),err)
 call h5sclose_f(conc_dspace_id(i),err)
end do

! Solvent group
do i = 1,7,1
 call h5dclose_f(solv_dset_id(i),err)
 call h5sclose_f(solv_dspace_id(i),err)
end do

! Dimer group
if ( numdims > 0 ) then
do i = 1,7,1
 call h5dclose_f(dimer_dset_id(i),err)
 call h5sclose_f(dimer_dspace_id(i),err)
end do
end if

! Close groups
if ( numdims > 0 ) then
 do i = 1,4,1
  call h5gclose_f(sub_group_id(i),err)
 end do
else
 do i = 1,3,1
  call h5gclose_f(sub_group_id(i),err)
 end do
end if

! Close file
call h5fclose_f(data_id,err)

call h5close_f(err)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_hdf5_write(cstep,Vol,Freq1,Freq2,Freq3,Freq4,NT,numdims,Ma,Mb,Ms,Mc,Mn,x_conc_slab,&
                             y_conc_slab,z_conc_slab,box,x_width,y_width,z_width,xsol,ysol,zsol,flag,&
                             vxsol,vysol,vzsol,cutoffC,Dcn,concA1d,concB1d,concA2d_xy,concB2d_xy,concA2d_xz,&
                             concB2d_xz,concA2d_yz,concB2d_yz,pxsol,pysol,pzsol,xdim,ydim,zdim,vxdim,vydim,vzdim,&
                             Tpotential)
implicit none

! External variables
integer(kind=sgl) , intent(in) :: cstep
integer(kind=sgl) , intent(in) :: Freq1
integer(kind=sgl) , intent(in) :: Freq2
integer(kind=sgl) , intent(in) :: Freq3
integer(kind=sgl) , intent(in) :: Freq4
integer(kind=sgl) , intent(in) :: NT
integer(kind=dbl) , intent(in) :: Vol
integer , intent(in) :: numdims
real(kind=dbl) , intent(in) :: Ma
real(kind=dbl) , intent(in) :: Mb
real(kind=dbl) , intent(in) :: Ms
real(kind=dbl) , intent(in) :: Mc
real(kind=dbl) , intent(in) :: Mn
real(kind=dbl) , intent(in) :: x_conc_slab
real(kind=dbl) , intent(in) :: y_conc_slab
real(kind=dbl) , intent(in) :: z_conc_slab
real(kind=dbl) , dimension(3) , intent(in) :: box
real(kind=dbl) , intent(in) :: cutoffC
real(kind=dbl) , intent(in) :: Dcn
real(kind=dbl) , intent(in) :: x_width
real(kind=dbl) , intent(in) :: y_width
real(kind=dbl) , intent(in) :: z_width
real(kind=dbl) , dimension(:) , intent(in) :: xsol
real(kind=dbl) , dimension(:) , intent(in) :: ysol
real(kind=dbl) , dimension(:) , intent(in) :: zsol
integer(kind=sgl) , dimension(:) , intent(in) :: flag
real(kind=dbl) , dimension(:) , intent(in) :: vxsol
real(kind=dbl) , dimension(:) , intent(in) :: vysol
real(kind=dbl) , dimension(:) , intent(in) :: vzsol
real(kind=dbl) , dimension(:) , intent(inout) :: concA1d
real(kind=dbl) , dimension(:) , intent(inout) :: concB1d
real(kind=dbl) , dimension(:,:) , intent(inout) :: concA2d_xy
real(kind=dbl) , dimension(:,:) , intent(inout) :: concB2d_xy
real(kind=dbl) , dimension(:,:) , intent(inout) :: concA2d_xz
real(kind=dbl) , dimension(:,:) , intent(inout) :: concB2d_xz
real(kind=dbl) , dimension(:,:) , intent(inout) :: concA2d_yz
real(kind=dbl) , dimension(:,:) , intent(inout) :: concB2d_yz
real(kind=dbl) , dimension(:) , intent(inout) :: pxsol
real(kind=dbl) , dimension(:) , intent(inout) :: pysol
real(kind=dbl) , dimension(:) , intent(inout) :: pzsol
real(kind=dbl) , dimension(:,:) , intent(in) :: xdim
real(kind=dbl) , dimension(:,:) , intent(in) :: ydim
real(kind=dbl) , dimension(:,:) , intent(in) :: zdim
real(kind=dbl) , dimension(:,:) , intent(in) :: vxdim
real(kind=dbl) , dimension(:,:) , intent(in) :: vydim
real(kind=dbl) , dimension(:,:) , intent(in) :: vzdim
real(kind=dbl) , intent(in) :: TPotential

! Internal Variables
integer(hsize_t) , dimension(2) :: solvent_extend
integer(hsize_t) , dimension(2) :: solvent_extend_offset
integer(hsize_t) , dimension(2) :: solvent_extend_count
integer(hid_t) :: solvent_memspace
integer(hsize_t) , dimension(1) :: concentration_extend
integer(hsize_t) , dimension(1) :: concentration_extend_offset
integer(hsize_t) , dimension(1) :: concentration_extend_count
integer(hsize_t) , dimension(2) :: conc1d_extend
integer(hsize_t) , dimension(2) :: conc1d_extend_offset
integer(hsize_t) , dimension(2) :: conc1d_extend_count
integer(hsize_t) , dimension(3) :: conc2d_extend_xy
integer(hsize_t) , dimension(3) :: conc2d_extend_offset_xy
integer(hsize_t) , dimension(3) :: conc2d_extend_count_xy
integer(hsize_t) , dimension(3) :: conc2d_extend_xz
integer(hsize_t) , dimension(3) :: conc2d_extend_offset_xz
integer(hsize_t) , dimension(3) :: conc2d_extend_count_xz
integer(hsize_t) , dimension(3) :: conc2d_extend_yz
integer(hsize_t) , dimension(3) :: conc2d_extend_offset_yz
integer(hsize_t) , dimension(3) :: conc2d_extend_count_yz
integer(hid_t) :: concentration_memspace
integer(hsize_t) , dimension(1) :: system_extend
integer(hsize_t) , dimension(1) :: system_extend_offset
integer(hsize_t) , dimension(1) :: system_extend_count
integer(hid_t) :: system_memspace
integer(hsize_t) , dimension(3) :: dimer_dimer_extend
integer(hsize_t) , dimension(3) :: dimer_dimer_extend_offset
integer(hsize_t) , dimension(3) :: dimer_dimer_extend_count
integer(hsize_t) , dimension(3) :: dimer_nucvec_extend
integer(hsize_t) , dimension(3) :: dimer_nucvec_extend_offset
integer(hsize_t) , dimension(3) :: dimer_nucvec_extend_count
integer(hid_t) :: dimer_memspace
real(kind=dbl) :: vol1d
real(kind=dbl) :: vol2d_xy
real(kind=dbl) :: vol2d_xz
real(kind=dbl) :: vol2d_yz
real(kind=dbl) :: TotA
real(kind=dbl) :: TotB
real(kind=dbl) :: TotS
real(kind=dbl) :: ConcA
real(kind=dbl) :: ConcB
real(kind=dbl) :: sumPX
real(kind=dbl) :: sumPY
real(kind=dbl) :: sumPZ
real(kind=dbl) , dimension(numdims,2) :: pxdim
real(kind=dbl) , dimension(numdims,2) :: pydim
real(kind=dbl) , dimension(numdims,2) :: pzdim
real(kind=dbl) :: KEsol
real(kind=dbl) :: KEdim
real(kind=dbl) :: Ekin
real(kind=dbl) :: Epot
real(kind=dbl) :: tot_energy
real(kind=dbl) :: tempsol
real(kind=dbl) :: tempdim
real(kind=dbl) , dimension(3) :: diffdim
real(kind=dbl) , dimension(numdims,3) :: nucvec
integer(kind=sgl) :: err

! Record solvent position,velocities,identity
if ( mod(cstep,Freq1) == 0 ) then
 
 if ( cstep == 0 ) then

 solvent_init_dims = (/NT,1/)
 call h5dwrite_f(solv_dset_id(1),h5t_native_double,xsol,solvent_init_dims,err)
 call h5dwrite_f(solv_dset_id(2),h5t_native_double,ysol,solvent_init_dims,err)
 call h5dwrite_f(solv_dset_id(3),h5t_native_double,zsol,solvent_init_dims,err)
 call h5dwrite_f(solv_dset_id(4),h5t_native_integer,flag,solvent_init_dims,err)
 call h5dwrite_f(solv_dset_id(5),h5t_native_double,vxsol,solvent_init_dims,err)
 call h5dwrite_f(solv_dset_id(6),h5t_native_double,vysol,solvent_init_dims,err)
 call h5dwrite_f(solv_dset_id(7),h5t_native_double,vzsol,solvent_init_dims,err)

else

 solvent_extend = (/NT,int(cstep/Freq1)+1/)
 solvent_extend_offset = (/0,int(cstep/Freq1)/)
 solvent_extend_count = (/NT,1/)

 do i = 1,7,1
  call h5dset_extent_f(solv_dset_id(i),solvent_extend,err)
  call h5screate_simple_f(rank(2),solvent_extend_count,solvent_memspace,err)
  call h5dget_space_f(solv_dset_id(i),solv_dspace_id(i),err)
  call h5sselect_hyperslab_f(solv_dspace_id(i),h5s_select_set_f,solvent_extend_offset,solvent_extend_count,err)
  if ( i == 1 ) then            ! xsol
   call h5dwrite_f(solv_dset_id(i),h5t_native_double,xsol,solvent_extend_count,err,solvent_memspace,solv_dspace_id(i))
  else if ( i == 2 ) then       ! ysol
   call h5dwrite_f(solv_dset_id(i),h5t_native_double,ysol,solvent_extend_count,err,solvent_memspace,solv_dspace_id(i))
  else if ( i == 3 ) then       ! zsol
   call h5dwrite_f(solv_dset_id(i),h5t_native_double,zsol,solvent_extend_count,err,solvent_memspace,solv_dspace_id(i))
  else if ( i == 4 ) then       ! flag
   call h5dwrite_f(solv_dset_id(i),h5t_native_integer,flag,solvent_extend_count,err,solvent_memspace,solv_dspace_id(i))
  else if ( i == 5 ) then       ! vxsol
   call h5dwrite_f(solv_dset_id(i),h5t_native_double,vxsol,solvent_extend_count,err,solvent_memspace,solv_dspace_id(i))
  else if ( i == 6 ) then       ! vysol
   call h5dwrite_f(solv_dset_id(i),h5t_native_double,vysol,solvent_extend_count,err,solvent_memspace,solv_dspace_id(i))
  else if ( i == 7 ) then       ! vzsol
   call h5dwrite_f(solv_dset_id(i),h5t_native_double,vzsol,solvent_extend_count,err,solvent_memspace,solv_dspace_id(i))
  end if
  call h5sclose_f(solvent_memspace,err)
 end do

end if

end if


if ( mod(cstep,Freq2) == 0 ) then

! System concentration
TotA = 0.
TotB = 0.

do i = 1,NT
 if ( flag(i) == 0 ) then
  totA = TotA + 1.
 else if ( flag(i) == 1 ) then
  totB = totB + 1.
 end if
end do

ConcA = TotA / Vol
ConcB = TotB / Vol

! 1D concentration (along x axis; 2 removed from upper and lower part of z BB
! boundary
concA1d = 0.
concB1d = 0.

vol1d = x_width*Box(2)*( 2.*cutoffC )
do i = 1,NT
 if ( zsol(i) <= ( z_conc_slab + cutoffC ) .and. zsol(i) >= ( z_conc_slab - cutoffC ) ) then
  j = ceiling(xsol(i)/x_width)
  if ( flag(i) == 0 ) then
   concA1d(j) = concA1d(j) + 1.
  else if ( flag(i) == 1 ) then
   concB1d(j) = concB1d(j) + 1.
  end if
 end if
end do

concA1d = concA1d / vol1d
concB1d = concB1d / vol1d

! 2D concentration (along xy plane)
concA2d_xy = 0.
concB2d_xy = 0.

vol2d_xy = x_width*y_width*( 2.*cutoffC  )
do i = 1,NT
 if ( zsol(i) <= ( z_conc_slab + cutoffC ) .and. zsol(i) >= ( z_conc_slab - cutoffC ) ) then
  j = ceiling(xsol(i)/x_width)
  k = ceiling(ysol(i)/y_width)
  if ( flag(i) == 0 ) then
   concA2d_xy(j,k) = concA2d_xy(j,k) + 1.
  else if ( flag(i) == 1 ) then
   concB2d_xy(j,k) = concB2d_xy(j,k) + 1.
  end if
 end if
end do

concA2d_xy = concA2d_xy / vol2d_xy
concB2d_xy = concB2d_xy / vol2d_xy

! 2D concentration (along xz plane)
concA2d_xz = 0.
concB2d_xz = 0.

vol2d_xz = x_width*z_width*( 2.*cutoffC  )
do i = 1,NT
 if ( ysol(i) <= ( y_conc_slab + cutoffC ) .and. ysol(i) >= ( y_conc_slab - cutoffC ) ) then
  j = ceiling(xsol(i)/x_width)
  k = ceiling(zsol(i)/z_width)
  if ( flag(i) == 0 ) then
   concA2d_xz(j,k) = concA2d_xz(j,k) + 1.
  else if ( flag(i) == 1 ) then
   concB2d_xz(j,k) = concB2d_xz(j,k) + 1.
  end if
 end if
end do

concA2d_xz = concA2d_xz / vol2d_xz
concB2d_xz = concB2d_xz / vol2d_xz

! 2D concentration (along yz plane)
concA2d_yz = 0.
concB2d_yz = 0.

vol2d_yz = y_width*z_width*( 2.*cutoffC  )
do i = 1,NT
 if ( xsol(i) <= ( x_conc_slab + cutoffC ) .and. xsol(i) >= ( x_conc_slab - cutoffC ) ) then
  j = ceiling(ysol(i)/y_width)
  k = ceiling(zsol(i)/z_width)
  if ( flag(i) == 0 ) then
   concA2d_yz(j,k) = concA2d_yz(j,k) + 1.
  else if ( flag(i) == 1 ) then
   concB2d_yz(j,k) = concB2d_yz(j,k) + 1.
  end if
 end if
end do

concA2d_yz = concA2d_yz / vol2d_yz
concB2d_yz = concB2d_yz / vol2d_yz

! Record total number and concentration of A and B species

if ( cstep == 0 ) then

conc1d_init_dims = (/int(box(1)/x_width),1/)
conc2d_init_dims_xy = (/int(box(1)/x_width),int(box(2)/y_width),1/)
conc2d_init_dims_xz = (/int(box(1)/x_width),int(box(3)/z_width),1/)
conc2d_init_dims_yz = (/int(box(2)/y_width),int(box(3)/z_width),1/)

call h5dwrite_f(conc_dset_id(1),h5t_native_double,TotA,conc_init_dims,err)
call h5dwrite_f(conc_dset_id(2),h5t_native_double,TotB,conc_init_dims,err)
call h5dwrite_f(conc_dset_id(3),h5t_native_double,ConcA,conc_init_dims,err)
call h5dwrite_f(conc_dset_id(4),h5t_native_double,ConcB,conc_init_dims,err)
call h5dwrite_f(conc_dset_id(5),h5t_native_double,concA1d,conc1d_init_dims,err)
call h5dwrite_f(conc_dset_id(6),h5t_native_double,concB1d,conc1d_init_dims,err)
call h5dwrite_f(conc_dset_id(7),h5t_native_double,concA2d_xy,conc2d_init_dims_xy,err)
call h5dwrite_f(conc_dset_id(8),h5t_native_double,concB2d_xy,conc2d_init_dims_xy,err)
call h5dwrite_f(conc_dset_id(9),h5t_native_double,concA2d_xz,conc2d_init_dims_xz,err)
call h5dwrite_f(conc_dset_id(10),h5t_native_double,concB2d_xz,conc2d_init_dims_xz,err)
call h5dwrite_f(conc_dset_id(11),h5t_native_double,concA2d_yz,conc2d_init_dims_yz,err)
call h5dwrite_f(conc_dset_id(12),h5t_native_double,concB2d_yz,conc2d_init_dims_yz,err)

else
concentration_extend = (/int(cstep/Freq2)+1/)
concentration_extend_offset = (/int(cstep/Freq2)/)
concentration_extend_count = (/1/)
conc1d_extend = (/int(box(1)/x_width),int(cstep/Freq2)+1/)
conc1d_extend_offset = (/0,int(cstep/Freq2)/)
conc1d_extend_count = (/int(box(1)/x_width),1/)
conc2d_extend_xy = (/int(box(1)/x_width),int(box(2)/y_width),int(cstep/Freq2)+1/)
conc2d_extend_offset_xy = (/0,0,int(cstep/Freq2)/)
conc2d_extend_count_xy = (/int(box(1)/x_width),int(box(2)/y_width),1/)
conc2d_extend_xz = (/int(box(1)/x_width),int(box(3)/z_width),int(cstep/Freq2)+1/)
conc2d_extend_offset_xz = (/0,0,int(cstep/Freq2)/)
conc2d_extend_count_xz = (/int(box(1)/x_width),int(box(3)/z_width),1/)
conc2d_extend_yz = (/int(box(2)/y_width),int(box(3)/z_width),int(cstep/Freq2)+1/)
conc2d_extend_offset_yz = (/0,0,int(cstep/Freq2)/)
conc2d_extend_count_yz = (/int(box(2)/y_width),int(box(3)/z_width),1/)

do i = 1,12,1
 if ( i > 0 .and. i <= 4 ) then
  call h5dset_extent_f(conc_dset_id(i),concentration_extend,err)
  call h5screate_simple_f(rank(1),concentration_extend_count,concentration_memspace,err)
  call h5dget_space_f(conc_dset_id(i),conc_dspace_id(i),err)
  call h5sselect_hyperslab_f(conc_dspace_id(i),h5s_select_set_f,concentration_extend_offset,concentration_extend_count,err)
 else if ( i >= 5 .and. i <= 6 ) then
  call h5dset_extent_f(conc_dset_id(i),conc1d_extend,err)
  call h5screate_simple_f(rank(2),conc1d_extend_count,concentration_memspace,err)
  call h5dget_space_f(conc_dset_id(i),conc_dspace_id(i),err)
  call h5sselect_hyperslab_f(conc_dspace_id(i),h5s_select_set_f,conc1d_extend_offset,conc1d_extend_count,err)
 else if ( i >= 7 .and. i <= 8 ) then
  call h5dset_extent_f(conc_dset_id(i),conc2d_extend_xy,err)
  call h5screate_simple_f(rank(3),conc2d_extend_count_xy,concentration_memspace,err)
  call h5dget_space_f(conc_dset_id(i),conc_dspace_id(i),err)
  call h5sselect_hyperslab_f(conc_dspace_id(i),h5s_select_set_f,conc2d_extend_offset_xy,conc2d_extend_count_xy,err)
 else if ( i >= 9 .and. i <= 10 ) then
  call h5dset_extent_f(conc_dset_id(i),conc2d_extend_xz,err)
  call h5screate_simple_f(rank(3),conc2d_extend_count_xz,concentration_memspace,err)
  call h5dget_space_f(conc_dset_id(i),conc_dspace_id(i),err)
  call h5sselect_hyperslab_f(conc_dspace_id(i),h5s_select_set_f,conc2d_extend_offset_xz,conc2d_extend_count_xz,err)
 else if ( i >= 11 .and. i <= 12 ) then
  call h5dset_extent_f(conc_dset_id(i),conc2d_extend_yz,err)
  call h5screate_simple_f(rank(3),conc2d_extend_count_yz,concentration_memspace,err)
  call h5dget_space_f(conc_dset_id(i),conc_dspace_id(i),err)
  call h5sselect_hyperslab_f(conc_dspace_id(i),h5s_select_set_f,conc2d_extend_offset_yz,conc2d_extend_count_yz,err)
 end if
 if ( i == 1 ) then             ! TotA
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,TotA,concentration_extend_count,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 2 ) then        ! TotB
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,TotB,concentration_extend_count,err,concentration_memspace,conc_dspace_id(i))
  else if ( i == 3 ) then       ! ConcA
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,ConcA,concentration_extend_count,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 4 ) then        ! ConcB
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,ConcB,concentration_extend_count,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 5 ) then        ! ConcA1d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concA1d,conc1d_extend_count,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 6 ) then        ! ConcB1d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concB1d,conc1d_extend_count,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 7 ) then        ! ConcA2d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concA2d_xy,conc2d_extend_count_xy,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 8 ) then       ! ConcB2d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concB2d_xy,conc2d_extend_count_xy,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 9 ) then        ! ConcA2d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concA2d_xz,conc2d_extend_count_xz,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 10 ) then       ! ConcB2d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concB2d_xz,conc2d_extend_count_xz,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 11 ) then        ! ConcA2d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concA2d_yz,conc2d_extend_count_yz,err,concentration_memspace,conc_dspace_id(i))
 else if ( i == 12 ) then       ! ConcB2d
  call h5dwrite_f(conc_dset_id(i),h5t_native_double,concB2d_yz,conc2d_extend_count_yz,err,concentration_memspace,conc_dspace_id(i))
 end if
 call h5sclose_f(concentration_memspace,err)
end do

end if

end if

if ( mod(cstep,Freq3) == 0 ) then

pxdim = 0.
pydim = 0.
pzdim = 0.
sumpx = 0.
sumpy = 0.
sumpz = 0.
KEdim = 0.
KEsol = 0.
TotA = 0.
TotB = 0.
TotS = 0.

! Energy and momentum calculation
do i = 1,NT
 if ( flag(i) == 0 ) then
  TotA = TotA + 1.
  PXsol(i) = Ma * VXsol(i)
  PYsol(i) = Ma * VYsol(i)
  PZsol(i) = Ma * VZsol(i)
 else if ( flag(i) == 1 ) then
  TotB = TotB + 1.
  PXsol(i) = Mb * VXsol(i)
  PYsol(i) = Mb * VYsol(i)
  PZsol(i) = Mb * VZsol(i)
 else if ( flag(i) == 2 ) then
  TotS = TotS + 1.
  PXsol(i) = Ms * VXsol(i)
  PYsol(i) = Ms * VYsol(i)
  PZsol(i) = Ms * VZsol(i)
 end if
 sumpx = sumpx + pxsol(i)
 sumpy = sumpy + pysol(i)
 sumpz = sumpz + pzsol(i)
end do

if ( numdims > 0 ) then

do i = 1,numdims
 pxdim(i,1) = Mc*vxdim(i,1)
 pydim(i,1) = Mc*vydim(i,1)
 pzdim(i,1) = Mc*vzdim(i,1)
 pxdim(i,2) = Mn*vxdim(i,2)
 pydim(i,2) = Mn*vydim(i,2)
 pzdim(i,2) = Mn*vzdim(i,2)
 sumpx = sumpx + pxdim(i,1) + pxdim(i,2)
 sumpy = sumpy + pydim(i,1) + pydim(i,2)
 sumpz = sumpz + pzdim(i,1) + pzdim(i,2)
 KEdim = KEdim + 0.5d0 * (1./Mc) * ( pxdim(i,1)**2 + pydim(i,1)**2 + pzdim(i,1)**2 ) &
               + 0.5d0 * (1./Mn) * ( pxdim(i,2)**2 + pydim(i,2)**2 + pzdim(i,2)**2 )
end do

end if

do i = 1,NT,1
 KEsol = KEsol + 0.5 * ( VXsol(i)*PXsol(i) + VYsol(i)*PYsol(i) + VZsol(i)*PZsol(i) )
end do

if  ( numdims > 0 ) then

Ekin = KEdim + KEsol
Epot = TPotential
tot_energy = Ekin + Epot

! Temperature calculation
tempsol = (KEsol/(TotA+TotB+TotS)) * (2./3.)
tempdim = KEdim * (2./3.)

else

Ekin = KEsol
Epot = TPotential
tot_energy = Ekin + Epot

! Temperature calculation
tempsol = (KEsol/(TotA+TotB+TotS)) * (2./3.)
tempdim = 0.

end if

! Record total momentum for each direction, total energy, and solvent
! temperature

if ( cstep == 0 ) then

 call h5dwrite_f(sys_dset_id(1),h5t_native_double,sumPX,system_init_dims,err)
 call h5dwrite_f(sys_dset_id(2),h5t_native_double,sumPY,system_init_dims,err)
 call h5dwrite_f(sys_dset_id(3),h5t_native_double,sumPZ,system_init_dims,err)
 call h5dwrite_f(sys_dset_id(4),h5t_native_double,tot_energy,system_init_dims,err)
 call h5dwrite_f(sys_dset_id(5),h5t_native_double,tempsol,system_init_dims,err)
 call h5dwrite_f(sys_dset_id(6),h5t_native_double,tempdim,system_init_dims,err)

else

 system_extend = (/int(cstep/Freq3)+1/)
 system_extend_offset = (/int(cstep/Freq3)/)
 system_extend_count = (/1/)

 do i = 1,6,1
  call h5dset_extent_f(sys_dset_id(i),system_extend,err)
  call h5screate_simple_f(rank(1),system_extend_count,system_memspace,err)
  call h5dget_space_f(sys_dset_id(i),sys_dspace_id(i),err)
  call h5sselect_hyperslab_f(sys_dspace_id(i),h5s_select_set_f,system_extend_offset,system_extend_count,err)
  if ( i == 1 ) then
   call h5dwrite_f(sys_dset_id(i),h5t_native_double,sumPX,system_extend_count,err,system_memspace,sys_dspace_id(i))
  else if ( i == 2 ) then
   call h5dwrite_f(sys_dset_id(i),h5t_native_double,sumPY,system_extend_count,err,system_memspace,sys_dspace_id(i))
  else if ( i == 3 ) then
   call h5dwrite_f(sys_dset_id(i),h5t_native_double,sumPZ,system_extend_count,err,system_memspace,sys_dspace_id(i))
  else if ( i == 4 ) then
   call h5dwrite_f(sys_dset_id(i),h5t_native_double,tot_energy,system_extend_count,err,system_memspace,sys_dspace_id(i))
  else if ( i == 5 ) then
   call h5dwrite_f(sys_dset_id(i),h5t_native_double,tempsol,system_extend_count,err,system_memspace,sys_dspace_id(i))
  else if ( i == 6 ) then
   call h5dwrite_f(sys_dset_id(i),h5t_native_double,tempdim,system_extend_count,err,system_memspace,sys_dspace_id(i))
  end if
  call h5sclose_f(system_memspace,err)
 end do

end if

end if

if ( numdims > 0 ) then

if ( mod(cstep,Freq4) == 0 ) then

! Difference in dimer positions
do j = 1,numdims
 diffdim(1) = Xdim(j,1) - Xdim(j,2)
 diffdim(2) = Ydim(j,1) - Ydim(j,2)
 diffdim(3) = Zdim(j,1) - Zdim(j,2)

! Apply periodic conditions
 do i = 1,2
  if ( diffdim(i) > 0.5d0*Box(i) ) then
   diffdim(i) = diffdim(i) - Box(i)
  else if ( diffdim(i) < -0.5d0*Box(i) ) then
   diffdim(i) = diffdim(i) + Box(i)
  end if
 end do

! Find distance between the two dimers and make unit vector
 do i = 1,3
  nucvec(j,i) = diffdim(i)/Dcn
 end do
end do

! Record total momentum for each direction, total energy, and solvent
! temperature

if ( cstep == 0 ) then

 call h5dwrite_f(dimer_dset_id(1),h5t_native_double,xdim,dimer_init_dims_dimer,err)
 call h5dwrite_f(dimer_dset_id(2),h5t_native_double,ydim,dimer_init_dims_dimer,err)
 call h5dwrite_f(dimer_dset_id(3),h5t_native_double,zdim,dimer_init_dims_dimer,err)
 call h5dwrite_f(dimer_dset_id(4),h5t_native_double,vxdim,dimer_init_dims_dimer,err)
 call h5dwrite_f(dimer_dset_id(5),h5t_native_double,vydim,dimer_init_dims_dimer,err)
 call h5dwrite_f(dimer_dset_id(6),h5t_native_double,vzdim,dimer_init_dims_dimer,err)
 call h5dwrite_f(dimer_dset_id(7),h5t_native_double,nucvec,dimer_init_dims_nucvec,err)

else

 dimer_dimer_extend = (/numdims,2,int(cstep/Freq4)+1/)
 dimer_dimer_extend_offset = (/0,0,int(cstep/Freq4)/)
 dimer_dimer_extend_count = (/numdims,2,1/)
 dimer_nucvec_extend = (/numdims,3,int(cstep/Freq4)+1/)
 dimer_nucvec_extend_offset = (/0,0,int(cstep/Freq4)/)
 dimer_nucvec_extend_count = (/numdims,3,1/)

 do i = 1,7,1
  if ( i /= 7 ) then
   call h5dset_extent_f(dimer_dset_id(i),dimer_dimer_extend,err)
   call h5screate_simple_f(rank(3),dimer_dimer_extend_count,dimer_memspace,err)
   call h5dget_space_f(dimer_dset_id(i),dimer_dspace_id(i),err)
   call h5sselect_hyperslab_f(dimer_dspace_id(i),h5s_select_set_f,dimer_dimer_extend_offset,dimer_dimer_extend_count,err)
  else if ( i == 7 ) then
   call h5dset_extent_f(dimer_dset_id(i),dimer_nucvec_extend,err)
   call h5screate_simple_f(rank(3),dimer_nucvec_extend_count,dimer_memspace,err)
   call h5dget_space_f(dimer_dset_id(i),dimer_dspace_id(i),err)
   call h5sselect_hyperslab_f(dimer_dspace_id(i),h5s_select_set_f,dimer_nucvec_extend_offset,dimer_nucvec_extend_count,err)
  end if
  if ( i == 1 ) then
   call h5dwrite_f(dimer_dset_id(i),h5t_native_double,xdim,dimer_dimer_extend_count,err,dimer_memspace,dimer_dspace_id(i))
  else if ( i == 2 ) then
   call h5dwrite_f(dimer_dset_id(i),h5t_native_double,ydim,dimer_dimer_extend_count,err,dimer_memspace,dimer_dspace_id(i))
  else if ( i == 3 ) then
   call h5dwrite_f(dimer_dset_id(i),h5t_native_double,zdim,dimer_dimer_extend_count,err,dimer_memspace,dimer_dspace_id(i))
  else if ( i == 4 ) then
   call h5dwrite_f(dimer_dset_id(i),h5t_native_double,vxdim,dimer_dimer_extend_count,err,dimer_memspace,dimer_dspace_id(i))
  else if ( i == 5 ) then
   call h5dwrite_f(dimer_dset_id(i),h5t_native_double,vydim,dimer_dimer_extend_count,err,dimer_memspace,dimer_dspace_id(i))
  else if ( i == 6 ) then
   call h5dwrite_f(dimer_dset_id(i),h5t_native_double,vzdim,dimer_dimer_extend_count,err,dimer_memspace,dimer_dspace_id(i))
  else if ( i == 7 ) then
   call h5dwrite_f(dimer_dset_id(i),h5t_native_double,nucvec,dimer_nucvec_extend_count,err,dimer_memspace,dimer_dspace_id(i))
  end if
  call h5sclose_f(dimer_memspace,err)
 end do

end if

end if

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_solv_type_num(NT,Vol,flag,iTotA,iTotB,iTotS)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
integer(kind=dp) , intent(in) :: Vol                            ! Volume of simulation box
integer(kind=sp) , dimension(:), intent(in) :: flag             ! Identity of solvent
integer(kind=dp) , intent(out) :: iTotA                         ! Total number of A particles
integer(kind=dp) , intent(out) :: iTotB                         ! Total number of B particles
integer(kind=dp) , intent(out) :: iTotS                         ! Total number of S particles

! Definitions for internal variables
integer(kind=dp) :: i                                           ! Loop variable

iTotA=0.
iTotB=0.
iTotS=0.

do i = 1,NT
 if ( flag(i) == 0 ) then
  iTotA = iTotA + 1
 else if ( flag(i) == 1 ) then
  iTotB = iTotB + 1
 else if ( flag(i) == 2 ) then
  iTotS = iTotS + 1
 end if
end do

end subroutine calc_solv_type_num

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module
