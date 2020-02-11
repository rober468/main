module output_mod
use hdf5
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)             ! single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! double precision

! global HDF5 variables
 ! overall file
character(len=7), parameter :: datafile = "data.h5"
integer(hid_t) :: data_id
integer(kind=sp) :: err
integer(kind=sp), dimension(3), parameter :: rank = (/1,2,3/)
! group initiation variables
integer(hid_t), dimension(2) :: sub_group_id
! angle subgroup
integer(hsize_t), dimension(2) :: ang_init_dims
integer(hid_t), dimension(1) :: ang_dspace_id
integer(hid_t), dimension(1) :: ang_dset_id
! restart subgroup
integer(hsize_t) , dimension(1) :: res_init_numrotors_dims
integer(hsize_t) , dimension(1), parameter :: res_init_1_dims = 1
integer(hid_t), dimension(3) :: res_dspace_id
integer(hid_t), dimension(3) :: res_dset_id

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_hdf5_initial(rotor_CM,timestep,kBT,Rc,Rn,d_CN,sep,dens,Ean,Ebn,D,D_R,k1,k_1,kappa,&
                               derj_len_N,p_plus,p_minus,Z,angle,F,w,numrotors,freq1,freq2)
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)             ! single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! double precision

! external variables
real(kind=dp), dimension(:,:), intent(in) :: rotor_CM   ! rotor centre-of-mass coordinates
real(kind=dp), intent(in) :: timestep                   ! timestep
real(kind=dp), intent(in) :: kBT                        ! temperature
real(kind=dp), intent(in) :: Rc                         ! radius of catalytic sphere
real(kind=dp), intent(in) :: Rn                         ! radius of noncatalytic sphere
real(kind=dp), intent(in) :: d_CN                       ! distance between C and N sphere for each rotor
real(kind=dp), intent(in) :: sep                        ! separation between rotors
real(kind=dp), intent(in) :: dens                       ! total number density of solvent
real(kind=dp), intent(in) :: Ean                        ! LJ potential of A with N sphere
real(kind=dp), intent(in) :: Ebn                        ! LJ potential of B with N sphere
real(kind=dp), intent(in) :: D                          ! diffusion coefficient of solvent
real(kind=dp), intent(in) :: D_R                        ! rotation diffusion coefficient of rotors
real(kind=dp), intent(in) :: k1                         ! solvent forward reaction rate A -> B
real(kind=dp), intent(in) :: k_1                        ! solvent reverse reaction rate B -> A
real(kind=dp), intent(in) :: kappa                      ! reaction exp. decay constant
real(kind=dp), intent(in) :: derj_len_N                 ! Derjaguin length for N sphere
real(kind=dp), intent(in) :: p_plus                     ! probability of forward reaction at motor
real(kind=dp), intent(in) :: p_minus                    ! probability of backward reaction at motor
real(kind=dp), intent(in) :: Z
real(kind=dp), dimension(:), intent(in) :: angle        ! angle of each rotor
real(kind=dp), dimension(:), intent(in) :: F            ! "force" on N sphere; contains sin_theta for torque
real(kind=dp), dimension(:), intent(in) :: w            ! standard normal random variable
integer, intent(in) :: numrotors                        ! total number of rotors
integer, intent(in) :: freq1                            ! frequency of angle recording
integer, intent(in) :: freq2                            ! frequency for restart data

! internal variables
integer :: i                                                    ! loop variables
integer(hsize_t), dimension(1) :: maxdims_1
integer(hsize_t), dimension(2) :: maxdims_2
! data file attributes
integer(hid_t), dimension(20) :: data_attr_id
integer(hsize_t), dimension(2) :: data_attr_rotor_CM_dims
integer(hsize_t), dimension(1), parameter :: data_attr_1_dims = 1
integer(hsize_t), dimension(1), parameter :: data_attr_2_dims = 2
integer(hid_t), dimension(20) :: data_attr_dspace_id
character(len=19), dimension(20), parameter :: data_attr_names = &
(/"rotor_CM           ","timestep           ","kB*temperature     ","C_sphere_radius    ",&
"N_sphere_radius    ","d_CN               ","separation         ","number_density     ",&
"LJ_energy_N_A      ","LJ_energy_N_B      ","solv_MPC_diffusion ","rotor_rot_diff     ",&
"solv_rxn_rate_k1   ","solv_rxn_rate_k_1  ","kappa              ","Derjaguin_length   ",&
"prob_forw_rotor_rxn","prob_back_rotor_rxn","Z                  ","number_rotors      "/)
! group initiation variables
character(len=7), dimension(2), parameter :: sub_group_names = (/"angle  ","restart"/)
! angle subgroup
integer(hid_t) :: ang_attr_id
integer(hsize_t), dimension(1), parameter :: ang_attr_1_dims = 1
integer(hsize_t), dimension(1), parameter :: ang_attr_2_dims = 2
integer(hid_t) :: ang_attr_dspace_id
character(len=12), parameter :: ang_attr_names = "io_frequency"
integer(hid_t) :: ang_crp_list
integer(hsize_t), dimension(2) :: ang_chunk
character(len=5), dimension(1), parameter :: ang_dset_names = (/"angle"/)
! restart subgroup
integer(hid_t) :: res_attr_id
integer(hsize_t), dimension(1), parameter :: res_attr_1_dims = 1
integer(hid_t) :: res_attr_dspace_id
character(len=12), parameter :: res_attr_names = "io_frequency"
integer(hid_t), dimension(1) :: res_crp_list
integer(hsize_t), dimension(1), parameter :: res_chunk_1 = (/75000/)
character(len=12), dimension(3), parameter :: res_dset_names = (/"force       ","st_norm_rand",&
"cstep       "/)

! open HDF5
call h5open_f(err)

! derived quantities
maxdims_1 = (/h5s_unlimited_f/)
maxdims_2 = (/h5s_unlimited_f,h5s_unlimited_f/)
data_attr_rotor_CM_dims = (/numrotors,2/)
ang_init_dims = (/numrotors,1/)
ang_chunk = (/numrotors,int(floor(100000./numrotors))/)
res_init_numrotors_dims = (/numrotors/)

! create file
call h5fcreate_f(datafile,h5f_acc_trunc_f,data_id,err)

! general simulation attributes
 ! 1d quantities - double precision
do i = 1,20,1

 if ( i == 1 ) then ! 2d quantities - double precision

  call h5screate_simple_f(rank(2),data_attr_rotor_CM_dims,data_attr_dspace_id(i),err)
  call h5acreate_f(data_id,data_attr_names(i),h5t_native_double,data_attr_dspace_id(i),data_attr_id(i),err)
  call h5awrite_f(data_attr_id(i),h5t_native_double,rotor_CM,data_attr_2_dims,err)

 else if ( i >= 2 .and. i <= 19 ) then ! 1d quantities - double precision

  call h5screate_simple_f(rank(1),data_attr_1_dims,data_attr_dspace_id(i),err)
  call h5acreate_f(data_id,data_attr_names(i),h5t_native_double,data_attr_dspace_id(i),data_attr_id(i),err)
  if ( i == 2 ) then             ! timestep
   call h5awrite_f(data_attr_id(i),h5t_native_double,timestep,data_attr_1_dims,err)
  else if ( i == 3 ) then        ! kBT
   call h5awrite_f(data_attr_id(i),h5t_native_double,kBT,data_attr_1_dims,err)
  else if ( i == 4 ) then        ! Rc
   call h5awrite_f(data_attr_id(i),h5t_native_double,Rc,data_attr_1_dims,err)
  else if ( i == 5 ) then        ! Rn
   call h5awrite_f(data_attr_id(i),h5t_native_double,Rn,data_attr_1_dims,err)
  else if ( i == 6 ) then        ! d_CN
   call h5awrite_f(data_attr_id(i),h5t_native_double,d_CN,data_attr_1_dims,err)
  else if ( i == 7 ) then        ! sep
   call h5awrite_f(data_attr_id(i),h5t_native_double,sep,data_attr_1_dims,err)
  else if ( i == 8 ) then        ! dens
   call h5awrite_f(data_attr_id(i),h5t_native_double,dens,data_attr_1_dims,err)
  else if ( i == 9 ) then        ! Ean
   call h5awrite_f(data_attr_id(i),h5t_native_double,Ean,data_attr_1_dims,err)
  else if ( i == 10 ) then       ! Ebn
   call h5awrite_f(data_attr_id(i),h5t_native_double,Ebn,data_attr_1_dims,err)
  else if ( i == 11 ) then       ! D
   call h5awrite_f(data_attr_id(i),h5t_native_double,D,data_attr_1_dims,err)
  else if ( i == 12 ) then       ! D_R
   call h5awrite_f(data_attr_id(i),h5t_native_double,D_R,data_attr_1_dims,err)
  else if ( i == 13 ) then       ! k1
   call h5awrite_f(data_attr_id(i),h5t_native_double,k1,data_attr_1_dims,err)
  else if ( i == 14 ) then       ! k_1
   call h5awrite_f(data_attr_id(i),h5t_native_double,k_1,data_attr_1_dims,err)
  else if ( i == 15 ) then       ! kappa
   call h5awrite_f(data_attr_id(i),h5t_native_double,kappa,data_attr_1_dims,err)
  else if ( i == 16 ) then       ! derj_len_N
   call h5awrite_f(data_attr_id(i),h5t_native_double,derj_len_N,data_attr_1_dims,err)
  else if ( i == 17 ) then       ! p_plus
   call h5awrite_f(data_attr_id(i),h5t_native_double,p_plus,data_attr_1_dims,err)
  else if ( i == 18 ) then       ! p_minus
   call h5awrite_f(data_attr_id(i),h5t_native_double,p_minus,data_attr_1_dims,err)
  else if ( i == 19 ) then       ! Z
   call h5awrite_f(data_attr_id(i),h5t_native_double,Z,data_attr_1_dims,err)
  end if

 else if ( i == 20 ) then ! 1d quantities - integer

  call h5screate_simple_f(rank(1),data_attr_1_dims,data_attr_dspace_id(i),err)
  call h5acreate_f(data_id,data_attr_names(i),h5t_native_integer,data_attr_dspace_id(i),data_attr_id(i),err)
  call h5awrite_f(data_attr_id(i),h5t_native_integer,numrotors,data_attr_1_dims,err)

 end if

 call h5aclose_f(data_attr_id(i),err)
end do


! create subgroups
do i = 1,2,1
 call h5gcreate_f(data_id,sub_group_names(i),sub_group_id(i),err)
end do

! create dataspaces, dataspace attributes, and properties
! angle subgroup
 ! attributes
call h5screate_simple_f(rank(1),ang_attr_1_dims,ang_attr_dspace_id,err)
call h5acreate_f(sub_group_id(1),ang_attr_names,h5t_native_integer,ang_attr_dspace_id,ang_attr_id,err)
call h5awrite_f(ang_attr_id,h5t_native_integer,freq1,ang_attr_1_dims,err)
call h5aclose_f(ang_attr_id,err)
 ! dataspace creation, properties, chunking, dataset creation
call h5screate_simple_f(rank(2),ang_init_dims,ang_dspace_id(1),err,maxdims_2)
call h5pcreate_f(h5p_dataset_create_f,ang_crp_list,err)
call h5pset_chunk_f(ang_crp_list,rank(2),ang_chunk,err)
call h5dcreate_f(sub_group_id(1),ang_dset_names(1),h5t_native_double,ang_dspace_id(1),ang_dset_id(1),err,ang_crp_list)

! restart subgroup
 ! attributes
call h5screate_simple_f(rank(1),res_attr_1_dims,res_attr_dspace_id,err)
call h5acreate_f(sub_group_id(2),res_attr_names,h5t_native_integer,res_attr_dspace_id,res_attr_id,err)
call h5awrite_f(res_attr_id,h5t_native_integer,freq2,res_attr_1_dims,err)
call h5aclose_f(res_attr_id,err)
 ! dataspace creation, properties, chunking, dataset creation
do i = 1,3,1
 if ( i >= 1 .and. i <= 2 ) then ! double precision
  call h5screate_simple_f(rank(1),res_init_numrotors_dims,res_dspace_id(i),err,maxdims_1)
  call h5pcreate_f(h5p_dataset_create_f,res_crp_list(i),err)
  call h5pset_chunk_f(res_crp_list(i),rank(1),res_chunk_1,err)
  call h5dcreate_f(sub_group_id(2),res_dset_names(i),h5t_native_double,res_dspace_id(i),res_dset_id(i),err,res_crp_list(i))
 else if ( i == 3 ) then ! integer
  call h5screate_simple_f(rank(1),res_init_1_dims,res_dspace_id(i),err,maxdims_1)
  call h5pcreate_f(h5p_dataset_create_f,res_crp_list(i),err)
  call h5pset_chunk_f(res_crp_list(i),rank(1),res_chunk_1,err)
  call h5dcreate_f(sub_group_id(2),res_dset_names(i),h5t_native_integer,res_dspace_id(i),res_dset_id(i),err,res_crp_list(i))
 end if
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_hdf5_close
implicit none

! internal variable
integer :: i

! angle subgroup
call h5dclose_f(ang_dset_id(1),err)
call h5sclose_f(ang_dspace_id(1),err)

! restart subgroup
do i = 1,3,1
 call h5dclose_f(res_dset_id(i),err)
 call h5sclose_f(res_dspace_id(i),err)
end do

! close subgroups
do i = 1,2,1
 call h5gclose_f(sub_group_id(i),err)
end do

! close file
call h5fclose_f(data_id,err)

call h5close_f(err)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_hdf5_write(cstep,numrotors,angle,F,w,freq1,freq2)
implicit none

! external variables
integer, intent(in) :: cstep
integer, intent(in) :: numrotors
real(kind=dp), dimension(:), intent(in) :: angle
real(kind=dp), dimension(:), intent(in) :: F
real(kind=dp), dimension(:), intent(in) :: w
integer, intent(in) :: freq1
integer, intent(in) :: freq2

! internal variables
integer(hsize_t), dimension(2) :: ang_extend
integer(hsize_t), dimension(2) :: ang_extend_offset
integer(hsize_t), dimension(2) :: ang_extend_count
integer(hid_t) :: ang_memspace

! initial or update angle
if ( mod(cstep,freq1) == 0 ) then
 if ( cstep == 0 ) then

  ang_init_dims = (/numrotors,1/)
  call h5dwrite_f(ang_dset_id(1),h5t_native_double,angle,ang_init_dims,err)

 else

  ang_extend = (/numrotors,int(cstep/freq1)+1/)
  ang_extend_offset = (/0,int(cstep/freq1)/)
  ang_extend_count = (/numrotors,1/)
  
  call h5dset_extent_f(ang_dset_id(1),ang_extend,err)
  call h5screate_simple_f(rank(2),ang_extend_count,ang_memspace,err)
  call h5dget_space_f(ang_dset_id(1),ang_dspace_id(1),err)
  call h5sselect_hyperslab_f(ang_dspace_id(1),h5s_select_set_f,ang_extend_offset,ang_extend_count,err)
  call h5dwrite_f(ang_dset_id(1),h5t_native_double,angle,ang_extend_count,err,ang_memspace,ang_dspace_id(1))

 end if
end if

! update F, w, timestep
if ( mod(cstep,freq2) == 0 ) then
  res_init_numrotors_dims = (/numrotors/)
  call h5dwrite_f(res_dset_id(1),h5t_native_double,F,res_init_numrotors_dims,err)
  call h5dwrite_f(res_dset_id(2),h5t_native_double,w,res_init_numrotors_dims,err)
  call h5dwrite_f(res_dset_id(3),h5t_native_integer,cstep,res_init_1_dims,err)

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module output_mod
