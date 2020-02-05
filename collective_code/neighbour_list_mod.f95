module neighbour_list_mod
implicit none

contains

subroutine neighbour_list(NT,Box,cutoff,skinrad,Xsol,Ysol,Zsol,Xdim,Ydim,Zdim,neighbour_num,neighbour)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
real(kind=dp) , dimension(:) , intent(in) :: Box                ! Length of simulation box
!integer , intent(in) :: sphere                                  ! Incoming sphere identity
real(kind=dp) , intent(in) :: cutoff                            ! LJ cutoff for monomer
real(kind=dp) , intent(in) :: skinrad                           ! Neighbour list skin radius
real(kind=dp) , dimension(:) , intent(in) :: Xsol               ! Force on solvent in x direction
real(kind=dp) , dimension(:) , intent(in) :: Ysol               ! Force on solvent in y direction
real(kind=dp) , dimension(:) , intent(in) :: Zsol               ! Force on solvent in z direction
real(kind=dp) , intent(in) :: Xdim             ! Force on solvent in x direction
real(kind=dp) , intent(in) :: Ydim             ! Force on solvent in y direction
real(kind=dp) , intent(in) :: Zdim             ! Force on solvent in z direction
integer(kind=dp) , intent(out) :: neighbour_num! Number of solvent particles in list
integer(kind=dp) , dimension(:) , intent(out) :: neighbour  ! Identify solvent number in list

! Definitions for internal variables
integer(kind=dp) :: i,j                                         ! Loop variables
real(kind=dp) , dimension(3) :: dR                              ! Solvent-monomer distance
real(kind=dp) :: full_rad
real(kind=dp) :: full_rad2

full_rad = cutoff + skinrad
full_rad2 = (cutoff+skinrad)**2

! Check to see if within skin distance and if so, add to neighbour list
do i = 1,NT

 dR(1) = Xsol(i) - Xdim
 dR(2) = Ysol(i) - Ydim
 dR(3) = Zsol(i) - Zdim

 do j = 1,2
  if ( dR(j) > 0.5d0*Box(j) ) then
   dR(j) = dR(j) - Box(j)
  else if ( dR(j) < -0.5d0*Box(j) ) then
   dR(j) = dR(j) + Box(j)
  end if
 end do 
 
 if ( dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) <= full_rad2 ) then
  neighbour_num = neighbour_num + 1
  neighbour(neighbour_num) = i
 end if

end do

end subroutine neighbour_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inout_particle_list(NT,box,cutoffC,cutoffN,skinrad,xsol,ysol,zsol,xdim,ydim,zdim,&
outside_particle_num,outside_particle,inside_particle_num,inside_particle,numdims)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
real(kind=dp) , dimension(:) , intent(in) :: box                ! Length of simulation box
real(kind=dp) , intent(in) :: cutoffC                           ! LJ cutoff for monomer
real(kind=dp) , intent(in) :: cutoffN                           ! LJ cutoff for monomer
real(kind=dp) , intent(in) :: skinrad                           ! Neighbour list skin radius
real(kind=dp) , dimension(:) , intent(in) :: xsol               ! Force on solvent in x direction
real(kind=dp) , dimension(:) , intent(in) :: ysol               ! Force on solvent in y direction
real(kind=dp) , dimension(:) , intent(in) :: zsol               ! Force on solvent in z direction
real(kind=dp) , dimension(:,:) , intent(in) :: xdim               ! Force on solvent in x direction
real(kind=dp) , dimension(:,:) , intent(in) :: ydim               ! Force on solvent in y direction
real(kind=dp) , dimension(:,:) , intent(in) :: zdim               ! Force on solvent in z direction
integer(kind=dp) , intent(inout) :: outside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: outside_particle
integer(kind=dp) , intent(inout) :: inside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: inside_particle
integer , intent(in) :: numdims                                              ! Number of dimers

! Definitions for internal variables
integer(kind=dp) :: i,j,k,m                                       ! Loop variables
real(kind=dp) , dimension(numdims,3,2) :: dR                      ! Solvent-monomer distance
integer :: bulk
integer(kind=dp) , dimension(:) , allocatable :: outside_particle_num_part
integer(kind=dp) , dimension(:,:) , allocatable :: outside_particle_part
integer(kind=dp) , dimension(:) , allocatable :: inside_particle_num_part
integer(kind=dp) , dimension(:,:) , allocatable :: inside_particle_part
integer(kind=dp) :: thread_num,omp_get_thread_num
integer(kind=dp) :: num_threads,omp_get_num_threads
integer , dimension(8) :: values
real(kind=dp) :: full_radC
real(kind=dp) :: full_radC2
real(kind=dp) :: full_radN
real(kind=dp) :: full_radN2

full_radN = cutoffN + skinrad
full_radN2 = (cutoffN+skinrad)**2
full_radC = cutoffC + skinrad
full_radC2 = (cutoffC+skinrad)**2

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED(num_threads)
num_threads = omp_get_num_threads()
!$OMP END PARALLEL

allocate(outside_particle_num_part(num_threads))
allocate(outside_particle_part(num_threads,NT))
allocate(inside_particle_num_part(num_threads))
allocate(inside_particle_part(num_threads,NT))

outside_particle_num_part = 0
inside_particle_num_part = 0

! Check to see if within skin distance and if so, add to neighbour list
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(NT,numdims,xsol,ysol,zsol,xdim,ydim,zdim,cutoffN,cutoffC,skinrad,box) &
!$OMP SHARED(outside_particle_num,outside_particle,inside_particle_num,inside_particle) &
!$OMP PRIVATE(i,j,m,k,dR,bulk,thread_num,num_threads)  &
!$OMP SHARED(outside_particle_num_part,outside_particle_part,inside_particle_num_part,inside_particle_part) 
do i = 1,NT

 thread_num = omp_get_thread_num()+1

 do m = 1,numdims
  do j = 1,2
   dR(m,1,j) = Xsol(i) - Xdim(m,j)
   dR(m,2,j) = Ysol(i) - Ydim(m,j)
   dR(m,3,j) = Zsol(i) - Zdim(m,j)
  end do

  do j = 1,2
   do k = 1,2
    if ( dR(m,j,k) > 0.5d0*Box(j) ) then
     dR(m,j,k) = dR(m,j,k) - Box(j)
    else if ( dR(m,j,k) < -0.5d0*Box(j) ) then
     dR(m,j,k) = dR(m,j,k) + Box(j)
    end if
   end do 
  end do
 end do

 bulk = 0
 do m = 1,numdims
  if ( dR(m,1,1)*dR(m,1,1) + dR(m,2,1)*dR(m,2,1) + dR(m,3,1)*dR(m,3,1) <= (cutoffC+skinrad)**2 .or. &
      dR(m,1,2)*dR(m,1,2) + dR(m,2,2)*dR(m,2,2) + dR(m,3,2)*dR(m,3,2) <= (cutoffN+skinrad)**2 ) then
   bulk = 1
  end if
 end do
 if ( bulk == 0 ) then
  outside_particle_num_part(thread_num) = outside_particle_num_part(thread_num) + 1
  outside_particle_part(thread_num,outside_particle_num_part(thread_num)) = i
 end if
 if ( bulk == 1 ) then
  inside_particle_num_part(thread_num) = inside_particle_num_part(thread_num) + 1
  inside_particle_part(thread_num,inside_particle_num_part(thread_num)) = i
 end if

end do
!$OMP END PARALLEL DO

do i = 1,num_threads
 outside_particle((outside_particle_num+1):(outside_particle_num+outside_particle_num_part(i))) = &
 outside_particle_part(i,1:outside_particle_num_part(i))
 outside_particle_num = outside_particle_num + outside_particle_num_part(i)
 inside_particle((inside_particle_num+1):(inside_particle_num+inside_particle_num_part(i))) = &
 inside_particle_part(i,1:inside_particle_num_part(i))
 inside_particle_num = inside_particle_num + inside_particle_num_part(i)
end do

deallocate(outside_particle_num_part)
deallocate(outside_particle_part)
deallocate(inside_particle_num_part)
deallocate(inside_particle_part)

end subroutine inout_particle_list

end module neighbour_list_mod
