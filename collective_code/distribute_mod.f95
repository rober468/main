module distribute_mod
implicit none

contains

subroutine distribute(NT,numdims,cutoff,Rc,Rn,flag,Xsol,Ysol,Zsol,Xdim,Ydim,Zdim,Box,tACellN,tBCellN,tSCellN, &
                      sACellN,sBCellN,sSCellN,cutoffC2,cutoffN2,cstep,outside_particle_num,outside_particle,&
                      inside_particle_num,inside_particle)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=sp) , intent(in) :: NT                             ! Total number of particles
integer , intent(in) :: numdims                                 ! Number of dimers
integer(kind=sp) , dimension(:), intent(in) :: flag             ! Identity of solvent
real(kind=dp) , dimension(:) , intent(in) :: Xsol               ! Solvent positions, x
real(kind=dp) , dimension(:) , intent(in) :: Ysol               ! Solvent positions, y
real(kind=dp) , dimension(:) , intent(in) :: Zsol               ! Solvent positions, z
real(kind=dp) , dimension(:,:) , intent(in) :: Xdim             ! Solvent positions, x
real(kind=dp) , dimension(:,:) , intent(in) :: Ydim             ! Solvent positions, y
real(kind=dp) , dimension(:,:) , intent(in) :: Zdim             ! Solvent positions, z
real(kind=dp) , dimension(:) , intent(in) :: Box                ! Length of simulation box;x=1,y=2,z=3
integer(kind=dp) , dimension(:), intent(out) :: tACellN         ! Total number of A solvent in cell N
integer(kind=dp) , dimension(:), intent(out) :: tBCellN         ! Total number of B solvent in cell N
integer(kind=dp) , dimension(:), intent(out) :: tSCellN         ! Total number of S solvent in cell N
integer(kind=dp) , dimension(:,:), intent(out) :: sACellN       ! Label A in Cell N with i
integer(kind=dp) , dimension(:,:), intent(out) :: sBCellN       ! Label B in Cell N with i
integer(kind=dp) , dimension(:,:), intent(out) :: sSCellN       ! Label S in Cell N with i
real(kind=dp) , intent(in) :: cutoff                            ! Cut. rad. C sphere LJ pot., squared
real(kind=dp) , intent(in) :: Rc                                ! Radius C sphere
real(kind=dp) , intent(in) :: Rn                                ! Radius N sphere
real(kind=dp) , intent(in) :: cutoffC2                          ! Cut. C sphere LJ pot. squared
real(kind=dp) , intent(in) :: cutoffN2                          ! Cut. N sphere LJ pot. squared
integer(kind=dp) , intent(inout) :: outside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: outside_particle
integer(kind=dp) , intent(inout) :: inside_particle_num
integer(kind=dp) , dimension(:) , intent(inout) :: inside_particle

! Definitions for internal variables
real(kind=dp) :: dRsq                                           ! Solvent-monomer distance, squared
real(kind=dp) , dimension(3,numdims,2) :: dR                    ! Solvent-monomer distance
integer(kind=dp) :: i,j,k,l                                     ! Loop variable
integer(kind=dp) :: CellN                                       ! Specify the cell in the sim. box
integer :: bulk                                                 ! If outside of monomer int. zone, = 0
integer :: cstep
integer , dimension(8) :: values
integer , dimension(:,:) , allocatable :: tacelln_part
integer , dimension(:,:) , allocatable :: tbcelln_part
integer , dimension(:,:) , allocatable :: tscelln_part
integer , dimension(:,:,:) , allocatable :: sacelln_part
integer , dimension(:,:,:) , allocatable :: sbcelln_part
integer , dimension(:,:,:) , allocatable :: sscelln_part
integer(kind=dp) :: thread_num,omp_get_thread_num
integer(kind=dp) :: num_threads,omp_get_num_threads
integer :: size_t
integer :: size_m

size_t = size(tACellN)
size_m = size(sACellN(1,:))

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED(num_threads)
num_threads = omp_get_num_threads()
!$OMP END PARALLEL

allocate(tacelln_part(num_threads,size_t))
allocate(tbcelln_part(num_threads,size_t))
allocate(tscelln_part(num_threads,size_t))
allocate(sacelln_part(num_threads,size_t,size_m))
allocate(sbcelln_part(num_threads,size_t,size_m))
allocate(sscelln_part(num_threads,size_t,size_m))

!!!!!!Distribute the molecules to the cells
!Clear previous values and start with all values zero
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(size_t,tacelln,tbcelln,tscelln,tacelln_part,tbcelln_part,tscelln_part) &
!$OMP PRIVATE(i)
do i = 1,size_t
 tACellN(i) = 0
 tacelln_part(:,i) = 0
 tBCellN(i) = 0
 tbcelln_part(:,i) = 0
 tSCellN(i) = 0
 tscelln_part(:,i) = 0
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(outside_particle_num,flag,tacelln_part,tbcelln_part,tscelln_part,sacelln_part,sbcelln_part,sscelln_part) &
!$OMP SHARED(xsol,ysol,zsol,outside_particle,box) &
!$OMP PRIVATE(i,celln,thread_num)
do i = 1,outside_particle_num
 
 thread_num = omp_get_thread_num()+1

 if ( flag(outside_particle(i)) == 0 ) then
  CellN = dint( floor(Xsol(outside_particle(i)))+Box(1)*floor(Ysol(outside_particle(i)))&
          +Box(1)*Box(2)*floor(Zsol(outside_particle(i)))+1 )  ! Identify cell
  tACellN_part(thread_num,CellN) = tACellN_part(thread_num,CellN) + 1             ! Total number of A in cell N
  sACellN_part(thread_num,CellN,tACellN_part(thread_num,CellN)) = outside_particle(i)                 ! Specify/label newly added A with i

 else if ( flag(outside_particle(i)) == 1 ) then
   CellN = dint( floor(Xsol(outside_particle(i)))+Box(1)*floor(Ysol(outside_particle(i)))&
           +Box(1)*Box(2)*floor(Zsol(outside_particle(i)))+1 )  ! Identify cell
   tBCellN_part(thread_num,CellN) = tBCellN_part(thread_num,CellN) + 1         ! Total number of B in cell N
   sBCellN_part(thread_num,CellN,tBCellN_part(thread_num,CellN)) = outside_particle(i) ! Specify/label newly added B with i

 else if ( flag(outside_particle(i)) == 2 ) then
   CellN = dint(floor(Xsol(outside_particle(i)))+Box(1)*floor(Ysol(outside_particle(i)))&
           +Box(1)*Box(2)*floor(Zsol(outside_particle(i)))+1 )  ! Identify cell
   tSCellN_part(thread_num,CellN) = tSCellN_part(thread_num,CellN) + 1         !Total number of S in cell N
   sSCellN_part(thread_num,CellN,tSCellN_part(thread_num,CellN)) = outside_particle(i) ! Specify/label newly added S with i

 end if

end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(inside_particle_num,numdims,xsol,ysol,zsol,xdim,ydim,zdim,cutoffC2,cutoffN2) &
!$OMP SHARED(flag,tacelln_part,tbcelln_part,tscelln_part,sacelln_part,sbcelln_part,sscelln_part,inside_particle,box) &
!$OMP PRIVATE(i,k,j,l,celln,bulk,dR,thread_num)
do i = 1,inside_particle_num

 thread_num = omp_get_thread_num()+1

 do k = 1,numdims
  do j = 1,2
   dR(1,k,j) = Xsol(inside_particle(i)) - Xdim(k,j)
   dR(2,k,j) = Ysol(inside_particle(i)) - Ydim(k,j)
   dR(3,k,j) = Zsol(inside_particle(i)) - Zdim(k,j)
  end do
 end do

 do j = 1,2
  do l = 1,numdims
   do k =1,2
    if ( dR(j,l,k) > 0.5d0*Box(j) ) then
     dR(j,l,k) = dR(j,l,k) - Box(j)
    else if ( dR(j,l,k) < -0.5d0*Box(j) ) then
     dR(j,l,k) = dR(j,l,k) + Box(j)
    end if
   end do
  end do
 end do

 bulk = 1
 do k = 1,numdims
  if ( dR(1,k,1)*dR(1,k,1) + dR(2,k,1)*dR(2,k,1) + dR(3,k,1)*dR(3,k,1) <=  cutoffC2 &
      .or. dR(1,k,2)*dR(1,k,2) + dR(2,k,2)*dR(2,k,2) + dR(3,k,2)*dR(3,k,2) <= cutoffN2 ) then
 bulk = 0
  end if
 end do

 if ( bulk == 1 ) then

 if ( flag(inside_particle(i)) == 0 ) then
  CellN = dint( floor(Xsol(inside_particle(i)))+Box(1)*floor(Ysol(inside_particle(i)))&
          +Box(1)*Box(2)*floor(Zsol(inside_particle(i)))+1 )  ! Identify cell
  tACellN_part(thread_num,CellN) = tACellN_part(thread_num,CellN) + 1 ! Total number of A in cell N
  sACellN_part(thread_num,CellN,tACellN_part(thread_num,CellN)) = inside_particle(i) ! Specify/label newly added A with i

 else if ( flag(inside_particle(i)) == 1 ) then
   CellN = dint( floor(Xsol(inside_particle(i)))+Box(1)*floor(Ysol(inside_particle(i)))&
           +Box(1)*Box(2)*floor(Zsol(inside_particle(i)))+1 )  ! Identify cell
   tBCellN_part(thread_num,CellN) = tBCellN_part(thread_num,CellN) + 1 ! Total number of B in cell N
   sBCellN_part(thread_num,CellN,tBCellN_part(thread_num,CellN)) = inside_particle(i) ! Specify/label newly added B with i

 else if ( flag(inside_particle(i)) == 2 ) then
   CellN = dint(floor(Xsol(inside_particle(i)))+Box(1)*floor(Ysol(inside_particle(i)))&
           +Box(1)*Box(2)*floor(Zsol(inside_particle(i)))+1 )  ! Identify cell
   tSCellN_part(thread_num,CellN) = tSCellN_part(thread_num,CellN) + 1 ! Totalnumber of B in cell N
   sSCellN_part(thread_num,CellN,tSCellN_part(thread_num,CellN)) = inside_particle(i) ! Specify/label newly added B with i

 end if

 end if

end do
!$OMP END PARALLEL DO

do i = 1,num_threads

 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(i,size_t,tacelln,tacelln_part,tbcelln,tbcelln_part,tscelln,tscelln_part) &
 !$OMP SHARED(sacelln,sacelln_part,sbcelln,sbcelln_part,sscelln,sscelln_part) &
 !$OMP PRIVATE(j)
 do j = 1,size_t
  sacelln(j,(tacelln(j)+1):(tacelln(j)+tacelln_part(i,j))) = sacelln_part(i,j,1:tacelln_part(i,j))
  tacelln(j) = tacelln(j) + tacelln_part(i,j) 
  sbcelln(j,(tbcelln(j)+1):(tbcelln(j)+tbcelln_part(i,j))) = sbcelln_part(i,j,1:tbcelln_part(i,j))
  tbcelln(j) = tbcelln(j) + tbcelln_part(i,j) 
  sscelln(j,(tscelln(j)+1):(tscelln(j)+tscelln_part(i,j))) = sscelln_part(i,j,1:tscelln_part(i,j))
  tscelln(j) = tscelln(j) + tscelln_part(i,j) 
 end do 
 !$OMP END PARALLEL DO

end do

deallocate(tacelln_part,tbcelln_part,tscelln_part,sacelln_part,sbcelln_part,sscelln_part)

end subroutine distribute

end module distribute_mod
