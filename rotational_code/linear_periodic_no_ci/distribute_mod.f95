module distribute_mod
implicit none

contains

subroutine distribute(NT,numdims,cutoff,Rc,Rn,flag,Xsol,Ysol,Zsol,Xdim,Ydim,Zdim,Box,tACellN,tBCellN, &
                      sACellN,sBCellN,cutoffC2,cutoffN2,cstep,outside_particle_num,outside_particle,&
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
integer(kind=dp) , dimension(:), intent(out) :: tACellN         ! Total number of F solvent in cell N
integer(kind=dp) , dimension(:), intent(out) :: tBCellN         ! Total number of I solvent in cell N
integer(kind=dp) , dimension(:,:), intent(out) :: sACellN       ! Label F in Cell N with i
integer(kind=dp) , dimension(:,:), intent(out) :: sBCellN       ! Label I in Cell N with i
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

!!!!!!Distribute the molecules to the cells
!Clear previous values and start with all values zero
sACellN = 0
sBCellN = 0
tACellN = 0
tBCellN = 0

do i = 1,outside_particle_num
 
 if ( flag(outside_particle(i)) == 0 ) then
  CellN = dint( floor(Xsol(outside_particle(i)))+Box(1)*floor(Ysol(outside_particle(i)))&
          +Box(1)*Box(2)*floor(Zsol(outside_particle(i)))+1 )  ! Identify cell
  tACellN(CellN) = tACellN(CellN) + 1                                                  ! Total number of A in cell N
  sACellN(CellN,tACellN(CellN)) = outside_particle(i)                                                    ! Specify/label newly added A with i

 else if ( flag(outside_particle(i)) == 1 ) then
   CellN = dint( floor(Xsol(outside_particle(i)))+Box(1)*floor(Ysol(outside_particle(i)))&
           +Box(1)*Box(2)*floor(Zsol(outside_particle(i)))+1 )  ! Identify cell
   tBCellN(CellN) = tBCellN(CellN) + 1                                                  ! Total number of B in cell N
   sBCellN(CellN,tBCellN(CellN)) = outside_particle(i)                                                    ! Specify/label newly added B with i
 end if

end do

do i = 1,inside_particle_num

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
  tACellN(CellN) = tACellN(CellN) + 1 ! Total number of A in cell N
  sACellN(CellN,tACellN(CellN)) = inside_particle(i) ! Specify/label newly added A with i

 else if ( flag(inside_particle(i)) == 1 ) then
   CellN = dint( floor(Xsol(inside_particle(i)))+Box(1)*floor(Ysol(inside_particle(i)))&
           +Box(1)*Box(2)*floor(Zsol(inside_particle(i)))+1 )  ! Identify cell
   tBCellN(CellN) = tBCellN(CellN) + 1 ! Total number of B in cell N
   sBCellN(CellN,tBCellN(CellN)) = inside_particle(i) ! Specify/label newly added B with i
 end if
 end if

end do

end subroutine distribute

end module distribute_mod
