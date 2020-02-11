module sort_hilbert_mod
implicit none

contains

subroutine sort_hilbert(vol,NT,box,xsol,ysol,zsol,flag,type_B,vxsol,vysol,vzsol)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for in/out variables
integer(kind=dp) , intent(in) :: vol
integer(kind=sp) , intent(in) :: NT
real(kind=dp) , dimension(:) , intent(in) :: box
real(kind=dp) , dimension(:) , intent(inout) :: xsol
real(kind=dp) , dimension(:) , intent(inout) :: ysol
real(kind=dp) , dimension(:) , intent(inout) :: zsol
integer , dimension(:) , intent(inout) :: flag
integer , dimension(:) , intent(inout) :: type_B
real(kind=dp) , dimension(:) , intent(inout) :: vxsol
real(kind=dp) , dimension(:) , intent(inout) :: vysol
real(kind=dp) , dimension(:) , intent(inout) :: vzsol

! Definitions for internal variables
integer :: i,j
integer :: new
integer , parameter :: cell_max = 100
integer :: celln
integer(kind=dp) , dimension(:) , allocatable :: tcellN
integer(kind=dp) , dimension(:,:) , allocatable :: scellN
integer(kind=dp) , dimension(:) , allocatable :: newvol
integer(kind=sp) , dimension(:) , allocatable :: newbox
real(kind=dp) , dimension(:) , allocatable :: newxsol
real(kind=dp) , dimension(:) , allocatable :: newysol
real(kind=dp) , dimension(:) , allocatable :: newzsol
integer , dimension(:) , allocatable :: newflag
integer , dimension(:) , allocatable :: newtype_B
real(kind=dp) , dimension(:) , allocatable :: newvxsol
real(kind=dp) , dimension(:) , allocatable :: newvysol
real(kind=dp) , dimension(:) , allocatable :: newvzsol

allocate(tcelln(vol),scelln(vol,cell_max),newvol(vol),newbox(3),newxsol(NT),&
         newysol(NT),newzsol(NT),newflag(NT),newtype_B(NT),newvxsol(NT),&
         newvysol(NT),newvzsol(NT))

! Distribute particles to cells
tcellN = 0
scellN = 0

do i = 1,NT 
 celln = int(floor(xsol(i)) + box(1)*floor(ysol(i)) + box(1)*box(2)*floor(zsol(i)) + 1) !Identify cell
 tcellN(celln) = tcelln(celln) + 1
 scelln(celln,tcelln(celln)) = i
end do

! Find the order of the Hilbert curve through the cells
open(unit=10,file="hilbert.txt")
 do i = 1,vol
  read(10,*)newbox(1),newbox(2),newbox(3)
  newvol(i) = int(newbox(1) + box(1)*newbox(2) + box(1)*box(2)*newbox(3) + 1)
 end do
close(10)

! Reorganize into new order
new = 0
do i = 1,vol
 do j=1,tcelln(newvol(i))
  new = new + 1
  newxsol(new) = xsol(scelln(newvol(i),j))
  newysol(new) = ysol(scelln(newvol(i),j))
  newzsol(new) = zsol(scelln(newvol(i),j))
  newflag(new) = flag(scelln(newvol(i),j))
  newtype_B(new) = type_B(scelln(newvol(i),j))
  newvxsol(new) = vxsol(scelln(newvol(i),j))
  newvysol(new) = vysol(scelln(newvol(i),j))
  newvzsol(new) = vzsol(scelln(newvol(i),j))
 end do
end do

! Equate back to original arrays
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(NT,xsol,ysol,zsol,flag,type_B,vxsol,vysol,vzsol) &
!$OMP SHARED(newxsol,newysol,newzsol,newflag,newtype_B,newvxsol,newvysol,newvzsol) &
!$OMP PRIVATE(i)
do i = 1,NT
 xsol(i) = newxsol(i)
 ysol(i) = newysol(i)
 zsol(i) = newzsol(i)
 flag(i) = newflag(i)
 type_B(i) = newtype_B(i)
 vxsol(i) = newvxsol(i)
 vysol(i) = newvysol(i)
 vzsol(i) = newvzsol(i)
end do
!OMP END PARALLEL DO

deallocate(tcelln,scelln,newvol,newbox,newxsol,newysol,newzsol,newflag,newtype_B,newvxsol,&
         newvysol,newvzsol)

end subroutine sort_hilbert

end module
