module reactive_grid_mod
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reactive_grid_setup(box,react_box,react_box_type)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for external variables
real(kind=dp) , dimension(:) , intent(in) :: box
integer , dimension(:) , intent(out) :: react_box
integer :: react_box_type

! Definitions for internal variables 
integer :: i
integer :: react_box_num

! Derived quantities
react_box_num = dint( box(1)*box(2) )


react_box = 0

if ( react_box_type == 1 ) then 

 ! Reactive grid 1 - Left half in x-direction

 do i = 1,react_box_num
  if ( mod(i,int(box(1))) > 0 .and. mod(i,int(box(1))) <= int(0.5*box(1)) ) then
   react_box(i) = 1
  end if
 end do

else if ( react_box_type == 2 ) then

 ! Reactive grid 2 - checkerboard 

 do i = 1,react_box_num
  if ( mod(i+floor(i/box(2)),2) == 0 ) then
   react_box(i) = 1
  end if
 end do

else if ( react_box_type == 3 ) then

 ! Reactive grid 3 - big circle

 do i = 1,react_box_num
  if ( (0.5*box(1)-mod(i-1,int(box(1)))+0.5)**2+(0.5*box(2)-floor(i/box(1))+0.5)**2 < 20.**2 ) then
   react_box(i) = 1
  end if
 end do

else if ( react_box_type == 4 ) then

 ! Reactive grid 4 - small circles

 do i = 1,react_box_num
  if ( ((1./6.)*box(1)-mod(i-1,int(box(1)))+0.5)**2+(0.25*box(2)-floor(i/box(1))+0.5)**2 < 10.**2 .and. &
       ((1./6.)*box(1)-mod(i-1,int(box(1)))+0.5)**2+(0.75*box(2)-floor(i/box(1))+0.5)**2 < 10.**2 .and. &
       (0.5*box(1)-mod(i-1,int(box(1)))+0.5)**2+(0.25*box(2)-floor(i/box(1))+0.5)**2 < 10.**2 .and. &
       (0.5*box(1)-mod(i-1,int(box(1)))+0.5)**2+(0.75*box(2)-floor(i/box(1))+0.5)**2 < 10.**2 .and. &
       ((5./6.)*box(1)-mod(i-1,int(box(1)))+0.5)**2+(0.25*box(2)-floor(i/box(1))+0.5)**2 < 10.**2 .and. &
       ((5./6.)*box(1)-mod(i-1,int(box(1)))+0.5)**2+(0.75*box(2)-floor(i/box(1))+0.5)**2 < 10.**2 ) then
   react_box(i) = 1
  end if
 end do

else if ( react_box_type == 5 ) then

 ! Reactive grid 5 - strip

 do i = 1,react_box_num
  if ( mod(i,int(box(1))) > floor(0.5*box(1)-4) .and. mod(i,int(box(1))) <= floor(0.5*box(1)+4) ) then
   react_box(i) = 1
  end if
 end do

else if ( react_box_type == 6 ) then

 ! Reactive grid 6 - whole wall

 do i = 1,react_box_num
  react_box(i) = 1
 end do

else if ( react_box_type == 7 ) then

 ! Reactive grid 7 - no reactive wall (do nothing

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reactive_grid(react_box,react_dir,box,xsol,ysol,vxsol,vysol,flag,time_aft_coll)
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Definitions for external variables
integer , dimension(:) , intent(in) :: react_box
integer , intent(in) :: react_dir
real(kind=dp) , dimension(:) , intent(in) :: box
real(kind=dp) , intent(in) :: xsol
real(kind=dp) , intent(in) :: ysol
real(kind=dp) , intent(in) :: vxsol
real(kind=dp) , intent(in) :: vysol
integer , intent(inout) :: flag
real(kind=dp) , intent(in) :: time_aft_coll

! Definitions for internal variables
real(kind=dp) :: adj_xsol
real(kind=dp) :: adj_ysol
integer :: hit_box


! Adjust position of particle back to wall ; incl. PBCs
adj_xsol = xsol - vxsol*time_aft_coll
adj_ysol = ysol - vysol*time_aft_coll
if ( adj_xsol >= box(1) ) then
 adj_xsol = adj_xsol - box(1)
else if ( adj_xsol < 0. ) then
 adj_xsol = adj_xsol + box(1)
end if
if ( adj_ysol >= box(2) ) then
 adj_ysol = adj_ysol - box(2)
else if ( adj_ysol < 0. ) then
 adj_ysol = adj_ysol + box(2)
end if

! Determine which reactive box patch particle hits and see if reaction occurs
hit_box = int( floor(adj_xsol) + box(1)*floor(adj_ysol) + 1 )
if ( react_box(hit_box) == 1 ) then
 if ( react_dir == 0 .and. flag == 0 ) then
  flag = 1
 else if ( react_dir == 1 .and. flag == 1 ) then
  flag = 0
 end if
end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module
