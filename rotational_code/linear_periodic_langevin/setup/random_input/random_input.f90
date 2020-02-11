program random_input
use mt95
implicit none

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Other variables
integer , parameter :: k = 3
integer , parameter :: sep = 4
integer , parameter :: eps = 2
integer , parameter :: runs = 50
integer :: tot
integer , dimension(:) , allocatable :: random
integer(kind=selected_int_kind(9)) :: seed
integer :: i

! derived quantities
tot = k * sep * eps * runs

allocate(random(tot))


! initialize random number generator with seed
seed = 60456845
call genrand_init( put=seed )

do i = 1,tot
 call genrand_int32(random(i))
end do

open(unit=10,file="random.txt",status="replace",position="append")
do i = 1,tot
 write(10,*)abs(random(i))
end do
close(10)

deallocate(random)

end program random_input 
