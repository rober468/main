module global_mod
implicit none

!This module contains all the necessary global variables that are used within the main program and appropriate subroutines. Some of the subroutines have their own variables that are not used within the main program.

! Precision
integer, parameter :: sp = selected_real_kind(6,37)             ! Single precision
integer, parameter :: dp = selected_real_kind(15,307)           ! Double precision

! Pseudo Random number generator
integer(kind=selected_int_kind(9)) :: seed                      ! Random seed

! Variables not specific to simulation details
integer(kind=dp) :: time                                        ! Time loop variable
integer , dimension(8) :: values                                ! Physical time in simulation
real(kind=dp) , parameter :: pi = 4.d0*datan(1.d0)              ! pi
integer :: m                                                    ! Loop variable
logical :: timestep_file

! Number of molecules
real(kind=dp) :: rhoa                                           ! Initial density of particle A
real(kind=dp) :: rhob                                           ! Initial density of particle B
integer(kind=dp) :: Na0                                         ! Initial number of particle A
integer(kind=dp) :: Nb0                                         ! Initial number of particle B
integer(kind=sp) :: NT                                          ! Total number of particles
integer(kind=dp) :: iTotA                                       ! Total number of A
integer(kind=dp) :: iTotB                                       ! Total number of B

! Simulation box
real(kind=dp) , dimension(1:3) :: Box                           ! Length of simulation box
integer(kind=dp) :: Vol                                         ! Volume of simulation box
real(kind=dp) :: Vact                                           ! Volume of sim. box w/o monomer vol
integer(kind=dp) :: mpcnumcell                                     ! Number of subcells
integer(kind=dp) :: rxnnumcell                                     ! Number of subcells
real(kind=dp) :: x_width                                        ! Width of x slice, concentration fields
real(kind=dp) :: y_width                                        ! Width of y slice, concentration fields

! Dimer and Solvent Characteristics
integer :: numdims                                              ! Number of dimers
real(kind=dp) :: Ma                                             ! Mass of one A particle
real(kind=dp) :: Mb                                             ! Mass of one B particle
real(kind=dp) :: Mc                                             ! Mass of C sphere
real(kind=dp) :: Mn                                             ! Mass of N sphere
real(kind=dp) :: Rc                                             ! Radius of C sphere
real(kind=dp) :: Rn                                             ! Radius of N sphere
real(kind=dp) :: Dcn                                            ! Internuclear distance
real(kind=dp) :: EaC                                            ! Reduced LJ Potential of A with C
real(kind=dp) :: EbC                                            ! Reduced LJ Potential of B with C
real(kind=dp) :: EaN                                            ! Reduced LJ Potential of A with N
real(kind=dp) :: EbN                                            ! Reduced LJ Potential of B with N
real(kind=dp) :: zdimeq                                         ! Equil. z distance from wall
real(kind=dp) :: k_harm                                         ! Wall potential force constant
real(kind=dp) :: gap                                            ! Gap between motor and wall
real(kind=dp) :: T                                              ! Reduced temperature
real(kind=dp) :: varA                                           ! Variance of F vel. Gauss. distr.
real(kind=dp) :: varB                                           ! Variance of I vel. Gauss. distr.
real(kind=dp) :: mean                                           ! Mean of vel. Gaussian distr.
integer , parameter , dimension(2) :: sphere = (/1,2/)          ! Sphere identity; C=1;N=2

! MPC and MD time steps, cell size, and rotation angle
real(kind=dp) :: a0                                             ! MPC cell size
real(kind=dp) :: rotang                                         ! MPC rotation angle
real(kind=dp) :: MDt                                            ! MD time step
real(kind=dp) :: MPCt                                           ! MPC collision time
integer(kind=sp) :: MPCtstep                                    ! MPC step happens every nth step
integer(kind=dp) :: tstep                                       ! Total number of time steps
integer(kind=sp) :: cstep                                       ! Current time step
integer :: sort_step

! Dimer and solvent positions, velocities, forces, identity
integer(kind=sp) , dimension(:) , allocatable :: flag           ! Identity of solvent
real(kind=dp) , dimension(:,:) , allocatable :: Xdim            ! Monomer position, x
real(kind=dp) , dimension(:,:) , allocatable :: Ydim            ! Monomer position, y
real(kind=dp) , dimension(:,:) , allocatable :: Zdim            ! Monomer position, z
real(kind=dp) , dimension(:,:) , allocatable :: VXdim           ! Monomer velocity, x
real(kind=dp) , dimension(:,:) , allocatable :: VYdim           ! Monomer velocity, y
real(kind=dp) , dimension(:,:) , allocatable :: VZdim           ! Monomer velocity, z
real(kind=dp) , dimension(:,:) , allocatable :: FXdim           ! Force on dimer, x
real(kind=dp) , dimension(:,:) , allocatable :: FYdim           ! Force on dimer, y
real(kind=dp) , dimension(:,:) , allocatable :: FZdim           ! Force on dimer, z
real(kind=dp) , dimension(:) , allocatable :: Xsol              ! Solvent position, x
real(kind=dp) , dimension(:) , allocatable :: Ysol              ! Solvent position, y
real(kind=dp) , dimension(:) , allocatable :: Zsol              ! Solvent position, z
real(kind=dp) , dimension(:) , allocatable :: VXsol             ! Solvent velocity, x
real(kind=dp) , dimension(:) , allocatable :: VYsol             ! Solvent velocity, y
real(kind=dp) , dimension(:) , allocatable :: VZsol             ! Solvent velocity, z
real(kind=dp) , dimension(:,:) , allocatable :: motorCM
real(kind=dp) :: sep

! Distribution of molecules to cells
integer(kind=dp) , parameter :: MaxA=100                        ! Predicted maximum of F and G in each cell
integer(kind=dp) , parameter :: MaxB=100
integer(kind=dp) , dimension(:) , allocatable :: tACellN        ! Total number of F solvent in cell N
integer(kind=dp) , dimension(:) , allocatable :: tBCellN        ! Total number of I solvent in cell N
integer(kind=dp) , dimension(:,:) , allocatable :: sACellN      ! Label F in Cell N with i
integer(kind=dp) , dimension(:,:) , allocatable :: sBCellN      ! Label I in Cell N with i

! Neighbour list
real(kind=dp) :: skinrad                                        ! Skin radius for neighbour list
real(kind=dp) , dimension(:) , allocatable :: max_solvent_dist  ! Max distance of solvent
integer :: time_penetrate                                       ! Time for outside particle to penetrate inside skin
integer :: time_count                                           ! Time counter to determine if less then time_penetrate
integer(kind=dp) , dimension(:,:), allocatable :: neighbour_num ! Number of solvent particles in list
integer(kind=dp) , dimension(:,:,:) , allocatable :: neighbour  ! Identify solvent number in list
integer(kind=dp) , parameter :: max_neighbour = 30000           ! Maximum number of neighbours
integer(kind=dp) :: outside_particle_num
integer(kind=dp) , dimension(:) , allocatable :: outside_particle
integer(kind=dp) :: inside_particle_num
integer(kind=dp) , dimension(:) , allocatable :: inside_particle

! Monomer sphere quantities
real(kind=dp) , parameter :: cutoff = 2.d0**(1.d0/6.d0)         ! Cutoff prefactor for LJ pot.
real(kind=dp) :: cutoffC                                        ! Cutoff for C sphere LJ pot.
real(kind=dp) :: cutoffN                                        ! Cutoff for N sphere LJ pot.
real(kind=dp) :: cutoffC2                                       ! Cut. C sphere LJ pot. squared
real(kind=dp) :: cutoffN2                                       ! Cut. N sphere LJ pot. squared
real(kind=dp) :: TPotential                                     ! Total potential

! Bulk and C sphere reactions
real(kind=dp) :: k1                                             ! Back reaction B -> A
real(kind=dp) :: RPab                                           ! Probability of reaction from A to B
real(kind=dp) :: RPba                                           ! Probability of reaction from B to A

! Freqeuncy of file writing
integer(kind=sp) :: Freq1                                       ! Files for dim. written at nth time step
integer(kind=sp) :: Freq2                                       ! Files for solv. F/G/I written at nth time step
integer(kind=sp) :: Freq3                                       ! Files for E and Mom. written at nth time step
integer(kind=sp) :: Freq4                                       ! Files for vel.on Dcn written at nth time step

! Arrays that are too large to be allocated over and over again (for performance) and others
real(kind=dp) , dimension(:,:,:) , allocatable :: Potential     ! Potential per particle
real(kind=dp) , dimension(:,:,:) , allocatable :: FpairX        ! LJ force in x direction
real(kind=dp) , dimension(:,:,:) , allocatable :: FpairY        ! LJ force in y direction
real(kind=dp) , dimension(:,:,:) , allocatable :: FpairZ        ! LJ force in z direction
real(kind=dp) , dimension(:) , allocatable :: random            ! Vector of random numbers
real(kind=dp) , dimension(:) , allocatable :: VXCM              ! Center of mass velocity, x
real(kind=dp) , dimension(:) , allocatable :: VYCM              ! Center of mass velocity, y
real(kind=dp) , dimension(:) , allocatable :: VZCM              ! Center of mass velocity, z
real(kind=dp) , dimension(:) , allocatable :: tCellShift        ! Tot.solv.in cellN after shift
integer(kind=dp) , dimension(:) , allocatable :: NewCellN       ! ith particle's new cell N
real(kind=dp) , dimension(:) , allocatable :: UVX               ! Unit vector, x
real(kind=dp) , dimension(:) , allocatable :: UVY               ! Unit vector, y
real(kind=dp) , dimension(:) , allocatable :: UVZ               ! Unit vector, z
real(kind=dp) , dimension(:) , allocatable :: PXsol             ! Momentum of solvent, x
real(kind=dp) , dimension(:) , allocatable :: PYsol             ! Momentum of solvent, y
real(kind=dp) , dimension(:) , allocatable :: PZsol             ! Momentum of solvent, z
real(kind=dp) , dimension(:) , allocatable :: concA1d
real(kind=dp) , dimension(:) , allocatable :: concB1d
real(kind=dp) , dimension(:,:) , allocatable :: concA2d
real(kind=dp) , dimension(:,:) , allocatable :: concB2d
end module global_mod
