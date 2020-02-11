module global_mod
implicit none

! precision
integer, parameter :: sp = selected_real_kind(6,37)       ! single precision
integer, parameter :: dp = selected_real_kind(15,307)     ! double precision

! parameters
real(kind=dp), parameter :: pi = 4.*atan(1.d0)

! integer variables
integer :: i                                            ! loop variable
integer :: numrotors                                    ! total number of rotors
integer :: numsteps                                     ! number of timesteps in each simulation window
integer :: cstep                                        ! total number of timesteps
integer :: seed                                         ! random number generator seed
integer :: freq1                                        ! frequency of angle recording
integer :: freq2                                        ! frequency of restart data recording
logical :: data_file                                    ! denotes if data file exists for restart

! double precision real variables
real(kind=dp) :: timestep                               ! timestep
real(kind=dp) :: kBT                                    ! temperature
real(kind=dp) :: Rc                                     ! radius of catalytic sphere
real(kind=dp) :: sigma_Rc                               ! radius of C sphere with interaction zone
real(kind=dp) :: Rn                                     ! radius of noncatalytic sphere
real(kind=dp) :: sigma_Rn                               ! radius of N sphere with interaction zone
real(kind=dp) :: d_CN                                   ! distance between C and N sphere for each rotor
real(kind=dp) :: r_CN                                   ! r_CN = 0.5 * d_CN
real(kind=dp) :: len_rotor                              ! length of each rotor, boundary layer included
real(kind=dp) :: sep                                    ! separation between rotors
real(kind=dp) :: sim_box_len                            ! simulation box length
real(kind=dp) :: half_sim_box_len                       ! half of the simulation box length
real(kind=dp) :: dens                                   ! total number density of solvent
real(kind=dp) :: mean                                   ! mean of standard normal distribution
real(kind=dp) :: var                                    ! variance of standard normal distribution
real(kind=dp) :: Ean                                    ! LJ potential of A with N sphere
real(kind=dp) :: Ebn                                    ! LJ potential of B with N sphere
real(kind=dp) :: a                                      ! MPC cell size
real(kind=dp) :: alpha                                  ! MPC rotation angle
real(kind=dp) :: D                                      ! diffusion coefficient of solvent
real(kind=dp) :: D_R                                    ! rotation diffusion coefficient of rotors
real(kind=dp) :: k1                                     ! solvent forward reaction rate A -> B
real(kind=dp) :: k_1                                    ! solvent reverse reaction rate B -> A
real(kind=dp) :: k_D                                    ! Smoluchowski diffusion rate constant
real(kind=dp) :: c_A_SS                                 ! steady-state concentration of A at infinity
real(kind=dp) :: c_B_SS                                 ! steady-state concentration of B at infinity
real(kind=dp) :: kappa                                  ! reaction exp. decay constant
real(kind=dp) :: m                                      ! mass of solvent particles
real(kind=dp) :: tau                                    ! MPC timestep
real(kind=dp) :: derj_len_N                             ! Derjaguin length for N sphere
real(kind=dp) :: p_plus                                 ! probability of forward reaction at motor
real(kind=dp) :: p_minus                                ! probability of backward reaction at motor
real(kind=dp) :: Z
real(kind=dp), dimension(:), allocatable :: angle       ! angle of each rotor
real(kind=dp), dimension(:), allocatable :: F           ! "force" on N sphere; contains sin_theta for torque
real(kind=dp), dimension(:), allocatable :: w           ! standard normal random variable
real(kind=dp), dimension(:,:), allocatable :: rotor_CM  ! rotor centre-of-mass coordinates

end module global_mod
