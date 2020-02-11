program main_langevin_rotors
use global_mod
use derived_quantities_mod
use initial_mod
use evolve_mod
use output_mod
implicit none


! read input parameters
open(unit=10,file="input.txt",status="old",action="read")
read(10,*)seed
read(10,*)numrotors             ! total number of rotors
read(10,*)numsteps              ! total number of timesteps
read(10,*)timestep              ! timestep
read(10,*)kBT                   ! temperature
read(10,*)Rc                    ! radius of catalytic sphere
read(10,*)Rn                    ! radius of noncatalytic sphere
read(10,*)d_CN
read(10,*)sep
read(10,*)dens                  ! total number density of solvent
read(10,*)D_R                   ! rotational diffusion coefficient
read(10,*)mean
read(10,*)var
read(10,*)Ean
read(10,*)Ebn
read(10,*)a
read(10,*)alpha
read(10,*)m
read(10,*)tau
read(10,*)k1
read(10,*)k_1
read(10,*)c_A_SS
read(10,*)c_B_SS
read(10,*)p_plus
read(10,*)p_minus
read(10,*)freq1
read(10,*)freq2
close(10)


! allocate arrays
allocate(angle(numrotors))
allocate(F(numrotors))
allocate(w(numrotors))
allocate(rotor_CM(numrotors,2))

! calculate derived parameters
!!!!!!!!!!! include half length of simulation box, space between rotors, etc.
! half distance between C and N sphere of each rotor
r_CN = 0.5 * d_CN
! length of one rotor
len_rotor = 2.**(1./6.)*(Rc+Rn) + d_CN
! simulation box length
sim_box_len = numrotors*(sep+len_rotor)
half_sim_box_len = 0.5*sim_box_len
! diffusion coefficient
call diffusion_coefficient_MPC(kBT,tau,alpha,dens,a,m,D)
! Smoluchowski diffusion reaction rate
k_D = 4.*pi*D*Rc
! Derjaguin length
call derjaguin_length_N(Rn,Ean,Ebn,derj_len_N)
! Kappa
call kappa_decay(k1,k_1,D,kappa)
! prefactor for e^{-kappa(r-R_c)}/r
call conc_exp_prefactor(p_plus,p_minus,Rc,kBT,m,k_D,D,kappa,c_A_SS,c_B_SS,Z)

! initialize or restart system from beginning or current timestep, respectively
! Create or read file with current timestep and numbers of particles
inquire(file="data.h5",exist=data_file)
if ( data_file ) then
! call restart()
else
 call initialize(seed,numrotors,r_CN,len_rotor,half_sim_box_len,sim_box_len,Rn,Rc,kappa,kBT,&
                derj_len_N,Z,mean,var,sep,cstep,rotor_CM,angle,F,w)
 call output_hdf5_initial(rotor_CM,timestep,kBT,Rc,Rn,d_CN,sep,dens,Ean,Ebn,D,D_R,k1,k_1,kappa,&
                         derj_len_N,p_plus,p_minus,Z,angle,F,w,numrotors,freq1,freq2)
 call output_hdf5_write(cstep,numrotors,angle,F,w,freq1,freq2)
end if

! time evolution
do i = 1,numsteps,1

 ! evolve system
 call evolve(timestep,numrotors,r_CN,D_R,kBT,rotor_CM,half_sim_box_len,sim_box_len,Rn,Rc,kappa,&
            derj_len_N,Z,mean,var,F,w,angle)
 cstep = cstep + 1
 ! output data
 call output_hdf5_write(cstep,numrotors,angle,F,w,freq1,freq2)

end do

! close HDF5 file
call output_HDF5_close

! deallocate arrays
deallocate(angle)
deallocate(F)
deallocate(w)
deallocate(rotor_CM)

end program main_langevin_rotors
