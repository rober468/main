program main_pattern
use mtmod
use global_mod
use neighbour_list_mod
use reactive_grid_mod
use initial_mod
use restart_mod
use evolve_mod
use distribute_mod
use mpc_rotation_mod
use reaction_mod
use output_mod
use sort_hilbert_mod
implicit none

integer :: omp_get_num_threads

! Read initial conditions
open(unit=10,file='input.txt',status='old',action='read')
read(10,*)ini_seed
read(10,*)tstep
read(10,*)rhoa
read(10,*)rhob
read(10,*)rhos
read(10,*)Box(1)
read(10,*)Box(2)
read(10,*)Box(3)
read(10,*)numdims
read(10,*)Ma
read(10,*)Mb
read(10,*)Ms
read(10,*)Rc
read(10,*)Rn
read(10,*)Dcn
read(10,*)EaC
read(10,*)EbC
read(10,*)EsC
read(10,*)EaN
read(10,*)EbN
read(10,*)EsN
read(10,*)Rw
read(10,*)Ew
read(10,*)Rm_delta
read(10,*)Em
read(10,*)RPab
read(10,*)RPba
read(10,*)T
read(10,*)skinrad
read(10,*)a0
read(10,*)rotang
read(10,*)MDt
read(10,*)MPCt
read(10,*)sort_step
read(10,*)k1
read(10,*)x_conc_slab
read(10,*)y_conc_slab
read(10,*)z_conc_slab
read(10,*)react_box_type
read(10,*)react_dir
read(10,*)Freq1
read(10,*)Freq2
read(10,*)Freq3
read(10,*)Freq4
close(10)

Vol = dint( Box(1)*Box(2)*Box(3) )
Vact = Vol - numdims*((4.*pi*Rc**3)/3. + (4.*pi*Rn**3)/3.)
Na0 = floor(rhoa * Vact)
Nb0 = floor(rhob * Vact)
Ns0 = floor(rhob * Vact)
NT=Na0+Nb0+Ns0
mpcnumcell = dint( (Box(1)/a0)*(Box(2)/a0)*(Box(3)/a0) + (Box(1)/a0)*(Box(2)/a0) )
rxnnumcell = dint( (Box(1)/a0)*(Box(2)/a0)*(Box(3)/a0) )
x_width = 0.5*a0
y_width = 0.5*a0
z_width = 0.5*a0

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(num_threads)
 num_threads = omp_get_num_threads()
!$OMP END PARALLEL

!Allocate arrays
allocate(flag(1:NT))
allocate(Xdim(1:numdims,1:2),Ydim(1:numdims,1:2),Zdim(1:numdims,1:2))
allocate(VXdim(1:numdims,1:2),VYdim(1:numdims,1:2),VZdim(1:numdims,1:2))
allocate(FXdim(1:numdims,1:2),FYdim(1:numdims,1:2),FZdim(1:numdims,1:2))
allocate(Xsol(1:NT),Ysol(1:NT),Zsol(1:NT))
allocate(VXsol(1:NT),VYsol(1:NT),VZsol(1:NT))
allocate(tACellN(1:Vol),sACellN(1:Vol,1:MaxA))
allocate(tBCellN(1:Vol),sBCellN(1:Vol,1:MaxB))
allocate(tSCellN(1:Vol),sSCellN(1:Vol,1:MaxS))
allocate(Potential(1:numdims,1:2,1:max_neighbour))
allocate(FpairX(1:numdims,1:2,1:max_neighbour),FpairY(1:numdims,1:2,1:max_neighbour),&
FpairZ(1:numdims,1:2,1:max_neighbour))
allocate(VXCM(1:mpcnumcell),VYCM(1:mpcnumcell),VZCM(1:mpcnumcell))
allocate(tCellShift(1:mpcnumcell))
allocate(NewCellN(1:NT))
allocate(UVX(1:mpcnumcell),UVY(1:mpcnumcell),UVZ(1:mpcnumcell))
allocate(PXsol(1:NT),PYsol(1:NT),PZsol(1:NT))
allocate(neighbour(1:numdims,1:2,1:max_neighbour),max_solvent_dist(1:NT),neighbour_num(1:numdims,1:2))
allocate(concA1d(1:int(box(1)/x_width)),concB1d(1:int(box(1)/x_width)))
allocate(concA2d_xy(1:int(box(1)/x_width),1:int(box(2)/y_width)),&
concB2d_xy(1:int(box(1)/x_width),1:int(box(2)/y_width)))
allocate(concA2d_xz(1:int(box(1)/x_width),1:int(box(3)/z_width)),&
concB2d_xz(1:int(box(1)/x_width),1:int(box(3)/z_width)))
allocate(concA2d_yz(1:int(box(2)/y_width),1:int(box(3)/z_width)),&
concB2d_yz(1:int(box(2)/y_width),1:int(box(3)/z_width)))
allocate(outside_particle(1:NT))
allocate(inside_particle(1:(numdims*2*max_neighbour)))
allocate(motorCM(numdims,3))
allocate(react_box(1:int(box(1))*int(box(2))))
allocate(mt(1:N,num_threads),mti(num_threads),seed(num_threads))


! Quantities derived from input parameters
MPCtstep = dint(MPCt/MDt)                               !Number of time step for MPC calc.
Mc = (NT/dble(Vol)) * ((4.*pi)/3)*Rc**3                 !Mass of C sphere
Mn = (NT/dble(Vol)) * ((4.*pi)/3)*Rn**3                 !Mass of N sphere
cutoffC = Rc * cutoff                                   !Cutoff for C sphere LJ pot.
cutoffN = Rn * cutoff                                   !Cutoff for N sphere LJ pot.
cutoffC2 = cutoffC**2                                   !Cut. C sphere LJ pot. squared
cutoffN2 = cutoffN**2                                   !Cut. N sphere LJ pot. squared
varA = T / Ma                                           !Variance of A velocity Gaussian distribution
varB = T / Mb                                           !Variance of B velocity Gauss. distr.
varS = T / Ms                                           !Variance of S velocity Gauss. distr.
mean = 0.                                               !Initial mean

! Create or read file with current timestep and numbers of particles
inquire(file='timestep.dat',exist=timestep_file)
if ( timestep_file ) then
 open(unit=10,file='timestep.dat',status='old',action='read')
  read(10,*)cstep
  read(10,*)iTotA
  read(10,*)iTotB
  read(10,*)iTotS
 close(10)
call restart(cstep,iTotA,iTotB,iTotS,NT,MDt,Box,vol,flag,sphere,numdims,cutoff,Rc,Rn,EaC,EbC,EsC,&
EaN,EbN,EsN,Rw,Ew,Rm_delta,Em,RPab,RPba,Xdim,Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,Xsol,&
Ysol,Zsol,VXsol,VYsol,VZsol,TPotential,Potential,FpairX,FpairY,FpairZ,motorCM,neighbour_num,&
neighbour,max_neighbour,skinrad,max_solvent_dist,time_count,time_penetrate,outside_particle_num,&
outside_particle,inside_particle_num,inside_particle,sort_step,Freq1,Freq4,data_id,sub_group_id,&
sys_dspace_id,sys_crp_list,sys_dset_id,conc_dspace_id,conc_crp_list,conc_dset_id,solv_dspace_id,&
solv_crp_list,solv_dset_id,dimer_dspace_id,dimer_crp_list,dimer_dset_id)
call reactive_grid_setup(box,react_box,react_box_type)
call mtgetf(mt,mti,"mt_state.dat","u")
else
 ! Random number generator seeding
 mti = N1
 call random_seed(size=seed_size)
 allocate(ini_seed_arr(seed_size))
 do i = 1,seed_size
  ini_seed_arr(i) = ini_seed + i
 end do
 call random_seed(put=ini_seed_arr)
 do i = 1,num_threads
     call random_number(seed_tr)
     seed(i) = floor(nn*(seed_tr+1))
     call sgrnd(seed(i),mt(:,i),mti(i))
 end do
 deallocate(ini_seed_arr)
 ! Initialize simulation
 cstep = 0
 call output_hdf5_initiate(cstep,Freq1,Freq2,Freq3,Freq4,box,T,MDt,MPCt,NT,numdims,x_width,y_width,z_width)
 call initial(Box,numdims,Dcn,Na0,Nb0,Ns0,NT,cutoff,Rc,Rn,Rw,varA,varB,varS,mean,flag,Xdim,Ydim,Zdim, &
                   VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,Xsol,Ysol,Zsol,VXsol,VYsol,VZsol, &
                   FpairX,FpairY,FpairZ,motorCM,mt(:,1),mti(1))
 call reactive_grid_setup(box,react_box,react_box_type)
 ! Sorting algorithm
 if ( mod(cstep,sort_step) == 0 ) then
  call sort_hilbert(vol,NT,box,xsol,ysol,zsol,flag,vxsol,vysol,vzsol)
 end if
 neighbour_num = 0
 outside_particle_num = 0
 inside_particle_num = 0
 do k = 1,numdims
  ! Catalytic sphere
  call neighbour_list(NT,Box,cutoffC,skinrad,Xsol,Ysol,Zsol,Xdim(k,1),Ydim(k,1),Zdim(k,1),neighbour_num(k,1),&
neighbour(k,1,:))
  ! Noncatalytic sphere
  call neighbour_list(NT,Box,cutoffN,skinrad,Xsol,Ysol,Zsol,Xdim(k,2),Ydim(k,2),Zdim(k,2),neighbour_num(k,2),&
neighbour(k,2,:))
 end do
 call inout_particle_list(NT,box,cutoffC,cutoffN,skinrad,xsol,ysol,zsol,xdim,ydim,zdim,&
outside_particle_num,outside_particle,inside_particle_num,inside_particle,numdims)
 max_solvent_dist = 0
 time_penetrate = floor((0.5*skinrad / (sqrt(maxval(VXsol*VXsol+VYsol*VYsol+VZsol*VZsol)))) / MDt)
 time_count = 0
 call output_hdf5_write(cstep,Vol,Freq1,Freq2,Freq3,Freq4,NT,numdims,Ma,Mb,Ms,Mc,Mn,x_conc_slab,&
                        y_conc_slab,z_conc_slab,box,x_width,y_width,z_width,xsol,ysol,zsol,flag,&
                        vxsol,vysol,vzsol,cutoffC,Dcn,concA1d,concB1d,concA2d_xy,concB2d_xy,concA2d_xz,&
                        concB2d_xz,concA2d_yz,concB2d_yz,pxsol,pysol,pzsol,xdim,ydim,zdim,vxdim,vydim,vzdim,&
                        Tpotential)
 call date_and_time(VALUES=values)
 open(unit=10,file='Time.dat',status='new',action='write')
  write(10,*) values(5),values(6),values(7), cstep
 close(10)

end if

! Main time loop
do time = 1,tstep
 cstep = cstep + 1
 time_count = time_count + 1.
  if ( mod(cstep,MPCtstep) /= 0 ) then
   call evolve(cstep,time,NT,MDt,MPCtstep,numdims,Ma,Mb,Ms,Mc,Mn,Dcn,Box,cutoff,Rc,Rn,EaC,EbC,EsC,EaN,EbN,EsN,&
                  Rw,Ew,Rm_delta,Em,RPab,RPba,sphere,flag,xsol,Ysol,Zsol,VXsol,VYsol,VZsol,Xdim,&
                  Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,motorCM,TPotential,Potential,FpairX,&
                  FpairY,FpairZ,skinrad,max_solvent_dist,neighbour_num,neighbour,time_count,&
                  time_penetrate,outside_particle_num,outside_particle,inside_particle_num,inside_particle,&
                  vol,sort_step,react_box,react_dir)
  else if ( mod(cstep,MPCtstep) == 0 ) then
   call evolve(cstep,time,NT,MDt,MPCtstep,numdims,Ma,Mb,Ms,Mc,Mn,Dcn,Box,cutoff,Rc,Rn,EaC,EbC,EsC,EaN,EbN,EsN,&
                  Rw,Ew,Rm_delta,Em,RPab,RPba,sphere,flag,xsol,Ysol,Zsol,VXsol,VYsol,VZsol,Xdim,&
                  Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,motorCM,TPotential,Potential,FpairX,&
                  FpairY,FpairZ,skinrad,max_solvent_dist,neighbour_num,neighbour,time_count,&
                  time_penetrate,outside_particle_num,outside_particle,inside_particle_num,inside_particle,&
                  vol,sort_step,react_box,react_dir)
   call distribute(NT,numdims,cutoff,Rc,Rn,flag,Xsol,Ysol,Zsol,Xdim,Ydim,Zdim,Box,tACellN,tBCellN,tSCellN, &
                   sACellN,sBCellN,sSCellN,cutoffC2,cutoffN2,cstep,outside_particle_num,outside_particle,&
                   inside_particle_num,inside_particle)
   call mpc_rotation(NT,rotang,mpcnumcell,a0,Box,Xsol,Ysol,Zsol,VXsol,VYsol,VZsol,VXCM,VYCM,VZCM, &
                     tCellShift,NewCellN,UVX,UVY,UVZ,mt,mti)
   call reaction_simple(MPCt,rxnnumcell,k1,flag,tBCellN,sBCellN,mt,mti)
  end if
! Output to file
  call output_hdf5_write(cstep,Vol,Freq1,Freq2,Freq3,Freq4,NT,numdims,Ma,Mb,Ms,Mc,Mn,x_conc_slab,&
                         y_conc_slab,z_conc_slab,box,x_width,y_width,z_width,xsol,ysol,zsol,flag,&
                         vxsol,vysol,vzsol,cutoffC,Dcn,concA1d,concB1d,concA2d_xy,concB2d_xy,concA2d_xz,&
                         concB2d_xz,concA2d_yz,concB2d_yz,pxsol,pysol,pzsol,xdim,ydim,zdim,vxdim,vydim,vzdim,&
                         Tpotential)
! Timing the simulation
  if ( mod(int(cstep),1000) == 0 ) then
   call date_and_time(VALUES=values)
   open(unit=10,file='Time.dat',status='old',action='write',position='append')
    write(10,*) values(5),values(6),values(7), cstep
   close(10)
  end if

end do

call calc_solv_type_num(NT,Vol,flag,iTotA,iTotB,iTotS)
open(unit=10,file='timestep.dat',status='replace')
 write(10,*)cstep
 write(10,*)iTotA
 write(10,*)iTotB
 write(10,*)iTotS
close(10)

call output_hdf5_close(numdims)
call mtsavef(mt,mti,"mt_state.dat","u")

!Deallocate arrays
deallocate(flag)
deallocate(Xdim,Ydim,Zdim)
deallocate(VXdim,VYdim,VZdim)
deallocate(FXdim,FYdim,FZdim)
deallocate(Xsol,Ysol,Zsol)
deallocate(VXsol,VYsol,VZsol)
deallocate(tACellN,sACellN)
deallocate(tBCellN,sBCellN)
deallocate(tSCellN,sSCellN)
deallocate(Potential,FpairX,FpairY,FpairZ)
deallocate(VXCM,VYCM,VZCM)
deallocate(tCellShift)
deallocate(NewCellN)
deallocate(UVX,UVY,UVZ)
deallocate(PXsol,PYsol,PZsol)
deallocate(neighbour,neighbour_num,max_solvent_dist)
deallocate(concA1d,concB1d)
deallocate(concA2d_xy,concB2d_xy)
deallocate(concA2d_xz,concB2d_xz)
deallocate(concA2d_yz,concB2d_yz)
deallocate(outside_particle)
deallocate(inside_particle)
deallocate(motorCM)
deallocate(react_box)
deallocate(mt,mti,seed)

end program main_pattern
