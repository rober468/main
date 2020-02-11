program main_rotation
use mt95
use global_mod
use neighbour_list_mod
use initial_mod
use restart_mod
use evolve_mod
use distribute_mod
use mpc_rotation_mod
use reaction_simple_mod
use output_mod
use sort_hilbert_mod
implicit none

! Read initial conditions
open(unit=10,file='input.txt',status='old',action='read')
read(10,*)seed
read(10,*)tstep
read(10,*)rhoa
read(10,*)rhob
read(10,*)Box(1)
read(10,*)Box(2)
read(10,*)Box(3)
read(10,*)numdims
read(10,*)Ma
read(10,*)Mb
read(10,*)Rc
read(10,*)Rn
read(10,*)Dcn
read(10,*)EaC
read(10,*)EbC
read(10,*)EaN
read(10,*)EbN
read(10,*)k_harm
read(10,*)gap
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
read(10,*)Freq1
read(10,*)Freq2
read(10,*)Freq3
read(10,*)Freq4
close(10)

sep = ( box(1) - numdims*(2.*Rc*cutoff+Dcn) ) / numdims
Vol = dint( Box(1)*Box(2)*Box(3) )
Vact = Vol - numdims*((4.*pi*Rc**3)/3. - (4.*pi*Rn**3)/3.)
Na0 = floor(rhoa * Vact)
Nb0 = floor(rhob * Vact)
NT=Na0+Nb0
mpcnumcell = dint( (Box(1)/a0)*(Box(2)/a0)*(Box(3)/a0) + (Box(1)/a0)*(Box(2)/a0) )
rxnnumcell = dint( (Box(1)/a0)*(Box(2)/a0)*(Box(3)/a0) )
zdimeq = Rc + gap
x_width = 0.5*a0
y_width = 0.5*a0

!Allocate arrays
allocate(flag(1:NT))
allocate(Xdim(1:numdims,1:2),Ydim(1:numdims,1:2),Zdim(1:numdims,1:2))
allocate(VXdim(1:numdims,1:2),VYdim(1:numdims,1:2),VZdim(1:numdims,1:2))
allocate(FXdim(1:numdims,1:2),FYdim(1:numdims,1:2),FZdim(1:numdims,1:2))
allocate(Xsol(1:NT),Ysol(1:NT),Zsol(1:NT))
allocate(VXsol(1:NT),VYsol(1:NT),VZsol(1:NT))
allocate(tACellN(1:Vol),sACellN(1:Vol,1:MaxA))
allocate(tBCellN(1:Vol),sBCellN(1:Vol,1:MaxB))
allocate(Potential(1:numdims,1:2,1:max_neighbour))
allocate(FpairX(1:numdims,1:2,1:max_neighbour),FpairY(1:numdims,1:2,1:max_neighbour),&
FpairZ(1:numdims,1:2,1:max_neighbour))
allocate(VXCM(1:mpcnumcell),VYCM(1:mpcnumcell),VZCM(1:mpcnumcell))
allocate(tCellShift(1:mpcnumcell))
allocate(NewCellN(1:NT))
allocate(UVX(1:mpcnumcell),UVY(1:mpcnumcell),UVZ(1:mpcnumcell))
allocate(PXsol(1:NT),PYsol(1:NT),PZsol(1:NT))
allocate(random(1:max_neighbour))
allocate(neighbour(1:numdims,1:2,1:max_neighbour),max_solvent_dist(1:NT),neighbour_num(1:numdims,1:2))
allocate(concA1d(1:int(box(1)/x_width)),concB1d(1:int(box(1)/x_width)))
allocate(concA2d(1:int(box(1)/x_width),1:int(box(2)/y_width)),&
concB2d(1:int(box(1)/x_width),1:int(box(2)/y_width)))
allocate(outside_particle(1:NT))
allocate(inside_particle(1:(numdims*2*max_neighbour)))
allocate(motorCM(numdims,3))
allocate(type_B(NT))

! Quantities derived from input parameters
MPCtstep = dint(MPCt/MDt)                               !Number of time step for MPC calc.
Mc = (NT/dble(Vol)) * ((4.*pi)/3)*Rc**3                 !Mass of C sphere
Mn = (NT/dble(Vol)) * ((4.*pi)/3)*Rn**3                 !Mass of N sphere
cutoffC = Rc * cutoff                                   !Cutoff for C sphere LJ pot.
cutoffN = Rn * cutoff                                   !Cutoff for N sphere LJ pot.
cutoffC2 = cutoffC**2                                   !Cut. C sphere LJ pot. squared
cutoffN2 = cutoffN**2                                   !Cut. N sphere LJ pot. squared
varA = T / Ma                                           !Variance of F velocity Gaussian distribution
varB = T / Mb                                           !Variance of inert velocity Gauss. distr.
mean = 0.                                               !Initial mean

! Initialize random number generator with seed
call genrand_init( put=seed )

! Create or read file with current timestep and numbers of particles
inquire(file='timestep.dat',exist=timestep_file)
if ( timestep_file ) then
 open(unit=10,file='timestep.dat',status='old',action='read')
  read(10,*)cstep
  read(10,*)iTotA
  read(10,*)iTotB
 close(10)
 call restart(cstep,iTotA,iTotB,NT,MDt,Box,vol,flag,type_B,sphere,numdims,cutoff,Rc,Rn,EaC,EbC,&
EaN,EbN,zdimeq,k_harm,RPab,RPba,Xdim,Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,Xsol,&
Ysol,Zsol,VXsol,VYsol,VZsol,TPotential,Potential,FpairX,FpairY,FpairZ,motorCM,random,neighbour_num,&
neighbour,max_neighbour,skinrad,max_solvent_dist,time_count,time_penetrate,outside_particle_num,&
outside_particle,inside_particle_num,inside_particle,sort_step,Freq1,Freq4,data_id,sub_group_id,&
sys_dspace_id,sys_crp_list,sys_dset_id,conc_dspace_id,conc_crp_list,conc_dset_id,solv_dspace_id,&
solv_crp_list,solv_dset_id,dimer_dspace_id,dimer_crp_list,dimer_dset_id)
else
 cstep = 0
 call output_hdf5_initiate(cstep,Freq1,Freq2,Freq3,Freq4,box,T,MDt,MPCt,NT,numdims)
 call initial(Box,numdims,Dcn,Na0,Nb0,NT,cutoff,Rc,Rn,zdimeq,varA,varB,mean,flag,type_B,Xdim,Ydim,Zdim, &
                   VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,Xsol,Ysol,Zsol,VXsol,VYsol,VZsol, &
                   FpairX,FpairY,FpairZ,motorCM,sep)
 ! Sorting algorithm
 if ( mod(cstep,sort_step) == 0 ) then
  call sort_hilbert(vol,NT,box,xsol,ysol,zsol,flag,type_B,vxsol,vysol,vzsol)
 end if
 neighbour_num = 0
 neighbour = 0
 outside_particle_num = 0
 outside_particle = 0
 inside_particle_num = 0
 inside_particle = 0
 do k = 1,numdims
  ! Catalytic sphere
  call neighbour_list(NT,Box,cutoffC,skinrad,Xsol,Ysol,Zsol,Xdim(k,1),Ydim(k,1),Zdim(k,1),neighbour_num(k,1),neighbour(k,1,:))
  ! Noncatalytic sphere
  call neighbour_list(NT,Box,cutoffN,skinrad,Xsol,Ysol,Zsol,Xdim(k,2),Ydim(k,2),Zdim(k,2),neighbour_num(k,2),neighbour(k,2,:))
 end do
 call inout_particle_list(NT,box,cutoffC,cutoffN,skinrad,xsol,ysol,zsol,xdim,ydim,zdim,&
outside_particle_num,outside_particle,inside_particle_num,inside_particle,numdims)
 max_solvent_dist = 0
 time_penetrate = floor((0.5*skinrad / (sqrt(maxval(VXsol*VXsol+VYsol*VYsol+VZsol*VZsol)))) / MDt)
 time_count = 0
 call output_hdf5_write(cstep,Vol,Freq1,Freq2,Freq3,Freq4,NT,numdims,Ma,Mb,Mc,Mn,zdimeq,box,x_width,y_width,&
xsol,ysol,zsol,flag,type_B,vxsol,vysol,vzsol,cutoffC,Dcn,concA1d,concB1d,concA2d,concB2d,pxsol,pysol,pzsol,&
xdim,ydim,zdim,vxdim,vydim,vzdim,Tpotential)
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
   call evolve(cstep,time,NT,MDt,MPCtstep,numdims,Ma,Mb,Mc,Mn,Dcn,Box,cutoff,Rc,Rn,EaC,EbC,EaN,EbN,&
           zdimeq,k_harm,RPab,RPba,sphere,flag,type_B,xsol,Ysol,Zsol,VXsol,VYsol,VZsol,Xdim,&
                  Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,motorCM,TPotential,Potential,FpairX,&
                  FpairY,FpairZ,random,skinrad,max_solvent_dist,neighbour_num,neighbour,time_count,&
                  time_penetrate,outside_particle_num,outside_particle,inside_particle_num,inside_particle,&
                  vol,sort_step,values)
  else if ( mod(cstep,MPCtstep) == 0 ) then
   call evolve(cstep,time,NT,MDt,MPCtstep,numdims,Ma,Mb,Mc,Mn,Dcn,Box,cutoff,Rc,Rn,EaC,EbC,EaN,EbN,&
           zdimeq,k_harm,RPab,RPba,sphere,flag,type_B,xsol,Ysol,Zsol,VXsol,VYsol,VZsol,Xdim,&
                  Ydim,Zdim,VXdim,VYdim,VZdim,FXdim,FYdim,FZdim,motorCM,TPotential,Potential,FpairX,&
                  FpairY,FpairZ,random,skinrad,max_solvent_dist,neighbour_num,neighbour,time_count,&
                  time_penetrate,outside_particle_num,outside_particle,inside_particle_num,inside_particle,&
                  vol,sort_step,values)
call distribute(NT,numdims,cutoff,Rc,Rn,flag,Xsol,Ysol,Zsol,Xdim,Ydim,Zdim,Box,tACellN,tBCellN, &
                      sACellN,sBCellN,cutoffC2,cutoffN2,cstep,outside_particle_num,outside_particle,&
                      inside_particle_num,inside_particle)
   call mpc_rotation(NT,rotang,mpcnumcell,a0,Box,Xsol,Ysol,Zsol,VXsol,VYsol,VZsol,VXCM,VYCM,VZCM, &
                     tCellShift,NewCellN,UVX,UVY,UVZ)
             call reaction_simple(MPCt,rxnnumcell,k1,flag,type_B,tBCellN,sBCellN)
  end if
! Output to file
 call output_hdf5_write(cstep,Vol,Freq1,Freq2,Freq3,Freq4,NT,numdims,Ma,Mb,Mc,Mn,zdimeq,box,x_width,y_width,&
         xsol,ysol,zsol,flag,type_B,vxsol,vysol,vzsol,cutoffC,Dcn,concA1d,concB1d,concA2d,concB2d,pxsol,pysol,pzsol,&
xdim,ydim,zdim,vxdim,vydim,vzdim,Tpotential)

! Timing the simulation
  if ( mod(int(cstep),1000) == 0 ) then
   call date_and_time(VALUES=values)
   open(unit=10,file='Time.dat',status='old',action='write',position='append')
    write(10,*) values(5),values(6),values(7), cstep
   close(10)
  end if

end do

call calc_solv_type_num(NT,Vol,flag,iTotA,iTotB)
open(unit=10,file='timestep.dat',status='replace')
 write(10,*)cstep
 write(10,*)iTotA
 write(10,*)iTotB
close(10)

!Deallocate arrays
deallocate(flag)
deallocate(Xdim,Ydim,Zdim)
deallocate(VXdim,VYdim,VZdim)
deallocate(FXdim,FYdim,FZdim)
deallocate(Xsol,Ysol,Zsol)
deallocate(VXsol,VYsol,VZsol)
deallocate(tACellN,sACellN)
deallocate(tBCellN,sBCellN)
deallocate(Potential,FpairX,FpairY,FpairZ)
deallocate(VXCM,VYCM,VZCM)
deallocate(tCellShift)
deallocate(NewCellN)
deallocate(UVX,UVY,UVZ)
deallocate(PXsol,PYsol,PZsol)
deallocate(random)
deallocate(neighbour,neighbour_num,max_solvent_dist)
deallocate(concA1d,concB1d)
deallocate(concA2d,concB2d)
deallocate(outside_particle)
deallocate(inside_particle)
deallocate(motorCM)
deallocate(type_B)

end program main_rotation
