PROGRAM Main_SelkovDimer
USE Global_SelkovDimer
USE mt95
IMPLICIT NONE

!Read initial conditions
OPEN(UNIT=10,FILE='Input_SelkovDimer',STATUS='OLD',ACTION='READ')
READ(10,*)seed
READ(10,*)Nf0
READ(10,*)Ng0
READ(10,*)Ni0
READ(10,*)Box(1)
READ(10,*)Box(2)
READ(10,*)Box(3)
READ(10,*)Mf
READ(10,*)Mg
READ(10,*)Mi
READ(10,*)Rc
READ(10,*)Rn
READ(10,*)Dcn
READ(10,*)EfC
READ(10,*)EgC
READ(10,*)EiC
READ(10,*)EfN
READ(10,*)EgN
READ(10,*)EiN
READ(10,*)RPfg
READ(10,*)RPgf
READ(10,*)T
READ(10,*)A0
READ(10,*)MDt
READ(10,*)MPCt
READ(10,*)k1
READ(10,*)k2
READ(10,*)k3
READ(10,*)k4
READ(10,*)k5
READ(10,*)k6
READ(10,*)Freq1
READ(10,*)Freq2
READ(10,*)Freq3
READ(10,*)Freq4
READ(10,*)Freq5
READ(10,*)Freq6
CLOSE(10)

NT=Nf0+Ng0+Ni0
space=NT

!Allocate arrays
ALLOCATE(Xdim(1:2),Ydim(1:2),Zdim(1:2))
ALLOCATE(prXdim(1:2),prYdim(1:2),prZdim(1:2))
ALLOCATE(Xsol(1:space),Ysol(1:space),Zsol(1:space))
ALLOCATE(FLAG(1:space))
ALLOCATE(R(1:2,1:3),dR(1:3))
ALLOCATE(VXdim(1:2),VYdim(1:2),VZdim(1:2))
ALLOCATE(prVXdim(1:2),prVYdim(1:2),prVZdim(1:2))
ALLOCATE(VXsol(1:space),VYsol(1:space),VZsol(1:space))
ALLOCATE(FXsol(1:space),FYsol(1:space),FZsol(1:space))
ALLOCATE(FXdim(1:2),FYdim(1:2),FZdim(1:2))
ALLOCATE(pFXsol(1:space),pFYsol(1:space),pFZsol(1:space))
ALLOCATE(pFXdim(1:2),pFYdim(1:2),pFZdim(1:2))
ALLOCATE(tFCellN(1:Vol),tGCellN(1:Vol),sFCellN(1:Vol,1:MaxF),sGCellN(1:Vol,1:MaxG))
ALLOCATE(QXdim(1:2),QYdim(1:2),QZdim(1:2))
ALLOCATE(affectedCellN(1:Vol))
ALLOCATE(PXsol(1:space),PYsol(1:space),PZsol(1:space))
ALLOCATE(VXCM(1:Vol),VYCM(1:Vol),VZCM(1:Vol))
ALLOCATE(tCellshift(1:Vol))
ALLOCATE(NewCellN(1:space))
ALLOCATE(UVX(1:Vol),UVY(1:Vol),UVZ(1:Vol))
ALLOCATE(Potential(1:space))
ALLOCATE(FpairX(1:space,1:2),FpairY(1:space,1:2),FpairZ(1:space,1:2))
ALLOCATE(SolCenterX(1:space),SolCenterY(1:space),SolCenterZ(1:space))
ALLOCATE(Rx(1:space),Ry(1:space),Rz(1:space),R2(1:space))
ALLOCATE(tICellN(1:Vol),sICellN(1:Vol,1:MaxI))
ALLOCATE(newXsol(1:space),newYsol(1:space),newZsol(1:space))
ALLOCATE(newFLAG(1:space))
ALLOCATE(newVXsol(1:space),newVYsol(1:space),newVZsol(1:space))
ALLOCATE(newFXsol(1:space),newFYsol(1:space),newFZsol(1:space))
ALLOCATE(newpFXsol(1:space),newpFYsol(1:space),newpFZsol(1:space))
ALLOCATE(Shift(1:3))
ALLOCATE(Xdimflow(1:2),Ydimflow(1:2),Zdimflow(1:2))
ALLOCATE(Xsolflow(1:space),Ysolflow(1:space),Zsolflow(1:space))
ALLOCATE(XsolrelCflow(1:space),YsolrelCflow(1:space),ZsolrelCflow(1:space))
ALLOCATE(Xsolrotflow(1:space),Ysolrotflow(1:space),Zsolrotflow(1:space))
ALLOCATE(VXsolrotflow(1:space),VYsolrotflow(1:space),VZsolrotflow(1:space))

!Read initial positions, velocities, and forces
MPCtstep=INT(MPCt/MDt)				!Number of time step for MPC calc.
T=T/5.d0					!If T is a fraction
Mc=(NT/DFLOAT(Vol))*((4.*pi)/3)*Rc**3			!Mass of C sphere
Mn=(NT/DFLOAT(Vol))*((4.*pi)/3)*Rn**3			!Mass of N sphere
 CUTOFFc=Rc*CUTOFF				!Cutoff for C sphere LJ pot.
 CUTOFFn=Rn*CUTOFF				!Cutoff for N sphere LJ pot.
 CUTOFFc2=CUTOFFc*CUTOFFc			!Cut. C sphere LJ pot. squared
 CUTOFFn2=CUTOFFn*CUTOFFn			!Cut. N sphere LJ pot. squared
HDcn=0.5d0*Dcn			!Half of internuclear distance
varF=T/Mf			!Variance of F velocity Gaussian distribution
varG=T/Mg			!Variance of G velocity Gaussian distribution
varI=T/Mi			!Variance of inert velocity Gauss. distr.
mean=0.				!Initial mean
DO i=1,3
 HBox(i)=0.5d0*Box(i)				!Center of box
 iBox(i)=INT(Box(i))
END DO
 CALL genrand_init( put=seed )

OPEN(UNIT=10,FILE='timestep.dat')
 READ(10,23)cstep
 READ(10,23)iTotF
 READ(10,23)iTotG
 READ(10,23)iTotI
CLOSE(10)

IF (cstep == 0) THEN
 CALL Initial_SelkovDimer
 cstep=1
ELSE IF (cstep /= 0) THEN
 CALL Restart
END IF

!Define initial values for the current and previous steps

!Main time loop
DO time=1,tstep
  IF (MOD(cstep,MPCtstep) /= 0) THEN
   CALL Motion_SelkovDimer
   CALL ForcesMD_SelkovDimer
  ELSE IF (MOD(cstep,MPCtstep) == 0) THEN
   IF (MOD(cstep,10000) == 0) THEN
    CALL Sort
   END IF
   CALL Motion_SelkovDimer
   CALL Distribute_SelkovDimer
   CALL ForcesMD_SelkovDimer
   CALL Affected_Cell
   CALL RMPC_SelkovDimer
  END IF
 ! CALL Output_SelkovDimer
 ! CALL Output_Concentration_Catalytic
  IF (MOD(cstep,100) == 0) THEN
   CALL itime(now)
   OPEN(UNIT=10,FILE='Time.dat',STATUS='replace',POSITION='APPEND')
    WRITE(10,94) now(1),now(2),now(3), cstep
   CLOSE(10)
   94 FORMAT(3(I2,1X),I9)
  END IF
cstep=cstep+1
END DO

cstep=cstep-1
OPEN(UNIT=10,FILE='timestep.dat',STATUS='REPLACE')
 WRITE(10,23)cstep
 WRITE(10,23)iTotF
 WRITE(10,23)iTotG
 WRITE(10,23)iTotI
CLOSE(10)

23 FORMAT(I9)

!Deallocate arrays
DEALLOCATE(Xdim,Ydim,Zdim)
DEALLOCATE(prXdim,prYdim,prZdim)
DEALLOCATE(Xsol,Ysol,Zsol)
DEALLOCATE(FLAG)
DEALLOCATE(R,dR)
DEALLOCATE(VXdim,VYdim,VZdim)
DEALLOCATE(prVXdim,prVYdim,prVZdim)
DEALLOCATE(VXsol,VYsol,VZsol)
DEALLOCATE(FXsol,FYsol,FZsol)
DEALLOCATE(FXdim,FYdim,FZdim)
DEALLOCATE(pFXsol,pFYsol,pFZsol)
DEALLOCATE(pFXdim,pFYdim,pFZdim)
DEALLOCATE(tFCellN,tGCellN,sFCellN,sGCellN)
DEALLOCATE(QXdim,QYdim,QZdim)
DEALLOCATE(affectedCellN)
DEALLOCATE(PXsol,PYsol,PZsol)
DEALLOCATE(VXCM,VYCM,VZCM)
DEALLOCATE(tCellshift)
DEALLOCATE(NewCellN)
DEALLOCATE(UVX,UVY,UVZ)
DEALLOCATE(Potential)
DEALLOCATE(FpairX,FpairY,FpairZ)
DEALLOCATE(SolCenterX,SolCenterY,SolCenterZ)
DEALLOCATE(Rx,Ry,Rz,R2)
DEALLOCATE(tICellN,sICellN)
DEALLOCATE(newXsol,newYsol,newZsol)
DEALLOCATE(newFLAG)
DEALLOCATE(newVXsol,newVYsol,newVZsol)
DEALLOCATE(newFXsol,newFYsol,newFZsol)
DEALLOCATE(newpFXsol,newpFYsol,newpFZsol)
DEALLOCATE(Shift)
DEALLOCATE(Xdimflow,Ydimflow,Zdimflow)
DEALLOCATE(Xsolflow,Ysolflow,Zsolflow)
DEALLOCATE(XsolrelCflow,YsolrelCflow,ZsolrelCflow)
DEALLOCATE(Xsolrotflow,Ysolrotflow,Zsolrotflow)
DEALLOCATE(VXsolrotflow,VYsolrotflow,VZsolrotflow)

END PROGRAM Main_SelkovDimer
