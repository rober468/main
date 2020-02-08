SUBROUTINE Output_SelkovDimer
USE Global_SelkovDimer
IMPLICIT NONE

!!!!!!Record dimer positions, velocities, and forces at cstep
IF (MOD(cstep,Freq1) == 0) THEN

 OPEN(UNIT=10,FILE='DimerPosC.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 OPEN(UNIT=20,FILE='DimerPosN.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 OPEN(UNIT=30,FILE='DimerVelC.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 OPEN(UNIT=40,FILE='DimerVelN.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 OPEN(UNIT=31,FILE='DimerForC.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 OPEN(UNIT=32,FILE='DimerForN.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 WRITE(10,15) cstep,Xdim(1),Ydim(1),Zdim(1)
 WRITE(20,15) cstep,Xdim(2),Ydim(2),Zdim(2)
 WRITE(30,15) cstep,VXdim(1),VYdim(1),VZdim(1)
 WRITE(40,15) cstep,VXdim(2),VYdim(2),VZdim(2)
 WRITE(31,15) cstep,FXdim(1),FYdim(1),FZdim(1)
 WRITE(32,15) cstep,FXdim(2),FYdim(2),FZdim(2)
 CLOSE(10)
 CLOSE(20)
 CLOSE(30)
 CLOSE(40)
 CLOSE(31)
 CLOSE(32)
15 FORMAT(I8,1X,3(E23.16,1X))

END IF

!!!!!!Record solvent F and G positions and velocities
IF (MOD(cstep,Freq2) == 0) THEN

 IF (cstep > 0 .AND. cstep < 100) THEN
  WRITE(SolventFPosFile,"(A11,I2,A4)") 'SolventFPos',cstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I2,A4)") 'SolventGPos',cstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I2,A4)") 'SolventIPos',cstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I2,A4)") 'SolventFVel',cstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I2,A4)") 'SolventGVel',cstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I2,A4)") 'SolventIVel',cstep,'.dat'
  WRITE(SolventFForFile,"(A11,I2,A4)") 'SolventFFor',cstep,'.dat'
  WRITE(SolventGForFile,"(A11,I2,A4)") 'SolventGFor',cstep,'.dat'
  WRITE(SolventIForFile,"(A11,I2,A4)") 'SolventIFor',cstep,'.dat'
 ELSE IF (cstep >= 100 .AND. cstep < 1000) THEN
  WRITE(SolventFPosFile,"(A11,I3,A4)") 'SolventFPos',cstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I3,A4)") 'SolventGPos',cstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I3,A4)") 'SolventIPos',cstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I3,A4)") 'SolventFVel',cstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I3,A4)") 'SolventGVel',cstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I3,A4)") 'SolventIVel',cstep,'.dat'
  WRITE(SolventFForFile,"(A11,I3,A4)") 'SolventFFor',cstep,'.dat'
  WRITE(SolventGForFile,"(A11,I3,A4)") 'SolventGFor',cstep,'.dat'
  WRITE(SolventIForFile,"(A11,I3,A4)") 'SolventIFor',cstep,'.dat'
 ELSE IF (cstep >= 1000 .AND. cstep < 10000) THEN
  WRITE(SolventFPosFile,"(A11,I4,A4)") 'SolventFPos',cstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I4,A4)") 'SolventGPos',cstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I4,A4)") 'SolventIPos',cstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I4,A4)") 'SolventFVel',cstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I4,A4)") 'SolventGVel',cstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I4,A4)") 'SolventIVel',cstep,'.dat'
  WRITE(SolventFForFile,"(A11,I4,A4)") 'SolventFFor',cstep,'.dat'
  WRITE(SolventGForFile,"(A11,I4,A4)") 'SolventGFor',cstep,'.dat'
  WRITE(SolventIForFile,"(A11,I4,A4)") 'SolventIFor',cstep,'.dat'
 ELSE IF (cstep >= 10000 .AND. cstep < 100000) THEN
  WRITE(SolventFPosFile,"(A11,I5,A4)") 'SolventFPos',cstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I5,A4)") 'SolventGPos',cstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I5,A4)") 'SolventIPos',cstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I5,A4)") 'SolventFVel',cstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I5,A4)") 'SolventGVel',cstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I5,A4)") 'SolventIVel',cstep,'.dat'
  WRITE(SolventFForFile,"(A11,I5,A4)") 'SolventFFor',cstep,'.dat'
  WRITE(SolventGForFile,"(A11,I5,A4)") 'SolventGFor',cstep,'.dat'
  WRITE(SolventIForFile,"(A11,I5,A4)") 'SolventIFor',cstep,'.dat'
 ELSE IF (cstep >= 100000 .AND. cstep < 1000000) THEN
  WRITE(SolventFPosFile,"(A11,I6,A4)") 'SolventFPos',cstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I6,A4)") 'SolventGPos',cstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I6,A4)") 'SolventIPos',cstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I6,A4)") 'SolventFVel',cstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I6,A4)") 'SolventGVel',cstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I6,A4)") 'SolventIVel',cstep,'.dat'
  WRITE(SolventFForFile,"(A11,I6,A4)") 'SolventFFor',cstep,'.dat'
  WRITE(SolventGForFile,"(A11,I6,A4)") 'SolventGFor',cstep,'.dat'
  WRITE(SolventIForFile,"(A11,I6,A4)") 'SolventIFor',cstep,'.dat'
 ELSE IF (cstep >= 1000000 .AND. cstep < 10000000) THEN
  WRITE(SolventFPosFile,"(A11,I7,A4)") 'SolventFPos',cstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I7,A4)") 'SolventGPos',cstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I7,A4)") 'SolventIPos',cstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I7,A4)") 'SolventFVel',cstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I7,A4)") 'SolventGVel',cstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I7,A4)") 'SolventIVel',cstep,'.dat'
  WRITE(SolventFForFile,"(A11,I7,A4)") 'SolventFFor',cstep,'.dat'
  WRITE(SolventGForFile,"(A11,I7,A4)") 'SolventGFor',cstep,'.dat'
  WRITE(SolventIForFile,"(A11,I7,A4)") 'SolventIFor',cstep,'.dat'
 ELSE IF (cstep >= 10000000 .AND. cstep < 100000000) THEN
  WRITE(SolventFPosFile,"(A11,I8,A4)") 'SolventFPos',cstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I8,A4)") 'SolventGPos',cstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I8,A4)") 'SolventIPos',cstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I8,A4)") 'SolventFVel',cstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I8,A4)") 'SolventGVel',cstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I8,A4)") 'SolventIVel',cstep,'.dat'
  WRITE(SolventFForFile,"(A11,I8,A4)") 'SolventFFor',cstep,'.dat'
  WRITE(SolventGForFile,"(A11,I8,A4)") 'SolventGFor',cstep,'.dat'
  WRITE(SolventIForFile,"(A11,I8,A4)") 'SolventIFor',cstep,'.dat'
 END IF

 OPEN(UNIT=10,FILE=SolventFPosFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=20,FILE=SolventGPosFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=30,FILE=SolventIPosFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=40,FILE=SolventFVelFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=31,FILE=SolventGVelFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=32,FILE=SolventIVelFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=33,FILE=SolventFForFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=34,FILE=SolventGForFile,STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=35,FILE=SolventIForFile,STATUS='REPLACE',ACTION='WRITE')
 DO i=1,NT
  IF (FLAG(i) == 0) THEN
   WRITE(10,16) Xsol(i),Ysol(i),Zsol(i), i
   WRITE(40,16) VXsol(i),VYsol(i),VZsol(i), i
   WRITE(33,16) FXsol(i),FYsol(i),FZsol(i), i
  ELSE IF (FLAG(i) == 1) THEN
   WRITE(20,16) Xsol(i),Ysol(i),Zsol(i), i
   WRITE(31,16) VXsol(i),VYsol(i),VZsol(i), i
   WRITE(34,16) FXsol(i),FYsol(i),FZsol(i), i
  ELSE IF (FLAG(i) == 3) THEN
   WRITE(30,16) Xsol(i),Ysol(i),Zsol(i), i
   WRITE(32,16) VXsol(i),VYsol(i),VZsol(i), i
   WRITE(35,16) FXsol(i),FYsol(i),FZsol(i), i
  END IF
 END DO
 CLOSE(10)
 CLOSE(20)
 CLOSE(30)
 CLOSE(40)
 CLOSE(31)
 CLOSE(32)
 CLOSE(33)
 CLOSE(34)
 CLOSE(35)

16 FORMAT(3(ES23.16,1X),I8)

END IF

!!!!!!Calculation of velocity along internuclear vector
IF(MOD(cstep,Freq4) == 0) THEN

 !Difference in dimer positions
 dXdim=Xdim(1)-Xdim(2)
 dYdim=Ydim(1)-Ydim(2)
 dZdim=Zdim(1)-Zdim(2)
 !Apply periodic conditions
 IF (dXdim > HBox(1)) THEN
  dXdim=dXdim-Box(1)
 ELSE IF (dXdim < -HBox(1)) THEN
  dXdim=dXdim+Box(1)
 END IF
 IF (dYdim > HBox(2)) THEN
  dYdim=dYdim-Box(2)
 ELSE IF (dYdim < -HBox(2)) THEN
  dYdim=dYdim+Box(2)
 END IF
 IF (dZdim > HBox(3)) THEN
  dZdim=dZdim-Box(3)
 ELSE IF (dZdim < -HBox(3)) THEN
  dZdim=dZdim+Box(3)
 END IF
 
 !Find distance between the two dimers and make unit vector
 dXYZ=SQRT(dXdim*dXdim+dYdim*dYdim+dZdim*dZdim)
 NucVecX=dXdim/dXYZ
 NucVecY=dYdim/dXYZ
 NucVecZ=dZdim/dXYZ

 VXdimav=(Mc*VXdim(1)+Mn*VXdim(2))/(Mn+Mc)
 VYdimav=(Mc*VYdim(1)+Mn*VYdim(2))/(Mn+Mc)
 VZdimav=(Mc*VZdim(1)+Mn*VZdim(2))/(Mn+Mc)

 DotVZ=VXdimav*NucVecX+VYdimav*NucVecY+VZdimav*NucVecZ

 !Record cstep and DotVZ
 OPEN(UNIT=10,FILE='DotVZ.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
  WRITE(10,18) cstep, DotVZ, NucVecX, NucVecY, NucVecZ
 CLOSE(10)

18 FORMAT(I8,4(1X,ES23.16))

END IF

!!!!!!Concentration of F and G
IF (MOD(cstep,Freq5) == 0) THEN

 iTotF=0
 iTotG=0
 iTotI=0
 TotF=0.
 TotG=0.
 TotI=0.

 DO i=1,NT
  IF (FLAG(i) == 0) THEN
   TotF=TotF+1.
   iTotF=iTotF+1
  ELSE IF (FLAG(i) == 1) THEN
   TotG=TotG+1.
   iTotG=iTotG+1
  ELSE IF (FLAG(i) == 3) THEN
   TotI=TotI+1.
   iTotI=iTotI+1
  END IF
 END DO

 ConcF=TotF/Vol
 ConcG=TotG/Vol
 ConcI=TotI/Vol

 tFBulk=0.
 tGBulk=0.
 tIBulk=0.
 TBulkCells=0.

 DO i=1,Vol
  IF(affectedCellN(i) /=1) THEN
   tFBulk=tFBulk+tFCellN(i)
   tGBulk=tGBulk+tGCellN(i)
   tIBulk=tIBulk+tICellN(i)
   TBulkCells=TBulkCells+1
  END IF
 END DO

 ConcFBulk=tFBulk/TBulkCells
 ConcGBulk=tGBulk/TBulkCells
 ConcIBulk=tIBulk/TBulkCells

 !Record concentrations of F and G
 OPEN(UNIT=10,FILE='Concentration.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
  WRITE(10,19) cstep, TotF, TotG, TotI, ConcF, ConcG, ConcI, ConcFBulk, ConcGBulk, ConcIBulk, TBulkCells
 CLOSE(10)

19 FORMAT(I8,10(1X,ES23.16))
END IF

!!!!!!Calculation of momentum and energy in system
IF (MOD(cstep,Freq3) == 0) THEN

 sumPXsol=0.d0
 sumPYsol=0.d0
 sumPZsol=0.d0
 sumPX=0.d0
 sumPY=0.d0
 sumPZ=0.d0
 KEsol=0.d0
 KEdim=0.d0 

 !Momentum calculation
 DO i=1,NT,1
  IF (FLAG(i) == 0) THEN
   PXsol(i)=Mf*VXsol(i)
   PYsol(i)=Mf*VYsol(i)
   PZsol(i)=Mf*VZsol(i)
  ELSE IF (FLAG(i) == 1) THEN
   PXsol(i)=Mg*VXsol(i)
   PYsol(i)=Mg*VYsol(i)
   PZsol(i)=Mg*VZsol(i)
  ELSE IF (FLAG(i) == 3) THEN
   PXsol(i)=Mi*VXsol(i)
   PYsol(i)=Mi*VYsol(i)
   PZsol(i)=Mi*VZsol(i)
  END IF
  IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
   sumPXsol=sumPXsol+PXsol(i)
   sumPYsol=sumPYsol+PYsol(i)
   sumPZsol=sumPZsol+PZsol(i)
  END IF
 END DO

 PXdim=(Mc+Mn)*Vxdimav
 PYdim=(Mc+Mn)*Vydimav
 PZdim=(Mc+Mn)*Vzdimav

 sumPX=sumPXsol+PXdim
 sumPY=sumPYsol+PYdim
 sumPZ=sumPZsol+PZdim

 !Energy Calculation
 DO i=1,NT,1
  IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
   KEsol=KEsol+0.5*(VXsol(i)*PXsol(i)+VYsol(i)*PYsol(i)+VZsol(i)*PZsol(i))
  END IF
 END DO
 KEdim=KEdim+0.5*(1/(Mc+Mn))*(PXdim*PXdim+PYdim*PYdim+PZdim*PZdim)

 Ekin=KEdim+KEsol
 Epot=TPotential
 Etot=Ekin+Epot

!Temperature Calculation
tempsol=(KEsol/(TotF+TotG+TotI))*(2./3.)
tempdim=KEdim*(2./3.)


 !Record cstep, sumPX, sum PY, sumPZ, KEsol, Ekin, Epot, Etot
 OPEN(UNIT=10,FILE='EnergyMomentum.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
  WRITE(10,17) cstep, sumPX, sumPY, sumPZ, Ekin, Epot, Etot, tempsol, tempdim
 CLOSE(10)

17 FORMAT(I8,8(1X,ES23.16))

END IF

!!!!!!Calculations for flow fields
!This converts the positions and velocities of the dimer and solvent particles so that the dimer is in the center. The data is converted so that the flow is visualized in a 2D slab, where the dimer's rigid bond is in the plane.
!|...................|   O = C sphere
!|.........__........|   __
!|......../  \.......|  /  \
!|.....O--|  |.......|  |  | = N sphere
!|........\__/.......|  \__/
!|...................|
!|___________________|  . = solvent

IF (MOD(cstep,Freq6) == 0) THEN

!Figure out translation distance in each direction
Shift(1)=(hBox(1)-0.5*Dcn)-Xdim(1)
Shift(2)=hBox(2)-Ydim(1)
Shift(3)=hBox(3)-Zdim(1)

!Apply translation to dimer and solvent (with periodic boundary conditions)
DO i=1,2
 Xdimflow(i)=Xdim(i)+Shift(1)
 Ydimflow(i)=Ydim(i)+Shift(2)
 Zdimflow(i)=Zdim(i)+Shift(3)
END DO
IF (Xdimflow(2) > Box(1)) THEN
 Xdimflow(2)=Xdimflow(2)-Box(1)
ELSE IF (Xdimflow(2) < 0.) THEN
 Xdimflow(2)=Xdimflow(2)+Box(1)
END IF
IF (Ydimflow(2) > Box(2)) THEN
 Ydimflow(2)=Ydimflow(2)-Box(2)
ELSE IF (Ydimflow(2) < 0.) THEN
 Ydimflow(2)=Ydimflow(2)+Box(2)
END IF
IF (Zdimflow(2) > Box(3)) THEN
 Zdimflow(2)=Zdimflow(2)-Box(3)
ELSE IF (Zdimflow(2) < 0.) THEN
 Zdimflow(2)=Zdimflow(2)+Box(3)
END IF

!$OMP PARALLEL DO &
!$OMP SHARED (Xsolflow,Ysolflow,Zsolflow,Shift,Box,FLAG) &
!$OMP PRIVATE (i)

DO i=1,NT
 IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
  Xsolflow(i)=Xsol(i)+Shift(1)
  Ysolflow(i)=Ysol(i)+Shift(2)
  Zsolflow(i)=Zsol(i)+Shift(3)
  IF (Xsolflow(i) > Box(1)) THEN
   Xsolflow(i)=Xsolflow(i)-Box(1)
  ELSE IF (Xsolflow(i) < 0.) THEN
   Xsolflow(i)=Xsolflow(i)+Box(1)
  END IF
  IF (Ysolflow(i) > Box(2)) THEN
   Ysolflow(i)=Ysolflow(i)-Box(2)
  ELSE IF (Ysolflow(i) < 0.) THEN
   Ysolflow(i)=Ysolflow(i)+Box(2)
  END IF
  IF (Zsolflow(i) > Box(3)) THEN
   Zsolflow(i)=Zsolflow(i)-Box(3)
  ELSE IF (Zsolflow(i) < 0.) THEN
   Zsolflow(i)=Zsolflow(i)+Box(3)
  END IF
 END IF
END DO

!$OMP END PARALLEL DO

!!!!!!Rotate everything around C dimer
!Figure out rotation matrix for polar angle and apply
XdimNrelCflow=Xdimflow(2)-Xdimflow(1)
YdimNrelCflow=Ydimflow(2)-Ydimflow(1)
RdimNrelCflow=SQRT(XdimNrelCflow**2+YdimNrelCflow**2)
IRdimNrelCflow=1./RdimNrelCflow
thetaflow=DACOS(XdimNrelCflow*IRdimNrelCflow)
IF (YdimNrelCflow > 0.) THEN
 thetaflow=thetaflow
ELSE IF (YdimNrelCflow < 0.) THEN
 thetaflow=-1*thetaflow
END IF
Xdimrotflow=DCOS(thetaflow)*XdimNrelCflow+DSIN(thetaflow)*YdimNrelCflow
Ydimrotflow=-1*DSIN(thetaflow)*XdimNrelCflow+DCOS(thetaflow)*YdimNrelCflow
Xdimflow(2)=Xdimflow(1)+Xdimrotflow
Ydimflow(2)=Ydimflow(1)+Ydimrotflow

!$OMP PARALLEL DO &
!$OMP SHARED (XsolrelCflow,YsolrelCflow,thetaflow,Xsolflow,Ysolflow) &
!$OMP SHARED (Xdimflow,Ydimflow,Xsolrotflow,Ysolrotflow) &
!$OMP SHARED (VXsolrotflow,VYsolrotflow,VXsol,VYsol,FLAG) &
!$OMP PRIVATE (i)

DO i=1,NT
 IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
  XsolrelCflow(i)=Xsolflow(i)-Xdimflow(1)
  YsolrelCflow(i)=Ysolflow(i)-Ydimflow(1)
  Xsolrotflow(i)=DCOS(thetaflow)*XsolrelCflow(i)+DSIN(thetaflow)*YsolrelCflow(i)
  Ysolrotflow(i)=-1*DSIN(thetaflow)*XsolrelCflow(i)+DCOS(thetaflow)*YsolrelCflow(i)
  Xsolflow(i)=Xdimflow(1)+Xsolrotflow(i)
  Ysolflow(i)=Ydimflow(1)+Ysolrotflow(i)
  VXsolrotflow(i)=DCOS(thetaflow)*VXsol(i)+DSIN(thetaflow)*VYsol(i)
  VYsolrotflow(i)=-1*DSIN(thetaflow)*VXsol(i)+DCOS(thetaflow)*VYsol(i)
 END IF
END DO

!$OMP END PARALLEL DO

!Figure out rotation matrix for azimuthal angle and apply
XdimNrelCflow=Xdimflow(2)-Xdimflow(1)
ZdimNrelCflow=Zdimflow(2)-Zdimflow(1)
RdimNrelCflow=SQRT(XdimNrelCflow**2+ZdimNrelCflow**2)
IRdimNrelCflow=1./RdimNrelCflow
phiflow=DACOS(XdimNrelCflow*IRdimNrelCflow)
IF (ZdimNrelCflow > 0) THEN
 phiflow=phiflow
ELSE IF (ZdimNrelCflow < 0) THEN
 phiflow=-1*phiflow
END IF
Xdimrotflow=DCOS(phiflow)*XdimNrelCflow+DSIN(phiflow)*ZdimNrelCflow
Zdimrotflow=-1*DSIN(phiflow)*XdimNrelCflow+DCOS(phiflow)*ZdimNrelCflow
Xdimflow(2)=Xdimflow(1)+Xdimrotflow
Zdimflow(2)=Zdimflow(1)+Zdimrotflow

!$OMP PARALLEL DO &
!$OMP SHARED (XsolrelCflow,ZsolrelCflow,phiflow,Xsolflow,Zsolflow) &
!$OMP SHARED (Xdimflow,Zdimflow,Xsolrotflow,Zsolrotflow) &
!$OMP SHARED (VXsolrotflow,VZsolrotflow,VZsol,FLAG) &
!$OMP PRIVATE (i)

DO i=1,NT
 IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
  XsolrelCflow(i)=Xsolflow(i)-Xdimflow(1)
  ZsolrelCflow(i)=Zsolflow(i)-Zdimflow(1)
  Xsolrotflow(i)=DCOS(phiflow)*XsolrelCflow(i)+DSIN(phiflow)*ZsolrelCflow(i)
  Zsolrotflow(i)=-1*DSIN(phiflow)*XsolrelCflow(i)+DCOS(phiflow)*ZsolrelCflow(i)
  Xsolflow(i)=Xdimflow(1)+Xsolrotflow(i)
  Zsolflow(i)=Zdimflow(1)+Zsolrotflow(i)
  VXsolrotflow(i)=DCOS(phiflow)*VXsolrotflow(i)+DSIN(phiflow)*VZsol(i)
  VZsolrotflow(i)=-1*DSIN(phiflow)*VXsolrotflow(i)+DCOS(phiflow)*VZsol(i)
 END IF
END DO

!$OMP END PARALLEL DO

!!!!!!Calculate average VXsolrot and VYsolrot in each cell of meshgrid
ALLOCATE(CellVelX(1:iBox(1),1:iBox(2)),CellVelY(1:iBox(1),1:iBox(2)),n(1:iBox(1),1:iBox(2)))

DO i=1,iBox(1)
 DO j=1,iBox(2)
  CellVelX(i,j)=0.
  CellVelY(i,j)=0.
  n(i,j)=0.
 END DO
END DO
DO i=1,NT
 IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
  IF (Zsolflow(i) >= 19.5 .AND. Zsolflow(i) < 20.5) THEN
   DO j=1,iBox(1)
    IF (Xsolflow(i) >= j-1.  .AND. Xsolflow(i) < j) THEN 
     DO k=1,iBox(2)
      IF (Ysolflow(i) >= k-1. .AND. Ysolflow(i) < k) THEN
       CellVelX(j,k)=CellVelX(j,k)+VXsolrotflow(i)
       CellVelY(j,k)=CellVelY(j,k)+VYsolrotflow(i)
       n(j,k)=n(j,k)+1.
      END IF
     END DO
    END IF
   END DO
  END IF
 END IF
END DO
DO j=1,iBox(1)
 DO k=1,iBox(2)
  IF (n(j,k) > 0.5) THEN
   CellVelX(j,k)=CellVelX(j,k)/n(j,k)
   CellVelY(j,k)=CellVelY(j,k)/n(j,k)
  ELSE IF (n(j,k) < 0.5) THEN
   CellVelX(j,k)=0.
   CellVelY(j,k)=0.
  END IF
 END DO
END DO

!Save velocity field in 2D slab
OPEN(UNIT=10,FILE='CellVelAvX.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
OPEN(UNIT=20,FILE='CellVelAvY.dat',STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 DO i=1,iBox(1)
  WRITE(10,*) (CellVelX(i,j), j=1,iBox(2))
  WRITE(20,*) (CellVelY(i,j), j=1,iBox(2))
 END DO
CLOSE(10)
CLOSE(20)

DEALLOCATE(CellVelX,CellVelY,n)

END IF

END SUBROUTINE Output_SelkovDimer
