SUBROUTINE Output_Concentration_Catalytic
USE Global_SelkovDimer
IMPLICIT NONE

!This is a subroutine to calculate the concentration in concentric spheres around the dimer's center of mass. The inner sphere is of radius 10, the first shell is between radii 10 and 15, and the outer shell is between radii 15 and 20.

IF (MOD(cstep,Freq5) == 0) THEN

!!!!!!Set initial values
Tot3F=0.
Tot3G=0.
Tot3I=0.
Tot4F=0.
Tot4G=0.
Tot4I=0.
Tot6F=0.
Tot6G=0.
Tot6I=0.
Tot8F=0.
Tot8G=0.
Tot8I=0.
Tot10F=0.
Tot10G=0.
Tot10I=0.
Tot12F=0.
Tot12G=0.
Tot12I=0.
Tot15F=0.
Tot15G=0.
Tot15I=0.
Tot20F=0.
Tot20G=0.
Tot20I=0.
Vol2=(4.d0/3.d0)*pi*(2.d0)**3
Vol3=(4.d0/3.d0)*pi*(3.d0)**3
Vol4=(4.d0/3.d0)*pi*(4.d0)**3
Vol6=(4.d0/3.d0)*pi*(6.d0)**3
Vol8=(4.d0/3.d0)*pi*(8.d0)**3
Vol10=(4.d0/3.d0)*pi*(10.d0)**3
Vol12=(4.d0/3.d0)*pi*(12.d0)**3
Vol15=(4.d0/3.d0)*pi*(15.d0)**3
Vol20=(4.d0/3.d0)*pi*(20.d0)**3
ShellVol3=Vol3-Vol2
ShellVol4=Vol4-Vol3-6.02139
ShellVol6=Vol6-Vol4-53.55668
ShellVol8=Vol8-Vol6-98.4366
ShellVol10=Vol10-Vol8-93.051
ShellVol12=Vol12-Vol10-17.0169
ShellVol15=Vol15-Vol12
ShellVol20=Vol20-Vol15

!!!!!!COM calc., translation of COM to center, and dimer/solvent translated accordingly
moveX=HBox(1)-Xdim(1)
moveY=HBox(2)-Ydim(1)
moveZ=HBox(3)-Zdim(1)

!$OMP PARALLEL DO SHARED(SolCenterX,Xsol,SolCenterY,Ysol,SolCenterZ) &
!$OMP SHARED(Zsol,Box,HBox,Rx,Ry,Rz,R,moveX,moveY,moveZ) &
!$OMP SCHEDULE(static)

DO i=1,NT
 IF(FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN

  SolCenterX(i)=Xsol(i)+moveX
  IF (SolCenterX(i) > Box(1))THEN
   SolCenterX(i)=SolCenterX(i)-Box(1)
  ELSE IF (SolCenterX(i) < 0.)THEN
   SolCenterX(i)=SolCenterX(i)+Box(1)
  END IF
  SolCenterY(i)=Ysol(i)+moveY
  IF (SolCenterY(i) > Box(2))THEN
   SolCenterY(i)=SolCenterY(i)-Box(2)
  ELSE IF (SolCenterY(i) < 0.)THEN
   SolCenterY(i)=SolCenterY(i)+Box(2)
  END IF
  SolCenterZ(i)=Zsol(i)+moveZ
  IF (SolCenterZ(i) > Box(3))THEN
   SolCenterZ(i)=SolCenterZ(i)-Box(3)
  ELSE IF (SolCenterZ(i) < 0.)THEN
   SolCenterZ(i)=SolCenterZ(i)+Box(3)
  END IF

  Rx(i)=SolCenterX(i)-HBox(1)
  Ry(i)=SolCenterY(i)-HBox(2)
  Rz(i)=SolCenterZ(i)-HBox(3)
  R2(i)=Rx(i)*Rx(i)+Ry(i)*Ry(i)+Rz(i)*Rz(i)

 END IF
END DO

!$OMP END PARALLEL DO

DO i=1,NT
 IF (FLAG(i) == 0) THEN
  IF (R2(i) >= 4.) THEN
   IF (R2(i) < 9.) THEN
     Tot3F=Tot3F+1.
   ELSE IF (R2(i) >= 9.0 .AND. R2(i) < 16.0)THEN
     Tot4F=Tot4F+1.
   ELSE IF (R2(i) >= 16.0 .AND. R2(i) < 36.0)THEN
     Tot6F=Tot6F+1.
   ELSE IF (R2(i) >= 36.0 .AND. R2(i) < 64.0)THEN
     Tot8F=Tot8F+1.
   ELSE IF (R2(i) >= 64.0 .AND. R2(i) < 100.0)THEN
     Tot10F=Tot10F+1.
   ELSE IF (R2(i) >= 100.0 .AND. R2(i) < 144.0)THEN
     Tot12F=Tot12F+1.
   ELSE IF (R2(i) >= 144.0 .AND. R2(i) < 225.0)THEN
     Tot15F=Tot15F+1.
   ELSE IF (R2(i) >= 225.0 .AND. R2(i) < 400.0)THEN
     Tot20F=Tot20F+1.
   END IF
  END IF
 ELSE IF (FLAG(i) == 1) THEN
  IF (R2(i) >= 4.) THEN
   IF (R2(i) < 9.)THEN
     Tot3G=Tot3G+1.
   ELSE IF (R2(i) >= 9.0 .AND. R2(i) < 16.0)THEN
     Tot4G=Tot4G+1.
   ELSE IF (R2(i) >= 16.0 .AND. R2(i) < 36.0)THEN
     Tot6G=Tot6G+1.
   ELSE IF (R2(i) >= 36.0 .AND. R2(i) < 64.0)THEN
     Tot8G=Tot8G+1.
   ELSE IF (R2(i) >= 64.0 .AND. R2(i) < 100.0)THEN
     Tot10G=Tot10G+1.
   ELSE IF (R2(i) >= 100.0 .AND. R2(i) < 144.0)THEN
     Tot12G=Tot12G+1.
   ELSE IF (R2(i) >= 144.0 .AND. R2(i) < 225.0)THEN
     Tot15G=Tot15G+1.
   ELSE IF (R2(i) >= 225.0 .AND. R2(i) < 400.0)THEN
     Tot20G=Tot20G+1.
   END IF
  END IF
 ELSE IF (FLAG(i) == 3) THEN
  IF (R2(i) >= 4.) THEN
   IF (R2(i) < 9.)THEN
     Tot3I=Tot3I+1.
   ELSE IF (R2(i) >= 9.0 .AND. R2(i) < 16.0)THEN
     Tot4I=Tot4I+1.
   ELSE IF (R2(i) >= 16.0 .AND. R2(i) < 36.0)THEN
     Tot6I=Tot6I+1.
   ELSE IF (R2(i) >= 36.0 .AND. R2(i) < 64.0)THEN
     Tot8I=Tot8I+1.
   ELSE IF (R2(i) >= 64.0 .AND. R2(i) < 100.0)THEN
     Tot10I=Tot10I+1.
   ELSE IF (R2(i) >= 100.0 .AND. R2(i) < 144.0)THEN
     Tot12I=Tot12I+1.
   ELSE IF (R2(i) >= 144.0 .AND. R2(i) < 225.0)THEN
     Tot15I=Tot15I+1.
   ELSE IF (R2(i) >= 225.0 .AND. R2(i) < 400.0)THEN
     Tot20I=Tot20I+1.
   END IF
  END IF
 END IF
END DO

Conc3F=Tot3F/ShellVol3
Conc3G=Tot3G/ShellVol3
Conc3I=Tot3I/ShellVol3
Conc4F=Tot4F/ShellVol4
Conc4G=Tot4G/ShellVol4
Conc4I=Tot4I/ShellVol4
Conc6F=Tot6F/ShellVol6
Conc6G=Tot6G/ShellVol6
Conc6I=Tot6I/ShellVol6
Conc8F=Tot8F/ShellVol8
Conc8G=Tot8G/ShellVol8
Conc8I=Tot8I/ShellVol8
Conc10F=Tot10F/ShellVol10
Conc10G=Tot10G/ShellVol10
Conc10I=Tot10I/ShellVol10
Conc12F=Tot12F/ShellVol12
Conc12G=Tot12G/ShellVol12
Conc12I=Tot12I/ShellVol12
Conc15F=Tot15F/ShellVol15
Conc15G=Tot15G/ShellVol15
Conc15I=Tot15I/ShellVol15
Conc20F=Tot20F/ShellVol20
Conc20G=Tot20G/ShellVol20
Conc20I=Tot20I/ShellVol20

OPEN(UNIT=40, FILE='ShellConcentration_Catalytic.dat', STATUS='replace',ACTION='WRITE',POSITION='APPEND')
 WRITE(40,22)cstep,Conc3F,Conc3G,Conc3I,Conc4F,Conc4G,Conc4I,Conc6F,Conc6G,Conc6I,Conc8F,Conc8G,Conc8I,Conc10F, &
             Conc10G,Conc10I,Conc12F,Conc12G,Conc12I,Conc15F,Conc15G,Conc15I,Conc20F,Conc20G,Conc20I
CLOSE(40)

22 FORMAT(I8,24(1X,ES23.16))

END IF

END SUBROUTINE Output_Concentration_Catalytic
