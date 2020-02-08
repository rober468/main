SUBROUTINE RMPC_SelkovDimer
USE Global_SelkovDimer
USE mt95
IMPLICIT NONE

!ThE first part of this subroutine is the MPC calculation.  This calculation assumes that the masses of the solvent F and G are the SAME (cannot be used for solvent particles of different masses).  The second part of this subroutine takes into account the reaction in the bulk phase of the solution.


!*********************************************************************************!
!******************************MPC************************************************!
!*********************************************************************************!

!!!!!!Shift particles between [-a0/2,a0/2] and calculate CM pre-collision velocity
!Set CM pre-collision velocities and total particles in cell N to 0
DO i=1,Vol
 VXCM(i)=0.
 VYCM(i)=0.
 VZCM(i)=0.
 tCellshift(i)=0
END DO

!Calculate random shift between [-a0/2,a0/2]
 CALL genrand_real1(Xshift)
Xshift=A0*(Xshift-0.5d0)			!Rand. # betw. [-a0/2,a0/2], x
 CALL genrand_real1(Yshift)
Yshift=A0*(Yshift-0.5d0)			!Rand. # betw. [-a0/2,a0/2], y
 CALL genrand_real1(Zshift)
Zshift=A0*(Zshift-0.5d0)			!Rand. # betw. [-a0/2,a0/2], z

!Apply shift and PBC to ith solvent particle
DO i=1,NT
 IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN

  Xsolshift=Xsol(i)+Xshift
  Xsolshift=Xsolshift-Box(1)*FLOOR(Xsolshift/Box(1))  !Periodic BC in x direction
  
  Ysolshift=Ysol(i)+Yshift
  Ysolshift=Ysolshift-Box(2)*FLOOR(Ysolshift/Box(2))  !Periodic BC in y direction
  
  Zsolshift=Zsol(i)+Zshift
  Zsolshift=Zsolshift-Box(3)*FLOOR(Zsolshift/Box(3))  !Periodic BC in z direction

  !Find new cell N for solvent and find total solv. in new cell N
  CellN=INT(FLOOR(Xsolshift)+Box(1)*FLOOR(Ysolshift)+Box(1)*Box(2)*FLOOR(Zsolshift)+1)
  NewCellN(i)=CellN
  tCellshift(CellN)=tCellshift(CellN)+1

  !Sum the velocities in each new cell N; used next for center of mass velocity calc.
  VXCM(CellN)=VXCM(CellN)+VXsol(i)
  VYCM(CellN)=VYCM(CellN)+VYsol(i)
  VZCM(CellN)=VZCM(CellN)+VZsol(i)
 END IF
END DO

DO i=1,Vol
 IF (tCellshift(i) /= 0) THEN
  ItCellshift=1.d0/tCellshift(i)
  VXCM(i)=VXCM(i)*ItCellshift
  VYCM(i)=VYCM(i)*ItCellshift
  VZCM(i)=VZCM(i)*ItCellshift
 END IF
END DO

!!!!!!Find the new velocities for ith particle (from Kapral, 2008, "long MPC article")
!Calculate random unit vector in Cart.coord. from spherical polar coord. (r=1)

DO i=1,Vol
 CALL genrand_real1(x)
 CALL genrand_real1(y)
 phi=2*pi*x					!Random angle between 0 and 2*pi
 theta=2*y-1.d0
 UVX(i)=DSQRT(1.d0-theta*theta)*DCOS(phi)			!X component of random unit vector
 UVY(i)=DSQRT(1.d0-theta*theta)*DSIN(phi)			!Y component of random unit vector
 UVZ(i)=theta					!Z component of random unit vector
END DO

!Find the rotated comp. of the difference between CM velocity and particle i velocity

!$OMP PARALLEL DO PRIVATE (dVXsolVXCM,dVYsolVYCM,dVZsolVZCM,i) &
!$OMP PRIVATE (UVXortho,UVYortho,UVZortho,UVXortho2,UVYortho2) &
!$OMP PRIVATE (UVZortho2,UVXproj,UVYproj,UVZproj,ProjVonUV) &
!$OMP PRIVATE (dVXsolVXCMrot,dVYsolVYCMrot,dVZsolVZCMrot) &
!$OMP SHARED (VXCM,VYCM,VZCM,UVX,UVY,UVZ,FLAG,VXsol,VYsol,VZsol) &
!$OMP SCHEDULE(static)

DO i=1,NT,1
 IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
  !Difference between CM velocity and particle i velocity, unrotated
  dVXsolVXCM=VXsol(i)-VXCM(NewCellN(i))
  dVYsolVYCM=VYsol(i)-VYCM(NewCellN(i))
  dVZsolVZCM=VZsol(i)-VZCM(NewCellN(i))

  !Projection of difference onto unit vector
  ProjVonUV=dVXsolVXCM*UVX(NewCellN(i))+dVYsolVYCM*UVY(NewCellN(i)) & 
            +dVZsolVZCM*UVZ(NewCellN(i))
  UVXproj=ProjVonUV*UVX(NewCellN(i))
  UVYproj=ProjVonUV*UVY(NewCellN(i))
  UVZproj=ProjVonUV*UVZ(NewCellN(i))

  !Find vector orthogonal to unit vector, on same plane as the unrotated vel. vector
  UVXortho=dVXsolVXCM-UVXproj
  UVYortho=dVYsolVYCM-UVYproj
  UVZortho=dVZsolVZCM-UVZproj

  !Cross product of UV_ortho and unit vector; magnitude is same as UV_ortho
  UVXortho2=(UVYortho*UVZ(NewCellN(i)))-(UVZortho*UVY(NewCellN(i)))
  UVYortho2=(UVZortho*UVX(NewCellN(i)))-(UVXortho*UVZ(NewCellN(i)))
  UVZortho2=(UVXortho*UVY(NewCellN(i)))-(UVYortho*UVX(NewCellN(i)))

  !Rotated vector; additions of contributions from UV_proj, UV_ortho, and UV_ortho2
  dVXsolVXCMrot=UVXproj+UVXortho*cosrot+UVXortho2*sinrot
  dVYsolVYCMrot=UVYproj+UVYortho*cosrot+UVYortho2*sinrot
  dVZsolVZCMrot=UVZproj+UVZortho*cosrot+UVZortho2*sinrot

  !New velocities after collision
  VXsol(i)=VXCM(NewCellN(i))+dVXsolVXCMrot
  VYsol(i)=VYCM(NewCellN(i))+dVYsolVYCMrot
  VZsol(i)=VZCM(NewCellN(i))+dVZsolVZCMrot
 END IF
END DO
!$OMP END PARALLEL DO

!***********************************************************************************!
!******************************Reactions********************************************!
!***********************************************************************************!

DO i=1,Vol
IF (affectedCellN(i) == 0) THEN

 !!!!!!Probability calculations
 !Probability factor a_
 a1=k1
 a2=k2*(tFCellN(i))
 a5=k5*(tGCellN(i))
 a6=k6
 IF (tGCellN(i) == 0) THEN
  a3=0
 ELSE IF (tGCellN(i) /= 0) THEN
  a3=k3*(tFCellN(i))*(tGCellN(i))*(tGCellN(i)-1)
 END IF
 IF (tGCellN(i) <= 1) THEN
  a4=0
 ELSE IF (tGCellN(i) > 1) THEN
  a4=k4*(tGCellN(i))*(tGCellN(i)-1)*(tGCellN(i)-2)
 END IF
 IF (tICellN(i) == 0) THEN
  a1=0
  a6=0
 END IF

 a00=a1+a2+a3+a4+a5+a6

 !Probability for a reaction in the MPC time, and prob. for reaction of A,B,F,G
 ProbMPCt=1-EXP(-a00*MPCt)

 CALL genrand_real1(x)
 a0rand=a00*x

 CALL genrand_real1(x)
 IF (ProbMPCt > x) THEN
 IF (tICellN(i) > 0) THEN
  IF (tFCellN(i) == 0) THEN
   IF (tGCellN(i) == 0) THEN
    IF (a0rand < a1) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=0
    ELSE IF (a0rand >= a1 .AND. a0rand < a00) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=1
    END IF
   ELSE IF (tGCellN(i) > 0 .AND. tGCellN(i) <= 2) THEN
    IF (a0rand < a1) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=0     
    ELSE IF (a0rand >= a1 .AND. a0rand < a1+a5) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a5 .AND. a0rand < a00) THEN
     CellN=i 
     FLAG(sICellN(CellN,1))=1
    END IF
   ELSE IF (tGCellN(i) > 2) THEN
    IF (a0rand < a1) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=0     
    ELSE IF (a0rand >= a1 .AND. a0rand < a1+a4) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=0
    ELSE IF (a0rand >= a1+a4 .AND. a0rand < a1+a4+a5) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a4+a5 .AND. a0rand < a00) THEN
     CellN=i 
     FLAG(sICellN(CellN,1))=1
    END IF
   END IF
  ELSE IF (tFCellN(i) > 0) THEN
   IF (tGCellN(i) == 0) THEN
    IF (a0rand < a1) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=0     
    ELSE IF (a0rand >= a1 .AND. a0rand < a1+a2) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2 .AND. a0rand < a00) THEN
     CellN=i 
     FLAG(sICellN(CellN,1))=1
    END IF
   ELSE IF (tGCellN(i) == 1) THEN
    IF (a0rand < a1) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=0     
    ELSE IF (a0rand >= a1 .AND. a0rand < a1+a2) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2 .AND. a0rand < a1+a2+a5) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2+a5 .AND. a0rand < a00) THEN
     CellN=i 
     FLAG(sICellN(CellN,1))=1
    END IF
   ELSE IF (tGCellN(i) == 2) THEN
    IF (a0rand < a1) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=0     
    ELSE IF (a0rand >= a1 .AND. a0rand < a1+a2) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2 .AND. a0rand < a1+a2+a3) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=1
    ELSE IF (a0rand >= a1+a2+a3 .AND. a0rand < a1+a2+a3+a5) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2+a3+a5 .AND. a0rand < a00) THEN
     CellN=i 
     FLAG(sICellN(CellN,1))=1
    END IF
   ELSE IF (tGCellN(i) > 2) THEN
    IF (a0rand < a1) THEN
     CellN=i
     FLAG(sICellN(CellN,1))=0     
    ELSE IF (a0rand >= a1 .AND. a0rand < a1+a2) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2 .AND. a0rand < a1+a2+a3) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=1
    ELSE IF (a0rand >= a1+a2+a3 .AND. a0rand < a1+a2+a3+a4) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=0
    ELSE IF (a0rand >= a1+a2+a3+a4 .AND. a0rand < a1+a2+a3+a4+a5) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2+a3+a4+a5 .AND. a0rand < a00) THEN
     CellN=i 
     FLAG(sICellN(CellN,1))=1
    END IF
   END IF
  END IF
 ELSE IF (tICellN(i) == 0) THEN
  IF (tFCellN(i) == 0) THEN
   IF (tGCellN(i) == 0) THEN
   ! The case where no particles in the cell; SHOULD never happen!
   ELSE IF (tGCellN(i) > 0 .AND. tGCellN(i) <= 2) THEN
    IF (a0rand >= a1 .AND. a0rand < a00) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    END IF
   ELSE IF (tGCellN(i) > 2) THEN
    IF (a0rand >= a1 .AND. a0rand < a1+a4) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=0
    ELSE IF (a0rand >= a1+a4 .AND. a0rand < a00) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    END IF
   END IF
  ELSE IF (tFCellN(i) > 0) THEN
   IF (tGCellN(i) == 0) THEN
    IF (a0rand >= a1 .AND. a0rand < a00) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    END IF
   ELSE IF (tGCellN(i) == 1) THEN
    IF (a0rand >= a1 .AND. a0rand < a1+a2) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2 .AND. a0rand < a00) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    END IF
   ELSE IF (tGCellN(i) == 2) THEN
    IF (a0rand >= a1 .AND. a0rand < a1+a2) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2 .AND. a0rand < a1+a2+a3) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=1
    ELSE IF (a0rand >= a1+a2+a3 .AND. a0rand < a00) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    END IF
   ELSE IF (tGCellN(i) > 2) THEN    
    IF (a0rand >= a1 .AND. a0rand < a1+a2) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=3
    ELSE IF (a0rand >= a1+a2 .AND. a0rand < a1+a2+a3) THEN
     CellN=i
     FLAG(sFCellN(CellN,1))=1
    ELSE IF (a0rand >= a1+a2+a3 .AND. a0rand < a1+a2+a3+a4) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=0
    ELSE IF (a0rand >= a1+a2+a3+a4 .AND. a0rand < a00) THEN
     CellN=i
     FLAG(sGCellN(CellN,1))=3
    END IF
   END IF
  END IF
 END IF
 END IF
END IF
END DO

END SUBROUTINE RMPC_SelkovDimer
