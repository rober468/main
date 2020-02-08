SUBROUTINE Initial_SelkovDimer
USE Global_SelkovDimer
USE mt95
IMPLICIT NONE

!!!!!!Initial position of dimer
Xdim(1)=HBox(1)+HDcn		!Catalytic sphere init. position, x
Ydim(1)=HBox(2)			!Catalytic sphere init. position, y
Zdim(1)=HBox(3)			!Catalytic sphere init. position, z
Xdim(2)=HBox(1)-HDcn		!Non-catalytic sphere init. position, x
Ydim(2)=HBox(2)			!Non-catalytic sphere init. position, y
Zdim(2)=HBox(3)			!Non-catalytic sphere init. position, z

!!!!!!Initial position of solvent F
DO i=1,Nf0
 FLAG(i)=0			!Solvent particle is initially F

 100 CALL genrand_real1(x)
     CALL genrand_real1(y)
     CALL genrand_real1(z)

     Xsol(i)=Box(1)*x		!Random position of solvent in box, x
     Ysol(i)=Box(2)*y		!Random position of solvent in box, y
     Zsol(i)=Box(3)*z		!Random position of solvent in box, z

 !Make sure F solvent isn't overlapping the C dimer sphere; if alright, go to N
 R(1,1)=Xsol(i)
 R(1,2)=Ysol(i)
 R(1,3)=Zsol(i)
 R(2,1)=Xdim(1)
 R(2,2)=Ydim(1)
 R(2,3)=Zdim(1)

 DO j=1,3
  dR(j)=R(1,j)-R(2,j)
  IF (dR(j) > HBox(j)) THEN		
   dR(j)=dR(j)-Box(j)
  ELSE IF (dR(j) < -HBox(j)) THEN
   dR(j)=dR(j)+Box(j)
  END IF				
 END DO
 dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3)

 IF (dRsq < CUTOFFc2) THEN 		!If outside of C sphere cutoff distance, then 
  GOTO 100				!move on to check with N; if not, restart
 ELSE

  !Make sure F solvent isn't overlapping N dimer sphere; if alright,then position good
  R(1,1)=Xsol(i)
  R(1,2)=Ysol(i)
  R(1,3)=Zsol(i)
  R(2,1)=Xdim(2)
  R(2,2)=Ydim(2)
  R(2,3)=Zdim(2)

  DO j=1,3
   dR(j)=R(1,j)-R(2,j)
   IF (dR(j) > HBox(j)) THEN		
    dR(j)=dR(j)-Box(j)
   ELSE IF (dR(j) < -HBox(j)) THEN
    dR(j)=dR(j)+Box(j)
   END IF				
  END DO
  dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3)

  IF (dRsq < CUTOFFn2) THEN 		!If outside of C sphere cutoff distance, then 
   GOTO 100				!move on to check with N; if not, restart
  END IF
 END IF
END DO

!!!!!!Initial position of solvent G
DO i=Nf0+1,Nf0+Ng0
 FLAG(i)=1			!Solvent particle is initially G

 200 CALL genrand_real1(x)
     CALL genrand_real1(y)
     CALL genrand_real1(z)

     Xsol(i)=Box(1)*x		!Random position of solvent in box, x
     Ysol(i)=Box(2)*y		!Random position of solvent in box, y
     Zsol(i)=Box(3)*z		!Random position of solvent in box, z

 !Make sure G solvent isn't overlapping the C dimer sphere; if alright, go to N
 R(1,1)=Xsol(i)
 R(1,2)=Ysol(i)
 R(1,3)=Zsol(i)
 R(2,1)=Xdim(1)
 R(2,2)=Ydim(1)
 R(2,3)=Zdim(1)

 DO j=1,3
  dR(j)=R(1,j)-R(2,j)
  IF (dR(j) > HBox(j)) THEN		
   dR(j)=dR(j)-Box(j)
  ELSE IF (dR(j) < -HBox(j)) THEN
   dR(j)=dR(j)+Box(j)
  END IF				
 END DO
 dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3)

 IF (dRsq < CUTOFFc2) THEN 		!If outside of C sphere cutoff distance, then 
  GOTO 200				!move on to check with N; if not, restart
 ELSE

  !Make sure G solvent isn't overlapping N dimer sphere; if alright,then position good
  R(1,1)=Xsol(i)
  R(1,2)=Ysol(i)
  R(1,3)=Zsol(i)
  R(2,1)=Xdim(2)
  R(2,2)=Ydim(2)
  R(2,3)=Zdim(2)

  DO j=1,3
   dR(j)=R(1,j)-R(2,j)
   IF (dR(j) > HBox(j)) THEN		
    dR(j)=dR(j)-Box(j)
   ELSE IF (dR(j) < -HBox(j)) THEN
    dR(j)=dR(j)+Box(j)
   END IF				
  END DO
  dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3)

  IF (dRsq < CUTOFFn2) THEN 		!If outside of C sphere cutoff distance, then 
   GOTO 200				!move on to check with N; if not, restart
  END IF
 END IF
END DO

!!!!!!Initial position of inert solvent
DO i=Nf0+Ng0+1,NT
 FLAG(i)=3			!Solvent particle is initially inert

 300 CALL genrand_real1(x)
     CALL genrand_real1(y)
     CALL genrand_real1(z)

     Xsol(i)=Box(1)*x		!Random position of solvent in box, x
     Ysol(i)=Box(2)*y		!Random position of solvent in box, y
     Zsol(i)=Box(3)*z		!Random position of solvent in box, z

 !Make sure inert solvent isn't overlapping the C dimer sphere; if alright, go to N
 R(1,1)=Xsol(i)
 R(1,2)=Ysol(i)
 R(1,3)=Zsol(i)
 R(2,1)=Xdim(1)
 R(2,2)=Ydim(1)
 R(2,3)=Zdim(1)

 DO j=1,3
  dR(j)=R(1,j)-R(2,j)
  IF (dR(j) > HBox(j)) THEN		
   dR(j)=dR(j)-Box(j)
  ELSE IF (dR(j) < -HBox(j)) THEN
   dR(j)=dR(j)+Box(j)
  END IF				
 END DO
 dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3)

 IF (dRsq < CUTOFFc2) THEN 		!If outside of C sphere cutoff distance, then 
  GOTO 300				!move on to check with N; if not, restart
 ELSE

  !Make sure inert solvent isn't overlapping N dimer sphere; if alright,then position good
  R(1,1)=Xsol(i)
  R(1,2)=Ysol(i)
  R(1,3)=Zsol(i)
  R(2,1)=Xdim(2)
  R(2,2)=Ydim(2)
  R(2,3)=Zdim(2)

  DO j=1,3
   dR(j)=R(1,j)-R(2,j)
   IF (dR(j) > HBox(j)) THEN		
    dR(j)=dR(j)-Box(j)
   ELSE IF (dR(j) < -HBox(j)) THEN
    dR(j)=dR(j)+Box(j)
   END IF				
  END DO
  dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3)

  IF (dRsq < CUTOFFn2) THEN 		!If outside of C sphere cutoff distance, then 
   GOTO 300				!move on to check with N; if not, restart
  END IF
 END IF
END DO

!!!!!!Initial velocity of dimer
DO i=1,2
 VXdim(i)=0.
 VYdim(i)=0.
 VZdim(i)=0.
END DO

sumVXsol=0.
sumVYsol=0.
sumVZsol=0.

!!!!!!Initial velocity of solvent F
!For x direction
DO i=1,Nf0,1					!Rand. F velocity from Gauss. distr.,x
 CALL Gaussian_Random (mean,varF,VXsol(i))
 sumVXsol=sumVXsol+VXsol(i)
END DO
!For y direction
DO i=1,Nf0,1					!Rand. F velocity from Gauss. distr.,y
 CALL Gaussian_Random (mean,varF,VYsol(i))
 sumVYsol=sumVYsol+VYsol(i)
END DO
!For z direction
DO i=1,Nf0,1					!Rand. F velocity from Gauss. distr.,z
 CALL Gaussian_Random (mean,varF,VZsol(i))
 sumVZsol=sumVZsol+VZsol(i)
END DO

!!!!!!Initial velocity of solvent G
!For x direction
DO i=Nf0+1,Nf0+Ng0				!Rand. G velocity from Gauss. distr.,x
 CALL Gaussian_Random (mean,varG,VXsol(i))
 sumVXsol=sumVXsol+VXsol(i)
END DO
!For y direction
DO i=Nf0+1,Nf0+Ng0				!Rand. G velocity from Gauss. distr.,y
 CALL Gaussian_Random (mean,varG,VYsol(i))
 sumVYsol=sumVYsol+VYsol(i)
END DO
!For z direction
DO i=Nf0+1,Nf0+Ng0				!Rand. G velocity from Gauss. distr.,z
 CALL Gaussian_Random (mean,varG,VZsol(i))
 sumVZsol=sumVZsol+VZsol(i)
END DO

!!!!!!Initial velocity of inert solvent
!For x direction
DO i=Nf0+Ng0+1,NT				!Rand. G velocity from Gauss. distr.,x
 CALL Gaussian_Random (mean,varI,VXsol(i))
 sumVXsol=sumVXsol+VXsol(i)
END DO
!For y direction
DO i=Nf0+Ng0+1,NT				!Rand. G velocity from Gauss. distr.,y
 CALL Gaussian_Random (mean,varI,VYsol(i))
 sumVYsol=sumVYsol+VYsol(i)
END DO
!For z direction
DO i=Nf0+Ng0+1,NT				!Rand. G velocity from Gauss. distr.,z
 CALL Gaussian_Random (mean,varI,VZsol(i))
 sumVZsol=sumVZsol+VZsol(i)
END DO

!!!!!!Make total momentum zero
!Average velocity in each direction
sumVXsol=sumVXsol/DFLOAT(NT)
sumVYsol=sumVYsol/DFLOAT(NT)
sumVZsol=sumVZsol/DFLOAT(NT)

!Shift and adjust
DO i=1,NT
 VXsol(i)=VXsol(i)-sumVXsol
 VYsol(i)=VYsol(i)-sumVYsol
 VZsol(i)=VZsol(i)-sumVZsol
END DO

!!!!!!Initial force on solvent and dimer
DO i=1,NT,1					!Initial force on solvent is zero
 FXsol(i)=0.
 FYsol(i)=0.
 FZsol(i)=0.
END DO
DO i=1,2,1					!Initial force on dimer is zero
 FXdim(i)=0.
 FYdim(i)=0.
 FZdim(i)=0.
END DO

!!!!!!Writing the initial information
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

OPEN(UNIT=10,FILE='SolventFPos0.dat',status='replace',ACTION='WRITE')
OPEN(UNIT=20,FILE='SolventGPos0.dat',status='replace',ACTION='WRITE')
OPEN(UNIT=30,FILE='SolventIPos0.dat',status='replace',ACTION='WRITE')
OPEN(UNIT=40,FILE='SolventFVel0.dat',status='replace',ACTION='WRITE')
OPEN(UNIT=31,FILE='SolventGVel0.dat',status='replace',ACTION='WRITE')
OPEN(UNIT=32,FILE='SolventIVel0.dat',status='replace',ACTION='WRITE')
DO i=1,NT
 IF (FLAG(i) == 0) THEN
  WRITE(10,14) Xsol(i),Ysol(i),Zsol(i), i
  WRITE(40,14) VXsol(i),VYsol(i),VZsol(i), i
 ELSE IF (FLAG(i) == 1) THEN
  WRITE(20,14) Xsol(i),Ysol(i),Zsol(i), i
  WRITE(31,14) VXsol(i),VYsol(i),VZsol(i), i
 ELSE IF (FLAG(i) == 3) THEN
  WRITE(30,14) Xsol(i),Ysol(i),Zsol(i), i
  WRITE(32,14) VXsol(i),VYsol(i),VZsol(i), i
 END IF
END DO

14 FORMAT(3(ES23.16,1X),I8)

CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
CLOSE(31)
CLOSE(32)

END SUBROUTINE Initial_SelkovDimer
