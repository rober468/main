SUBROUTINE Restart
USE Global_SelkovDimer
IMPLICIT NONE

pstep=cstep
cstep=cstep+1
rcstep=(pstep/Freq1)+1

!!!!!!Record total number of F and G molecules
NT=iTotF+iTotG+iTotI

 IF (cstep > 0 .AND. cstep < 100) THEN
  WRITE(SolventFPosFile,"(A11,I2,A4)") 'SolventFPos',pstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I2,A4)") 'SolventGPos',pstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I2,A4)") 'SolventIPos',pstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I2,A4)") 'SolventFVel',pstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I2,A4)") 'SolventGVel',pstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I2,A4)") 'SolventIVel',pstep,'.dat'
  WRITE(SolventFForFile,"(A11,I2,A4)") 'SolventFFor',pstep,'.dat'
  WRITE(SolventGForFile,"(A11,I2,A4)") 'SolventGFor',pstep,'.dat'
  WRITE(SolventIForFile,"(A11,I2,A4)") 'SolventIFor',pstep,'.dat'
 ELSE IF (cstep >= 100 .AND. cstep < 1000) THEN
  WRITE(SolventFPosFile,"(A11,I3,A4)") 'SolventFPos',pstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I3,A4)") 'SolventGPos',pstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I3,A4)") 'SolventIPos',pstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I3,A4)") 'SolventFVel',pstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I3,A4)") 'SolventGVel',pstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I3,A4)") 'SolventIVel',pstep,'.dat'
  WRITE(SolventFForFile,"(A11,I3,A4)") 'SolventFFor',pstep,'.dat'
  WRITE(SolventGForFile,"(A11,I3,A4)") 'SolventGFor',pstep,'.dat'
  WRITE(SolventIForFile,"(A11,I3,A4)") 'SolventIFor',pstep,'.dat'
 ELSE IF (cstep >= 1000 .AND. cstep < 10000) THEN
  WRITE(SolventFPosFile,"(A11,I4,A4)") 'SolventFPos',pstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I4,A4)") 'SolventGPos',pstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I4,A4)") 'SolventIPos',pstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I4,A4)") 'SolventFVel',pstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I4,A4)") 'SolventGVel',pstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I4,A4)") 'SolventIVel',pstep,'.dat'
  WRITE(SolventFForFile,"(A11,I4,A4)") 'SolventFFor',pstep,'.dat'
  WRITE(SolventGForFile,"(A11,I4,A4)") 'SolventGFor',pstep,'.dat'
  WRITE(SolventIForFile,"(A11,I4,A4)") 'SolventIFor',pstep,'.dat'
 ELSE IF (cstep >= 10000 .AND. cstep < 100000) THEN
  WRITE(SolventFPosFile,"(A11,I5,A4)") 'SolventFPos',pstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I5,A4)") 'SolventGPos',pstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I5,A4)") 'SolventIPos',pstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I5,A4)") 'SolventFVel',pstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I5,A4)") 'SolventGVel',pstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I5,A4)") 'SolventIVel',pstep,'.dat'
  WRITE(SolventFForFile,"(A11,I5,A4)") 'SolventFFor',pstep,'.dat'
  WRITE(SolventGForFile,"(A11,I5,A4)") 'SolventGFor',pstep,'.dat'
  WRITE(SolventIForFile,"(A11,I5,A4)") 'SolventIFor',pstep,'.dat'
 ELSE IF (cstep >= 100000 .AND. cstep < 1000000) THEN
  WRITE(SolventFPosFile,"(A11,I6,A4)") 'SolventFPos',pstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I6,A4)") 'SolventGPos',pstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I6,A4)") 'SolventIPos',pstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I6,A4)") 'SolventFVel',pstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I6,A4)") 'SolventGVel',pstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I6,A4)") 'SolventIVel',pstep,'.dat'
  WRITE(SolventFForFile,"(A11,I6,A4)") 'SolventFFor',pstep,'.dat'
  WRITE(SolventGForFile,"(A11,I6,A4)") 'SolventGFor',pstep,'.dat'
  WRITE(SolventIForFile,"(A11,I6,A4)") 'SolventIFor',pstep,'.dat'
 ELSE IF (cstep >= 1000000 .AND. cstep < 10000000) THEN
  WRITE(SolventFPosFile,"(A11,I7,A4)") 'SolventFPos',pstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I7,A4)") 'SolventGPos',pstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I7,A4)") 'SolventIPos',pstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I7,A4)") 'SolventFVel',pstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I7,A4)") 'SolventGVel',pstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I7,A4)") 'SolventIVel',pstep,'.dat'
  WRITE(SolventFForFile,"(A11,I7,A4)") 'SolventFFor',pstep,'.dat'
  WRITE(SolventGForFile,"(A11,I7,A4)") 'SolventGFor',pstep,'.dat'
  WRITE(SolventIForFile,"(A11,I7,A4)") 'SolventIFor',pstep,'.dat'
 ELSE IF (cstep >= 10000000 .AND. cstep < 100000000) THEN
  WRITE(SolventFPosFile,"(A11,I8,A4)") 'SolventFPos',pstep,'.dat'
  WRITE(SolventGPosFile,"(A11,I8,A4)") 'SolventGPos',pstep,'.dat'
  WRITE(SolventIPosFile,"(A11,I8,A4)") 'SolventIPos',pstep,'.dat'
  WRITE(SolventFVelFile,"(A11,I8,A4)") 'SolventFVel',pstep,'.dat'
  WRITE(SolventGVelFile,"(A11,I8,A4)") 'SolventGVel',pstep,'.dat'
  WRITE(SolventIVelFile,"(A11,I8,A4)") 'SolventIVel',pstep,'.dat'
  WRITE(SolventFForFile,"(A11,I8,A4)") 'SolventFFor',pstep,'.dat'
  WRITE(SolventGForFile,"(A11,I8,A4)") 'SolventGFor',pstep,'.dat'
  WRITE(SolventIForFile,"(A11,I8,A4)") 'SolventIFor',pstep,'.dat'
 END IF

!!!!!!Record positions and velocity of F solvent
OPEN(UNIT=10, FILE=SolventFPosFile, STATUS='OLD')
OPEN(UNIT=20, FILE=SolventFVelFile, STATUS='OLD')
OPEN(UNIT=30, FILE=SolventFForFile, STATUS='OLD')
DO j=1,iTotF
 READ(10,21)ab,ac,ad
 READ(20,21)VXsol(j),VYsol(j),VZsol(j)
 READ(30,21)FXsol(j),FYsol(j),FZsol(j)
 Xsol(j)=ab
 Ysol(j)=ac
 Zsol(j)=ad
 FLAG(j)=0
END DO
CLOSE(10)
CLOSE(20)
CLOSE(30)

21 FORMAT(3(ES23.16,1X),I7)

!!!!!!Record positions and velocity of G solvent
OPEN(UNIT=10, FILE=SolventGPosFile, STATUS='OLD')
OPEN(UNIT=20, FILE=SolventGVelFile, STATUS='OLD')
OPEN(UNIT=30, FILE=SolventGForFile, STATUS='OLD')
DO j=iTotF+1,iTotF+iTotG
 READ(10,21)ab,ac,ad
 READ(20,21)VXsol(j),VYsol(j),VZsol(j)
 READ(30,21)FXsol(j),FYsol(j),FZsol(j)
 Xsol(j)=ab
 Ysol(j)=ac
 Zsol(j)=ad
 FLAG(j)=1
END DO
CLOSE(10)
CLOSE(20)
CLOSE(30)

!!!!!!Record positions and velocity of I solvent
OPEN(UNIT=10, FILE=SolventIPosFile, STATUS='OLD')
OPEN(UNIT=20, FILE=SolventIVelFile, STATUS='OLD')
OPEN(UNIT=30, FILE=SolventIForFile, STATUS='OLD')
DO j=iTotF+iTotG+1,NT
 READ(10,21)ab,ac,ad
 READ(20,21)VXsol(j),VYsol(j),VZsol(j)
 READ(30,21)FXsol(j),FYsol(j),FZsol(j)
 Xsol(j)=ab
 Ysol(j)=ac
 Zsol(j)=ad
 FLAG(j)=3
END DO
CLOSE(10)
CLOSE(20)
CLOSE(30)

!!!!!!Record positions and velocity of dimer
OPEN(UNIT=10, FILE='DimerPosC.dat', STATUS='OLD')
OPEN(UNIT=20, FILE='DimerPosN.dat', STATUS='OLD')
OPEN(UNIT=30, FILE='DimerVelC.dat', STATUS='OLD')
OPEN(UNIT=40, FILE='DimerVelN.dat', STATUS='OLD')
OPEN(UNIT=31, FILE='DimerForC.dat', STATUS='OLD')
OPEN(UNIT=32, FILE='DimerForN.dat', STATUS='OLD')
DO i=1,rcstep
 IF (i < rcstep) THEN
  READ(10,22)
  READ(20,22)
  READ(30,22)
  READ(40,22)
  READ(31,22)
  READ(32,22)
 ELSE IF (i == rcstep) THEN
  READ(10,22)j,Xdim(1),Ydim(1),Zdim(1)
  READ(20,22)j,Xdim(2),Ydim(2),Zdim(2)
  READ(30,22)j,VXdim(1),VYdim(1),VZdim(1)
  READ(40,22)j,VXdim(2),VYdim(2),VZdim(2)
  READ(31,22)j,FXdim(1),FYdim(1),FZdim(1)
  READ(32,22)j,FXdim(2),FYdim(2),FZdim(2)
 END IF
END DO
CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
CLOSE(31)
CLOSE(32)

22 FORMAT(I8,3(1X,E23.16))

END SUBROUTINE Restart
