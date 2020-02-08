SUBROUTINE Sort
USE Global_SelkovDimer
IMPLICIT NONE

!This is a subroutine within "Distribute" that sorts the memory locations of the solvent particles into a new order that is reflective of their current positions relative to each other in the physical system. The sorting is done taking the current solvent particles (and their indices), assigning them to a specific box, and numbering them within that box.  This is already done with "Distribute". They are then saved into a new matrix with new indices according to the box order; this ensures particles in the same box are saved in closed proximity in memory.

new=0

!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SHARED(tFCellN,tGCellN,sFCellN,sGCellN,tICellN,sICellN)

DO i=1,Vol
 tFCellN(i)=0
 tGCellN(i)=0
 tICellN(i)=0
 DO j=1,MaxF
  sFCellN(i,j)=0
 END DO
 DO j=1,MaxG
  sGCellN(i,j)=0
 END DO
 DO j=1,MaxI
  sICellN(i,j)=0
 END DO
END DO

!$OMP END PARALLEL DO

DO i=1,NT 
 IF (FLAG(i) == 0) THEN
  CellN=INT(FLOOR(Xsol(i))+Box(1)*FLOOR(Ysol(i))+Box(1)*Box(2)*FLOOR(Zsol(i))+1) !Identify cell
  tFCellN(CellN)=tFCellN(CellN)+1	!Total number of F in cell N
  sFCellN(CellN,tFCellN(CellN))=i	!Specify/label newly added F with i
 ELSE IF (FLAG(i) == 1) THEN
  CellN=INT(FLOOR(Xsol(i))+Box(1)*FLOOR(Ysol(i))+Box(1)*Box(2)*FLOOR(Zsol(i))+1) !Identify cell
  tGCellN(CellN)=tGCellN(CellN)+1	!Total number of G in cell N
  sGCellN(CellN,tGCellN(CellN))=i	!Specify/label newly added G with i
 ELSE IF (FLAG(i) == 3) THEN
  CellN=INT(FLOOR(Xsol(i))+Box(1)*FLOOR(Ysol(i))+Box(1)*Box(2)*FLOOR(Zsol(i))+1) !Identify cell
  tICellN(CellN)=tICellN(CellN)+1	!Total number of I in cell N
  sICellN(CellN,tICellN(CellN))=i	!Specify/label newly added I with i
 END IF
END DO

DO m=1,Vol			!Loop over all cells
 IF (tFCellN(m) == 0) THEN
 ELSE IF (tFCellN(m) == 1) THEN
  new=new+1
  newXsol(new)=Xsol(sFCellN(m,1))
  newYsol(new)=Ysol(sFCellN(m,1))
  newZsol(new)=Zsol(sFCellN(m,1))
  newFLAG(new)=FLAG(sFCellN(m,1))
  newVXsol(new)=VXsol(sFCellN(m,1))
  newVYsol(new)=VYsol(sFCellN(m,1))
  newVZsol(new)=VZsol(sFCellN(m,1))
  newFXsol(new)=FXsol(sFCellN(m,1))
  newFYsol(new)=FYsol(sFCellN(m,1))
  newFZsol(new)=FZsol(sFCellN(m,1))
  newpFXsol(new)=pFXsol(sFCellN(m,1))
  newpFYsol(new)=pFYsol(sFCellN(m,1))
  newpFZsol(new)=pFZsol(sFCellN(m,1))
 ELSE IF (tFCellN(m) > 1) THEN
  DO p=1,tFCellN(m)		!Loop over each F solvent in mth cell
   new=new+1
   newXsol(new)=Xsol(sFCellN(m,p))
   newYsol(new)=Ysol(sFCellN(m,p))
   newZsol(new)=Zsol(sFCellN(m,p))
   newFLAG(new)=FLAG(sFCellN(m,p))
   newVXsol(new)=VXsol(sFCellN(m,p))
   newVYsol(new)=VYsol(sFCellN(m,p))
   newVZsol(new)=VZsol(sFCellN(m,p))
   newFXsol(new)=FXsol(sFCellN(m,p))
   newFYsol(new)=FYsol(sFCellN(m,p))
   newFZsol(new)=FZsol(sFCellN(m,p))
   newpFXsol(new)=pFXsol(sFCellN(m,p))
   newpFYsol(new)=pFYsol(sFCellN(m,p))
   newpFZsol(new)=pFZsol(sFCellN(m,p))
  END DO
 END IF
 IF (tGCellN(m) == 0) THEN
 ELSE IF (tGCellN(m) == 1) THEN
  new=new+1
  newXsol(new)=Xsol(sGCellN(m,1))
  newYsol(new)=Ysol(sGCellN(m,1))
  newZsol(new)=Zsol(sGCellN(m,1))
  newFLAG(new)=FLAG(sGCellN(m,1))
  newVXsol(new)=VXsol(sGCellN(m,1))
  newVYsol(new)=VYsol(sGCellN(m,1))
  newVZsol(new)=VZsol(sGCellN(m,1))
  newFXsol(new)=FXsol(sGCellN(m,1))
  newFYsol(new)=FYsol(sGCellN(m,1))
  newFZsol(new)=FZsol(sGCellN(m,1))
  newpFXsol(new)=pFXsol(sGCellN(m,1))
  newpFYsol(new)=pFYsol(sGCellN(m,1))
  newpFZsol(new)=pFZsol(sGCellN(m,1))
 ELSE IF (tGCellN(m) > 1) THEN
  DO p=1,tGCellN(m)		!Loop over each G solvent in mth cell
   new=new+1
   newXsol(new)=Xsol(sGCellN(m,p))
   newYsol(new)=Ysol(sGCellN(m,p))
   newZsol(new)=Zsol(sGCellN(m,p))
   newFLAG(new)=FLAG(sGCellN(m,p))
   newVXsol(new)=VXsol(sGCellN(m,p))
   newVYsol(new)=VYsol(sGCellN(m,p))
   newVZsol(new)=VZsol(sGCellN(m,p))
   newFXsol(new)=FXsol(sGCellN(m,p))
   newFYsol(new)=FYsol(sGCellN(m,p))
   newFZsol(new)=FZsol(sGCellN(m,p))
   newpFXsol(new)=pFXsol(sGCellN(m,p))
   newpFYsol(new)=pFYsol(sGCellN(m,p))
   newpFZsol(new)=pFZsol(sGCellN(m,p))
  END DO
 END IF
 IF (tICellN(m) == 0) THEN
 ELSE IF (tICellN(m) == 1) THEN
  new=new+1
  newXsol(new)=Xsol(sICellN(m,1))
  newYsol(new)=Ysol(sICellN(m,1))
  newZsol(new)=Zsol(sICellN(m,1))
  newFLAG(new)=FLAG(sICellN(m,1))
  newVXsol(new)=VXsol(sICellN(m,1))
  newVYsol(new)=VYsol(sICellN(m,1))
  newVZsol(new)=VZsol(sICellN(m,1))
  newFXsol(new)=FXsol(sICellN(m,1))
  newFYsol(new)=FYsol(sICellN(m,1))
  newFZsol(new)=FZsol(sICellN(m,1))
  newpFXsol(new)=pFXsol(sICellN(m,1))
  newpFYsol(new)=pFYsol(sICellN(m,1))
  newpFZsol(new)=pFZsol(sICellN(m,1))
 ELSE IF (tICellN(m) > 1) THEN
  DO p=1,tICellN(m)		!Loop over each I solvent in mth cell
   new=new+1
   newXsol(new)=Xsol(sICellN(m,p))
   newYsol(new)=Ysol(sICellN(m,p))
   newZsol(new)=Zsol(sICellN(m,p))
   newFLAG(new)=FLAG(sICellN(m,p))
   newVXsol(new)=VXsol(sICellN(m,p))
   newVYsol(new)=VYsol(sICellN(m,p))
   newVZsol(new)=VZsol(sICellN(m,p))
   newFXsol(new)=FXsol(sICellN(m,p))
   newFYsol(new)=FYsol(sICellN(m,p))
   newFZsol(new)=FZsol(sICellN(m,p))
   newpFXsol(new)=pFXsol(sICellN(m,p))
   newpFYsol(new)=pFYsol(sICellN(m,p))
   newpFZsol(new)=pFZsol(sICellN(m,p))
  END DO
 END IF
END DO

DO m=1,NT
 IF (m <= new) THEN
  Xsol(m)=newXsol(m)
  Ysol(m)=newYsol(m)
  Zsol(m)=newZsol(m)
  FLAG(m)=newFLAG(m)
  VXsol(m)=newVXsol(m)
  VYsol(m)=newVYsol(m)
  VZsol(m)=newVZsol(m)
  FXsol(m)=newFXsol(m)
  FYsol(m)=newFYsol(m)
  FZsol(m)=newFZsol(m)
  pFXsol(m)=newpFXsol(m)
  pFYsol(m)=newpFYsol(m)
  pFZsol(m)=newpFZsol(m)
 ELSE IF (m > new) THEN
  Xsol(m)=0.
  Ysol(m)=0.
  Zsol(m)=0.
  FLAG(m)=2
  VXsol(m)=0.
  VYsol(m)=0.
  VZsol(m)=0.
  FXsol(m)=0.
  FYsol(m)=0.
  FZsol(m)=0.
  pFXsol(m)=0.
  pFYsol(m)=0.
  pFZsol(m)=0.
 END IF
END DO

NT=new

END SUBROUTINE Sort
