SUBROUTINE Distribute_SelkovDimer
USE Global_SelkovDimer
IMPLICIT NONE

!!!!!!Distribute the molecules to the cells
!Clear previous values and start with all values zero

!$OMP PARALLEL DO &
!$OMP PRIVATE(i) &
!$OMP SHARED(tFCellN,tGCellN,sFCellN,sGCellN)

DO i=1,Vol
 tFCellN(i)=0
 tGCellN(i)=0
 tICellN(i)=0
END DO

!$OMP END PARALLEL DO

!Assign molecules to cells and designate with appropriate labels
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
  tICellN(CellN)=tICellN(CellN)+1       !Total number of I in cell N
  sICellN(CellN,tICellN(CellN))=i       !Specify/label newly added I with i
 END IF
END DO

END SUBROUTINE Distribute_SelkovDimer
