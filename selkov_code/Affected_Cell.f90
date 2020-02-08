SUBROUTINE Affected_Cell
USE GLobal_SelkovDimer
IMPLICIT NONE

DO i=1,Vol
 affectedCellN(i)=0
END DO

lowCXdim=CEILING(Xdim(1)-(CUTOFFc+1.d0))!Lower limit of affected cells in x dir.
upCXdim=FLOOR(Xdim(1)+(CUTOFFc+1.d0))   !Upper limit of affected cells in x dir.
lowCYdim=CEILING(Ydim(1)-(CUTOFFc+1.d0))!Lower limit of affected cells in y dir.
upCYdim=FLOOR(Ydim(1)+(CUTOFFc+1.d0))   !Upper limit of affected cells in y dir.
lowCZdim=CEILING(Zdim(1)-(CUTOFFc+1.d0))!Lower limit of affected cells in z dir.
upCZdim=FLOOR(Zdim(1)+(CUTOFFc+1.d0))   !Upper limit of affected cells in z dir.

!$OMP PARALLEL DO &
!$OMP PRIVATE(dimXCell,dimYCell,dimZCell,CellN,i,j,k) &
!$OMP SHARED(iBox)
DO i=lowCXdim,upCXdim
 DO j=lowCYdim,upCYdim
  DO k=lowCZdim,upCZdim
   IF ((ABS(REAL(i)-Xdim(1))-1.)**2 + (ABS(REAL(j)-Ydim(1))-1.)**2 &
        + (ABS(REAL(k)-Zdim(1))-1.)**2 < CUTOFFc2) THEN			!If part of box is in cutoff region of C

    dimXCell=i				!Find ith affected cell,x direction, PBC
    IF (dimXCell < 0) THEN
     dimXCell=dimXCell+iBox(1)
    ELSE IF (dimXCell >= iBox(1)) THEN
     dimXCell=dimXCell-iBox(1)
    END IF
    dimYCell=j				!Find jth affected cell,y direction, PBC
    IF (dimYCell < 0) THEN
     dimYCell=dimYCell+iBox(2)
    ELSE IF (dimYCell >= iBox(2)) THEN
     dimYCell=dimYCell-iBox(2)
    END IF
    dimZCell=k				!Find kth affected cell,z direction, PBC
    IF (dimZCell < 0) THEN
     dimZCell=dimZCell+iBox(3)
    ELSE IF (dimZCell >= iBox(3)) THEN
     dimZCell=dimZCell-iBox(3)
    END IF

    CellN=dimXCell+iBox(1)*dimYCell+iBox(1)*iBox(2)*dimZCell+1
    !CellN denotes the affected cell given by i=dimXCell, j=dimYCell, k=dimZCell
    affectedCellN(CellN)=1

   END IF
  END DO
 END DO
END DO
!$OMP END PARALLEL DO


lowNXdim=CEILING(Xdim(2)-(CUTOFFn+1.d0))!Lower limit of affected cells in x dir.
upNXdim=FLOOR(Xdim(2)+(CUTOFFn+1.d0))   !Upper limit of affected cells in x dir.
lowNYdim=CEILING(Ydim(2)-(CUTOFFn+1.d0))!Lower limit of affected cells in y dir.
upNYdim=FLOOR(Ydim(2)+(CUTOFFn+1.d0))   !Upper limit of affected cells in y dir.
lowNZdim=CEILING(Zdim(2)-(CUTOFFn+1.d0))!Lower limit of affected cells in z dir.
upNZdim=FLOOR(Zdim(2)+(CUTOFFn+1.d0))   !Upper limit of affected cells in z dir.

!$OMP PARALLEL DO &
!$OMP PRIVATE(dimXCell,dimYCell,dimZCell,CellN,i,j,k) &
!$OMP SHARED(iBox)

DO i=lowNXdim,upNXdim,1
 DO j=lowNYdim,upNYdim,1
  DO k=lowNZdim,upNZdim,1
   IF ((ABS(REAL(i)-Xdim(2))-1.)**2 + (ABS(REAL(j)-Ydim(2))-1.)**2 &
        + (ABS(REAL(k)-Zdim(2))-1.)**2 < CUTOFFn2) THEN			!If part of box is in cutoff region of N

    dimXCell=i				!Find ith affected cell,x direction, PBC
    IF (dimXCell < 0) THEN
     dimXCell=dimXCell+iBox(1)
    ELSE IF (dimXCell >= iBox(1)) THEN
     dimXCell=dimXCell-iBox(1)
    END IF
    dimYCell=j				!Find jth affected cell,y direction, PBC
    IF (dimYCell < 0) THEN
     dimYCell=dimYCell+iBox(2)
    ELSE IF (dimYCell >= iBox(2)) THEN
     dimYCell=dimYCell-iBox(2)
    END IF
    dimZCell=k				!Find kth affected cell,z direction, PBC
    IF (dimZCell < 0) THEN
     dimZCell=dimZCell+iBox(3)
    ELSE IF (dimZCell >= iBox(3)) THEN
     dimZCell=dimZCell-iBox(3)
    END IF

    CellN=dimXCell+iBox(1)*dimYCell+iBox(1)*iBox(2)*dimZCell+1
    !CellN denotes the affected cell given by i=dimXCell, j=dimYCell, k=dimZCell
    affectedCellN(CellN)=1

   END IF
  END DO
 END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE Affected_Cell
