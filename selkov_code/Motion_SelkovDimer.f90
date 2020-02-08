SUBROUTINE Motion_SelkovDimer
USE Global_SelkovDimer
IMPLICIT NONE

!The first of two subroutines to describe the MD motion of the solvent (Velocity Verlet) and dimer (RATTLE algorithm). This subroutine deals with the position of the solvent at the next time step r(t+h) and the position of the dimer at the next time step under the constraint of constant bond length (first Lagrange multiplier determined).  The next subroutine, Forces, calculates the force at r(t+h) and the velocities of both the solvent and dimer (and second Lagrange multiplier).

!!!!!!Define variables within subroutine
hIMfMDt2=(0.5d0*MDt*MDt)/Mf			!Value used a lot for calculations
hIMgMDt2=(0.5d0*MDt*MDt)/Mg			!Value used a lot for calculations
hIMiMDt2=(0.5d0*MDt*MDt)/Mi			!Value used a lot for calculations
hIMcMDt=(0.5d0*MDt)/Mc				!Value used a lot for calculations
hIMnMDt=(0.5d0*MDt)/Mn				!Value used a lot for calculations
IMc=1.d0/Mc					!Inverse Mc
IMn=1.d0/Mn					!Inverse Mn
IIMcIMn=1.d0/((IMc)+(IMn))			!Value used for first Lagr. mult.
IMDt=1.d0/MDt					!Inverse MDt

!!!!!!Set forces at time t to be previous forces
DO i=1,NT
  pFXsol(i)=FXsol(i)
  pFYsol(i)=FYsol(i)
  pFZsol(i)=FZsol(i)
END DO
DO i=1,2
 pFXdim(i)=FXdim(i)
 pFYdim(i)=FYdim(i)
 pFZdim(i)=FZdim(i)
END DO

!!!!!!Set old dimer pos. and vel. to previous time step
DO i=1,2
  prXdim(i)=Xdim(i)
  prYdim(i)=Ydim(i)
  prZdim(i)=Zdim(i)
  prVXdim(i)=VXdim(i)
  prVYdim(i)=VYdim(i)
  prVZdim(i)=VZdim(i)
END DO

!!!!!!Motion of the solvent molecules
!$OMP PARALLEL DO &
!$OMP SHARED(MDt,hIMfMDt2,hIMgMDt2,hIMiMDt2,Box,Xsol,Ysol,Zsol) &
!$OMP SHARED(VXsol,VYsol,VZsol,FXsol,FYsol,FZsol,FLAG) &
!$OMP PRIVATE(i) SCHEDULE(static)

DO i=1,NT
 IF (FLAG(i) == 0) THEN
  Xsol(i)=Xsol(i)+VXsol(i)*MDt+pFXsol(i)*hIMfMDt2	!Solvent F x pos. at next step
  Ysol(i)=Ysol(i)+VYsol(i)*MDt+pFYsol(i)*hIMfMDt2	!Solvent F y pos. at next step
  Zsol(i)=Zsol(i)+VZsol(i)*MDt+pFZsol(i)*hIMfMDt2	!Solvent F z pos. at next step
 ELSE IF (FLAG(i) == 1) THEN
  Xsol(i)=Xsol(i)+VXsol(i)*MDt+pFXsol(i)*hIMgMDt2	!Solvent G x pos. at next step
  Ysol(i)=Ysol(i)+VYsol(i)*MDt+pFYsol(i)*hIMgMDt2	!Solvent G y pos. at next step
  Zsol(i)=Zsol(i)+VZsol(i)*MDt+pFZsol(i)*hIMgMDt2	!Solvent G z pos. at next step
 ELSE IF (FLAG(i) == 3) THEN
  Xsol(i)=Xsol(i)+VXsol(i)*MDt+pFXsol(i)*hIMiMDt2	!In. solvent x pos. at next step
  Ysol(i)=Ysol(i)+VYsol(i)*MDt+pFYsol(i)*hIMiMDt2	!In. solvent G y pos. at next step
  Zsol(i)=Zsol(i)+VZsol(i)*MDt+pFZsol(i)*hIMiMDt2	!IN. solvent G z pos. at next step
 END IF
  Xsol(i)=Xsol(i)-Box(1)*FLOOR(Xsol(i)/Box(1))	!Periodic BC in x direction
  Ysol(i)=Ysol(i)-Box(2)*FLOOR(Ysol(i)/Box(2))	!Periodic BC in y direction
  Zsol(i)=Zsol(i)-Box(3)*FLOOR(Zsol(i)/Box(3))	!Periodic BC in z direction
END DO
!$OMP END PARALLEL DO

!!!!!!Motion of the dimer (following Andersen, 1983, Appendix C)
!Difference in position at prev. step (used for first Lagrange multiplier later)
dXdimP=prXdim(1)-prXdim(2)
dYdimP=prYdim(1)-prYdim(2)
dZdimp=prZdim(1)-prZdim(2)

!Periodic BC for diff. in pos. at pstep
IF (dXdimP > HBox(1)) THEN
 dXdimP=dXdimP-Box(1)
ELSE IF (dXdimP < -HBox(1)) THEN
 dXdimP=dXdimP+Box(1)
END IF
IF (dYdimP > HBox(2)) THEN
 dYdimP=dYdimP-Box(2)
ELSE IF (dYdimP < -HBox(2)) THEN
 dYdimP=dYdimP+Box(2)
END IF
IF (dZdimP > HBox(3)) THEN
 dZdimP=dZdimP-Box(3)
ELSE IF (dZdimP < -HBox(3)) THEN
 dZdimP=dZdimP+Box(3)
END IF

!The value for q without added Lagrange multiplier term
QXdim(1)=prVXdim(1)+hIMcMDt*pFXdim(1)
QYdim(1)=prVYdim(1)+hIMcMDt*pFYdim(1)
QZdim(1)=prVZdim(1)+hIMcMDt*pFZdim(1)
QXdim(2)=prVXdim(2)+hIMnMDt*pFXdim(2)
QYdim(2)=prVYdim(2)+hIMnMDt*pFYdim(2)
QZdim(2)=prVZdim(2)+hIMnMDt*pFZdim(2)

!Advance position by one time step
DO i=1,2
 Xdim(i)=prXdim(i)+MDt*QXdim(i)
 Ydim(i)=prYdim(i)+MDt*QYdim(i)
 Zdim(i)=prZdim(i)+MDt*QZdim(i)

 !Apply periodic BC (if particle exits sim. box, bring it back into box)
 Xdim(i)=Xdim(i)-Box(1)*FLOOR(Xdim(i)/Box(1))
 Ydim(i)=Ydim(i)-Box(2)*FLOOR(Ydim(i)/Box(2))
 Zdim(i)=Zdim(i)-Box(3)*FLOOR(Zdim(i)/Box(3))
END DO

Con1Acc=0

DO
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

 !Find distance between the two dimers and compare to bond length
 dXYZ2=dXdim*dXdim+dYdim*dYdim+dZdim*dZdim
 Con1diff=dXYZ2-Dcn*Dcn

 !If distance is within constraint, exit and move on; if not, repeat with Lagrange mul
 IF(DABS(Con1diff) < Con1Conv .OR. Con1Acc == Con1Max) EXIT
  !Langrange multiplier formula for G
  G=(0.5d0*Con1diff*IMDt*IIMcIMn)/(dXdim*dXdimP+dYdim*dYdimP+dZdim*dZdimP)
  !Modify q for constraint by subtracting or adding constraint force
  QXdim(1)=QXdim(1)-G*IMc*dXdimP
  QYdim(1)=QYdim(1)-G*IMc*dYdimP
  QZdim(1)=QZdim(1)-G*IMc*dZdimP
  QXdim(2)=QXdim(2)+G*IMn*dXdimP
  QYdim(2)=QYdim(2)+G*IMn*dYdimP
  QZdim(2)=QZdim(2)+G*IMn*dZdimP
  !Define new value for advancing to next time step (with constraint)
  DO i=1,2
   Xdim(i)=prXdim(i)+MDt*QXdim(i)
   Ydim(i)=prYdim(i)+MDt*QYdim(i)
   Zdim(i)=prZdim(i)+MDt*QZdim(i)
   !Apply periodic BC (if particle exits sim. box, bring it back into box)
   Xdim(i)=Xdim(i)-Box(1)*FLOOR(Xdim(i)/Box(1))
   Ydim(i)=Ydim(i)-Box(2)*FLOOR(Ydim(i)/Box(2))
   Zdim(i)=Zdim(i)-Box(3)*FLOOR(Zdim(i)/Box(3))
  END DO
 Con1Acc=Con1Acc+1
END DO

END SUBROUTINE Motion_SelkovDimer
