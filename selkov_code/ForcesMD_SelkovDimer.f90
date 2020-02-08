SUBROUTINE ForcesMD_SelkovDimer
USE Global_SelkovDimer
USE mt95
IMPLICIT NONE

!The second of two subroutines to describe to describe the MD motion of the solvent (Velocity Verlet) and dimer (RATTLE algorithm). This subroutine deals with calculating the forces on the dimer and solvent at the next time step. If the F solvent is also within the vicinity of the dimer (within cutoff region), then a reaction may occur and convert F to G. Lastly, the velocity at the next time step for the solvent particles is calculated using the constraint on the velocities.

!!!!!!Define variables within subroutine
Rn6=Rn**6
Rn12=Rn**12
r6Rn6=6.*Rn6
r12Rn12=12.*Rn12
Rc6=Rc**6
Rc12=Rc**12
r6Rc6=6.*Rc6
r12Rc12=12.*Rc12
hIMfMDt=(0.5d0*MDt)/Mf
hIMgMDt=(0.5d0*MDt)/Mg
hIMiMDt=(0.5d0*MDt)/Mi
TPotential=0.d0

!!!!!!Set forces to 0 at time t+h
FXsol=0.
FYsol=0.
FZsol=0.
Potential=0.d0
FpairX=0.
FpairY=0.
FpairZ=0.
FXdim=0.
FYdim=0.
FZdim=0.

!!!!!!LJ forces between F and G solvent particles with non-catalytic sphere

!$OMP PARALLEL DO &
!$OMP PRIVATE(i,p,R,dR,dRsq,IdRsq,IdR6,IdR8,IdR12,IdR14) &
!$OMP SHARED(Xdim,Ydim,Zdim,HBox,Box,CUTOFFn2,r12Rn12) &
!$OMP SHARED(r6Rn6,Rn12,Rn6,EfN,EgN,EiN) PRIVATE (Fbrack,Ubrack) &
!$OMP SCHEDULE(dynamic,1000)

DO i=1,NT
 R(1,1)=Xsol(i)
 R(1,2)=Ysol(i)
 R(1,3)=Zsol(i)
 R(2,1)=Xdim(2)
 R(2,2)=Ydim(2)
 R(2,3)=Zdim(2)

 DO p=1,3					!Periodic BC
  dR(p)=R(1,p)-R(2,p)				!dR vector is from solvent to dimer
  IF (dR(p) > HBox(p)) THEN
   dR(p)=dR(p)-Box(p)
  ELSE IF (dR(p) < -HBox(p)) THEN		
   dR(p)=dR(p)+Box(p)		
  END IF
 END DO 
 dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) 

 IF (dRsq <= CUTOFFn2) THEN			!If within cutoff region of LJ pot.
  IF (FLAG(i) == 0) THEN			!For F molecules, then:
   IdRsq=1.d0/dRsq				!Inverse of the square of dR
   IdR6=IdRsq**3				!Inverse of (dR)**6
   IdR8=IdR6*IdRsq				!Inverse of (dR)**8
   IdR12=IdR6*IdR6				!Inverse of (dR)**12
   IdR14=IdR12*IdRsq				!Inverse of (dR)**14
   Fbrack=4.d0*(r12Rn12*IdR14-r6Rn6*IdR8)	!LJ force calc., part with radii
   Ubrack=4.d0*(Rn12*IdR12-Rn6*IdR6)+1.d0	!LJ pot. calc., part with radii
   FpairX(i,2)=EfN*Fbrack*dR(1)			!LJ force in x direction
   FpairY(i,2)=EfN*Fbrack*dR(2)			!LJ force in y direction
   FpairZ(i,2)=EfN*Fbrack*dR(3)			!LJ force in z direction
   Potential(i)=EfN*Ubrack			!Addition of LJ Potential
   
   FXsol(i)=FpairX(i,2)			!Addition of LJ force on F solv., x
   FYsol(i)=FpairY(i,2)			!Addition of LJ force on F solv., y
   FZsol(i)=FpairZ(i,2)			!Addition of LJ force on F solv., z

  ELSE IF (FLAG(i) == 1) THEN			!For G molecules, then:
   IdRsq=1.d0/dRsq				!Inverse of the square of dR
   IdR6=IdRsq**3				!Inverse of (dR)**6
   IdR8=IdR6*IdRsq				!Inverse of (dR)**8
   IdR12=IdR6*IdR6				!Inverse of (dR)**12
   IdR14=IdR12*IdRsq				!Inverse of (dR)**14
   Fbrack=4.d0*(r12Rn12*IdR14-r6Rn6*IdR8)	!LJ force calc., part with radii
   Ubrack=4.d0*(Rn12*IdR12-Rn6*IdR6)+1.d0	!LJ pot. calc., part with radii
   FpairX(i,2)=EgN*Fbrack*dR(1)			!LJ force in x direction
   FpairY(i,2)=EgN*Fbrack*dR(2)			!LJ force in y direction
   FpairZ(i,2)=EgN*Fbrack*dR(3)			!LJ force in z direction
   Potential(i)=EgN*Ubrack			!Addition of LJ Potential
      
   FXsol(i)=FpairX(i,2)			!Addition of LJ force on F solv., x
   FYsol(i)=FpairY(i,2)			!Addition of LJ force on F solv., y
   FZsol(i)=FpairZ(i,2)			!Addition of LJ force on F solv., z

  ELSE IF (FLAG(i) == 3) THEN			!For inert molecules, then:
   IdRsq=1.d0/dRsq				!Inverse of the square of dR
   IdR6=IdRsq**3				!Inverse of (dR)**6
   IdR8=IdR6*IdRsq				!Inverse of (dR)**8
   IdR12=IdR6*IdR6				!Inverse of (dR)**12
   IdR14=IdR12*IdRsq				!Inverse of (dR)**14
   Fbrack=4.d0*(r12Rn12*IdR14-r6Rn6*IdR8)	!LJ force calc., part with radii
   Ubrack=4.d0*(Rn12*IdR12-Rn6*IdR6)+1.d0	!LJ pot. calc., part with radii
   FpairX(i,2)=EiN*Fbrack*dR(1)			!LJ force in x direction
   FpairY(i,2)=EiN*Fbrack*dR(2)			!LJ force in y direction
   FpairZ(i,2)=EiN*Fbrack*dR(3)			!LJ force in z direction
   Potential(i)=EiN*Ubrack			!Addition of LJ Potential
      
   FXsol(i)=FpairX(i,2)			!Addition of LJ force on F solv., x
   FYsol(i)=FpairY(i,2)			!Addition of LJ force on F solv., y
   FZsol(i)=FpairZ(i,2)			!Addition of LJ force on F solv., z
  END IF
 END IF
END DO
!$OMP END PARALLEL DO

!!!!!!LJ forces between F and G solvent particles with catalytic sphere

!$OMP PARALLEL DO &
!$OMP PRIVATE(i,p,R,dR,dRsq,IdRsq,IdR6,IdR8,IdR12,IdR14,Fbrack) &
!$OMP SHARED(Xdim,Ydim,Zdim,HBox,Box,CUTOFFc2,r12Rc12,r6Rc6,Rc12) &
!$OMP SHARED(Rc6,EfC,EgC,EiC,RPfg,RPgf) PRIVATE(Ubrack,x) &
!$OMP SCHEDULE(dynamic,1000)
DO i=1,NT
 R(1,1)=Xsol(i)
 R(1,2)=Ysol(i)
 R(1,3)=Zsol(i)
 R(2,1)=Xdim(1)
 R(2,2)=Ydim(1)
 R(2,3)=Zdim(1)

 DO p=1,3					!Periodic BC
  dR(p)=R(1,p)-R(2,p)				!dR vector is from solvent to dimer
  IF (dR(p) > HBox(p)) THEN
   dR(p)=dR(p)-Box(p)
  ELSE IF (dR(p) < -HBox(p)) THEN		
   dR(p)=dR(p)+Box(p)		
  END IF
 END DO
 dRsq= dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) 
 
 IF (dRsq <= CUTOFFc2) THEN			!If within cutoff region of LJ pot.
  IF (FLAG(i) == 0) THEN			!For F molecules:
   CALL genrand_real1(x)			!Call random number
   IF (RPfg > x) THEN				!If reaction F -> G occurs, then:
    IdRsq=1.d0/dRsq				!Inverse of the square of dR
    IdR6=IdRsq**3				!Inverse of (dR)**6
    IdR8=IdR6*IdRsq				!Inverse of (dR)**8
    IdR12=IdR6*IdR6				!Inverse of (dR)**12
    IdR14=IdR12*IdRsq				!Inverse of (dR)**14
    Fbrack=4.d0*(r12Rc12*IdR14-r6Rc6*IdR8)	!LJ force calc., part with radii
    Ubrack=4.d0*(Rc12*IdR12-Rc6*IdR6)+1.d0	!LJ pot. calc., part with radii
    FpairX(i,1)=EfC*Fbrack*dR(1)			!LJ force in x direction
    FpairY(i,1)=EfC*Fbrack*dR(2)			!LJ force in y direction
    FpairZ(i,1)=EfC*Fbrack*dR(3)			!LJ force in z direction
    Potential(i)=EfC*Ubrack		!Addition of LJ Potential
   
    FXsol(i)=FpairX(i,1)			!Addition of LJ force on F solv., x
    FYsol(i)=FpairY(i,1)			!Addition of LJ force on F solv., y
    FZsol(i)=FpairZ(i,1)			!Addition of LJ force on F solv., z

    FLAG(i)=1					!Convert F -> G
   ELSE IF (RPfg <= x) THEN			!If F -> G doesn't occur, then:
    IdRsq=1.d0/dRsq				!Inverse of the square of dR
    IdR6=IdRsq**3				!Inverse of (dR)**6
    IdR8=IdR6*IdRsq				!Inverse of (dR)**8
    IdR12=IdR6*IdR6				!Inverse of (dR)**12
    IdR14=IdR12*IdRsq				!Inverse of (dR)**14
    Fbrack=4.d0*(r12Rc12*IdR14-r6Rc6*IdR8)	!LJ force calc., part with radii
    Ubrack=4.d0*(Rc12*IdR12-Rc6*IdR6)+1.d0	!LJ pot. calc., part with radii
    FpairX(i,1)=EfC*Fbrack*dR(1)			!LJ force in x direction
    FpairY(i,1)=EfC*Fbrack*dR(2)			!LJ force in y direction
    FpairZ(i,1)=EfC*Fbrack*dR(3)			!LJ force in z direction
    Potential(i)=EfC*Ubrack		!Addition of LJ Potential
   
    FXsol(i)=FpairX(i,1)			!Addition of LJ force on F solv., x
    FYsol(i)=FpairY(i,1)			!Addition of LJ force on F solv., y
    FZsol(i)=FpairZ(i,1)			!Addition of LJ force on F solv., z
   END IF
  ELSE IF (FLAG(i) == 1) THEN			!If molecule is G, then:
   CALL genrand_real1(x)			!Call a random number
   IF (RPgf > x) THEN				!If reaction G -> F occurs, then:
    IdRsq=1.d0/dRsq				!Inverse of the square of dR
    IdR6=IdRsq**3				!Inverse of (dR)**6
    IdR8=IdR6*IdRsq				!Inverse of (dR)**8
    IdR12=IdR6*IdR6				!Inverse of (dR)**12
    IdR14=IdR12*IdRsq				!Inverse of (dR)**14
    Fbrack=4.d0*(r12Rc12*IdR14-r6Rc6*IdR8)	!LJ force calc., part with radii
    Ubrack=4.d0*(Rc12*IdR12-Rc6*IdR6)+1.d0	!LJ pot. calc., part with radii
    FpairX(i,1)=EgC*Fbrack*dR(1)			!LJ force in x direction
    FpairY(i,1)=EgC*Fbrack*dR(2)			!LJ force in y direction
    FpairZ(i,1)=EgC*Fbrack*dR(3)			!LJ force in z direction
    Potential(i)=EgC*Ubrack			!Addition of LJ Potential
   
    FXsol(i)=FpairX(i,1)			!Addition of LJ force on F solv., x
    FYsol(i)=FpairY(i,1)			!Addition of LJ force on F solv., y
    FZsol(i)=FpairZ(i,1)			!Addition of LJ force on F solv., z

    FLAG(i)=0 					!Convert G -> F
   ELSE IF (RPgf <= x) THEN			!If G -> F doesn't occur, then:
    IdRsq=1.d0/dRsq				!Inverse of the square of dR
    IdR6=IdRsq**3				!Inverse of (dR)**6
    IdR8=IdR6*IdRsq				!Inverse of (dR)**8
    IdR12=IdR6*IdR6				!Inverse of (dR)**12
    IdR14=IdR12*IdRsq				!Inverse of (dR)**14
    Fbrack=4.d0*(r12Rc12*IdR14-r6Rc6*IdR8)	!LJ force calc., part with radii
    Ubrack=4.d0*(Rc12*IdR12-Rc6*IdR6)+1.d0	!LJ pot. calc., part with radii
    FpairX(i,1)=EgC*Fbrack*dR(1)			!LJ force in x direction
    FpairY(i,1)=EgC*Fbrack*dR(2)			!LJ force in y direction
    FpairZ(i,1)=EgC*Fbrack*dR(3)			!LJ force in z direction
    Potential(i)=EgC*Ubrack			!Addition of LJ Potential
   
    FXsol(i)=FpairX(i,1)			!Addition of LJ force on F solv., x
    FYsol(i)=FpairY(i,1)			!Addition of LJ force on F solv., y
    FZsol(i)=FpairZ(i,1)			!Addition of LJ force on F solv., z
   END IF
  ELSE IF (FLAG(i) == 3) THEN
   IdRsq=1.d0/dRsq				!Inverse of the square of dR
   IdR6=IdRsq**3				!Inverse of (dR)**6
   IdR8=IdR6*IdRsq				!Inverse of (dR)**8
   IdR12=IdR6*IdR6				!Inverse of (dR)**12
   IdR14=IdR12*IdRsq				!Inverse of (dR)**14
   Fbrack=4.d0*(r12Rc12*IdR14-r6Rc6*IdR8)	!LJ force calc., part with radii
   Ubrack=4.d0*(Rc12*IdR12-Rc6*IdR6)+1.d0	!LJ pot. calc., part with radii
   FpairX(i,1)=EiC*Fbrack*dR(1)			!LJ force in x direction
   FpairY(i,1)=EiC*Fbrack*dR(2)			!LJ force in y direction
   FpairZ(i,1)=EiC*Fbrack*dR(3)			!LJ force in z direction
   Potential(i)=EiC*Ubrack			!Addition of LJ Potential
      
   FXsol(i)=FpairX(i,1)			!Addition of LJ force on F solv., x
   FYsol(i)=FpairY(i,1)			!Addition of LJ force on F solv., y
   FZsol(i)=FpairZ(i,1)			!Addition of LJ force on F solv., z
  END IF
 END IF
END DO
!$OMP END PARALLEL DO

DO i=1,NT
 IF (FLAG(i) == 0 .OR. FLAG(i) == 1 .OR. FLAG(i) == 3) THEN
  TPotential=TPotential+Potential(i)
  DO j=1,2
   FXdim(j)=FXdim(j)-FpairX(i,j)		!Addition of LJ force on N dim., x
   FYdim(j)=FYdim(j)-FpairY(i,j)		!Addition of LJ force on N dim., y
   FZdim(j)=FZdim(j)-FpairZ(i,j)		!Addition of LJ force on N dim., z
  END DO
 END IF
END DO

!!!!!!Velocity of solvent molecules
!Find velocities according to velocity Verlet algorithm

!$OMP PARALLEL DO &
!$OMP SHARED(hIMfMDT,hIMgMDT,hIMiMDT,NT,FLAG,VXsol,VYsol,VZsol) &
!$OMP SHARED(pFXsol,FXsol,pFYsol,FYsol,pFZsol,FZsol) &
!$OMP PRIVATE(i) 
DO i=1,NT
 IF (FLAG(i) == 0) THEN				!If F molecule, then:
  VXsol(i)=VXsol(i)+(pFXsol(i)+FXsol(i))*hIMfMDt
  VYsol(i)=VYsol(i)+(pFYsol(i)+FYsol(i))*hIMfMDt
  VZsol(i)=VZsol(i)+(pFZsol(i)+FZsol(i))*hIMfMDt
 ELSE IF (FLAG(i) == 1) THEN			!If G molecule, then:
  VXsol(i)=VXsol(i)+(pFXsol(i)+FXsol(i))*hIMgMDt
  VYsol(i)=VYsol(i)+(pFYsol(i)+FYsol(i))*hIMgMDt
  VZsol(i)=VZsol(i)+(pFZsol(i)+FZsol(i))*hIMgMDt
 ELSE IF (FLAG(i) == 3) THEN			!If inert mol., then:
  VXsol(i)=VXsol(i)+(pFXsol(i)+FXsol(i))*hIMiMDt
  VYsol(i)=VYsol(i)+(pFYsol(i)+FYsol(i))*hIMiMDt
  VZsol(i)=VZsol(i)+(pFZsol(i)+FZsol(i))*hIMiMDt
 END IF
END DO
!$OMP END PARALLEL DO

!!!!!!Velocity of Dimer(following Andersen, 1983, Appendix C)
!Find value of velocity without Langrange multiplier
VXdim(1)=QXdim(1)+hIMcMDt*FXdim(1)
VYdim(1)=QYdim(1)+hIMcMDt*FYdim(1)
VZdim(1)=QZdim(1)+hIMcMDt*FZdim(1)
VXdim(2)=QXdim(2)+hIMnMDt*FXdim(2)
VYdim(2)=QYdim(2)+hIMnMDt*FYdim(2)
VZdim(2)=QZdim(2)+hIMnMDt*FZdim(2)

Con2Acc=0						!Begin with zeroeth time

DO
 !Calculate difference in dimer velocities
 dVXdim=VXdim(1)-VXdim(2)
 dVYdim=VYdim(1)-VYdim(2)
 dVZdim=VZdim(1)-VZdim(2)

 !Check if dot product=0 is within acceptable range
 Con2diff=dVXdim*dXdim+dVYdim*dYdim+dVZdim*dZdim

 IF (DABS(Con2diff) < Con2Conv .OR. Con2Acc == Con2Max) EXIT
  O=(IIMcIMn*Con2diff)/(Dcn*Dcn)				!Calculate Lagrange mult.

  VXdim(1)=VXdim(1)-O*dXdim*IMc
  VYdim(1)=VYdim(1)-O*dYdim*IMc
  VZdim(1)=VZdim(1)-O*dZdim*IMc
  VXdim(2)=VXdim(2)+O*dXdim*IMn
  VYdim(2)=VYdim(2)+O*dYdim*IMn
  VZdim(2)=VZdim(2)+O*dZdim*IMn

  Con2Acc=Con2Acc+1
END DO

END SUBROUTINE ForcesMD_SelkovDimer
