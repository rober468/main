MODULE Global_SelkovDimer
IMPLICIT NONE

!This module contains all the necessary global variables that are used within the main program and appropriate subroutines. Some of the subroutines have their own variables that are not used within the main program.

!Define kind value for precision
INTEGER, PARAMETER::SGL=4
INTEGER, PARAMETER::DBL=8
INTEGER(KIND=selected_int_kind(9))::seed
INTEGER(KIND=DBL)::i,A,B,C,D,E,F,H	!Loop index
INTEGER(KIND=DBL)::j			!Loop index
INTEGER(KIND=DBL)::k			!Loop index
INTEGER(KIND=DBL)::m			!Loop index
INTEGER(KIND=DBL)::p			!Loop index
REAL(KIND=DBL)::numb
REAL(KIND=DBL)::tempsol
REAL(KIND=DBL)::tempdim
REAL(KIND=DBL)::scl,ab,ac,ad
INTEGER(KIND=DBL)::space
REAL(KIND=DBL),PARAMETER::pi=4*DATAN(1.d0)

!Number of molecules
INTEGER(KIND=DBL)::Nf0			!Initial number of particle F
INTEGER(KIND=DBL)::Ng0			!Initial number of particle G
INTEGER(KIND=DBL)::Ni0			!Number of inert particles
INTEGER(KIND=DBL)::NT			!Total number of F and G particles

!Simulation box
REAL(KIND=DBL),DIMENSION(1:3)::Box	!Real length of simulation box;x=1,y=2,z=3
INTEGER(KIND=DBL),PARAMETER::Vol=64000	!Volume of simulation box
REAL(KIND=DBL),DIMENSION(1:3)::HBox	!Half length of simulation box;x=1,y=2,z=3
INTEGER(KIND=DBL),DIMENSION(1:3)::iBox	!Integer length of simulation box

!Dimer and Solvent Characteristics
REAL(KIND=DBL)::Mf			!Mass of one F particle
REAL(KIND=DBL)::Mg			!Mass of one G particle
REAL(KIND=DBL)::Mi			!Mass of one inert particle
REAL(KIND=DBL)::Mc			!Mass of C sphere
REAL(KIND=DBL)::Mn			!Mass of N sphere
REAL(KIND=DBL)::Rc			!Radius of C sphere
REAL(KIND=DBL)::Rn			!Radius of N sphere
REAL(KIND=DBL)::Dcn			!Internuclear distance
REAL(KIND=DBL)::HDcn			!Half of internuclear distance
REAL(KIND=DBL)::EfC			!Reduced LJ Potential of F with C
REAL(KIND=DBL)::EgC			!Reduced LJ Potential of G with C
REAL(KIND=DBL)::EiC			!Reduced LJ Potential of G with C
REAL(KIND=DBL)::EfN			!Reduced LJ Potential of F with N
REAL(KIND=DBL)::EgN			!Reduced LJ Potential of G with N
REAL(KIND=DBL)::EiN			!Reduced LJ Potential of F with N
REAL(KIND=DBL)::T			!Reduced temperature

!MPC and MD time steps, cell size, and rotation angle
REAL(KIND=DBL)::A0				!MPC cell size
REAL(KIND=DBL)::Ang				!MPC rotation angle
REAL(KIND=DBL)::MDt				!MD time step
REAL(KIND=DBL)::MPCt				!MPC collision time
INTEGER(KIND=DBL)::MPCtstep			!MPC step happens every nth step
INTEGER(KIND=DBL),PARAMETER::tstep=1000	!Total number of time steps
INTEGER(KIND=DBL)::cstep			!Current time step
INTEGER(KIND=DBL)::pstep			!Previous time step
INTEGER(KIND=DBL)::rcstep			!=(pstep/100)+1
INTEGER(KIND=DBL)::time				!Number of times main loop is run

!Dimer and solvent positions
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::Xdim	!Dimer spheres positions, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::Ydim	!Dimer spheres positions, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::Zdim	!Dimer spheres positions, z
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::prXdim	!Previous dimer spheres positions, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::prYdim	!Previous dimer spheres positions, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::prZdim	!Previous dimer spheres positions, z
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::Xsol	!Solvent positions, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::Ysol	!Solvent positions, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::Zsol	!Solvent positions, z
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::FLAG!Identity of solvent
REAL(KIND=DBL)::dXdim				!Difference in dimer positions, x
REAL(KIND=DBL)::dYdim				!Difference in dimer positions, y
REAL(KIND=DBL)::dZdim				!Difference in dimer positions, z
REAL(KIND=DBL)::dXdimp				!Prev.Difference in dimer positions, x
REAL(KIND=DBL)::dYdimp				!Prev.Difference in dimer positions, y
REAL(KIND=DBL)::dZdimp				!Prev.Difference in dimer positions, z
REAL(KIND=DBL)::sumVXsol			!Sum of initial velocities, x
REAL(KIND=DBL)::sumVYsol			!Sum of initial velocities, y
REAL(KIND=DBL)::sumVZsol			!Sum of initial velocities, z

!Initial position additional quantities used
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::R	!Pos.of solv. (1,?), pos.of dim. (2,?)
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::dR	!Diff. in dim. and solv. positions
REAL(KIND=DBL)::dRsq				!Square of dXYZ (radial dis. squared)
REAL(KIND=DBL)::x				!Random number for x
REAL(KIND=DBL)::y				!Random number for y
REAL(KIND=DBL)::z				!Random number for z

!Dimer and solvent velocities
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VXdim	!Dimer spheres velocities, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VYdim	!Dimer spheres velocities, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VZdim	!Dimer spheres velocities, z
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::prVXdim!Previous dimer spheres velocities, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::prVYdim!Previous dimer spheres velocities, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::prVZdim!Previous dimer spheres velocities, z
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VXsol!Solvent velocities, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VYsol!Solvent velocities, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VZsol!Solvent velocities, z

!Initial velocity additional quantities used
REAL(KIND=DBL)::varF				!Variance of F vel. Gaussian distr.
REAL(KIND=DBL)::varG				!Variance of G vel. Gaussian distr.
REAL(KIND=DBL)::VarI				!Variance of inert vel. Gauss. distr.
REAL(KIND=DBL)::mean				!Mean of vel. Gaussian distr.
REAL(KIND=DBL)::dVXsol				!Variable used for iSumdVXsol
REAL(KIND=DBL)::dVYsol				!Variable used for iSumdVYsol
REAL(KIND=DBL)::dVZsol				!Variable used for iSumdVZsol

!Forces on the solvent and dimer
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::FXsol!Force on solvent in x direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::FYsol!Force on solvent in y direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::FZsol!Force on solvent in z direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::FXdim	!Force on solvent in x direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::FYdim	!Force on solvent in y direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::FZdim	!Force on solvent in z direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::pFXsol	!Force on solvent in x direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::pFYsol	!Force on solvent in y direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::pFZsol	!Force on solvent in z direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::pFXdim	!Force on solvent in x direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::pFYdim	!Force on solvent in y direction
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::pFZdim	!Force on solvent in z direction

!Distribution of molecules to cells
INTEGER(KIND=DBL)::CellN			!Specify the cell in the sim. box
INTEGER(KIND=SGL),PARAMETER::MaxF=200	!Predicted maximum of F and G in each cell
INTEGER(KIND=SGL),PARAMETER::MaxG=200	!If over this amt., will mess up calculation
INTEGER(KIND=SGL),PARAMETER::MaxI=200
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::tFCellN !Total number of F solvent in cell N
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::tGCellN !Total number of G solvent in cell N
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::tICellN !Total number of I solvent in cell N
INTEGER(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::sFCellN !Label F in Cell N with i
INTEGER(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::sGCellN !Label G in Cell N with i
INTEGER(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::sICellN !Label I in Cell N with i

!Verlocity Verlet/RATTLE motion additional quantities
REAL(KIND=DBL)::hIMfMDt2			!Equals (0.5*(MDt)^2)/Mf
REAL(KIND=DBL)::hIMgMDt2			!Equals (0.5*(MDt)^2)/Mg
REAL(KIND=DBL)::hIMiMDt2			!Equals (0.5*(MDt)^2)/Mi
REAL(KIND=DBL)::hIMcMDt				!Equals (0.5*MDt)/Mc
REAL(KIND=DBL)::hIMnMDt				!Equals (0.5*MDt)/Mn
REAL(KIND=DBL)::IMDt				!Inverse MDt
REAL(KIND=DBL)::IMc				!Inverse Mc
REAL(KIND=DBL)::IMn				!Inverse Mn
REAL(KIND=DBL)::hIMfMDt				!Equals (0.5*MDt)/Mf
REAL(KIND=DBL)::hIMgMDt				!Equals (0.5*MDt)/Mg
REAL(KIND=DBL)::hIMiMDt				!Equals (0.5*MDt)/Mi
REAL(KIND=DBL)::IIMcIMn				!Value used for first Lagr. mult.
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::QXdim	!The q(i) in Andersen paper, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::QYdim	!The q(i) in Andersen paper, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::QZdim	!The q(i) in Andersen paper, z
INTEGER(KIND=DBL)::Con1Acc			!Number of iter. for first constraint
INTEGER(KIND=DBL),PARAMETER::Con1Max=200000000	!Max numb. of iter. for first constr.
REAL(KIND=DBL),PARAMETER::Con1Conv=1.d-6	!Tolerance from bond length
REAL(KIND=DBL)::Con1diff			!Difference of new and old bond length
INTEGER(KIND=DBL)::Con2Acc			!Number of iter. for first constraint
INTEGER(KIND=DBL),PARAMETER::Con2Max=200000000	!Max numb. of iter. for first constr.
REAL(KIND=DBL),PARAMETER::Con2Conv=1.d-8	!Tolerance from dot prod. constr.
REAL(KIND=DBL)::Con2diff			!Difference of new and old velocities
REAL(KIND=DBL)::dXYZ				!Distance between two dimers
REAL(KIND=DBL)::dXYZ2				!Distance between two dimers, squared
REAL(KIND=DBL)::G				!Lagrange multiplier for first constr.
REAL(KIND=DBL)::O				!Lagrange multiplier for sec. constr.
REAL(KIND=DBL)::dVXdim				!Difference in dimer velocities, x
REAL(KIND=DBL)::dVYdim				!Difference in dimer velocities, y
REAL(KIND=DBL)::dVZdim				!Difference in dimer velocities, z

!LJ potential additional quantities
REAL(KIND=DBL),PARAMETER::CUTOFF=2.d0**(1.d0/6.d0)	!Cutoff prefactor for LJ pot.
REAL(KIND=DBL)::CUTOFFc				!Cutoff for C sphere LJ pot.
REAL(KIND=DBL)::CUTOFFn				!Cutoff for N sphere LJ pot.
REAL(KIND=DBL)::CUTOFFc2			!Cut. C sphere LJ pot. squared
REAL(KIND=DBL)::CUTOFFn2			!Cut. N sphere LJ pot. squared
INTEGER(KIND=DBL)::lowNXdim			!Lower limit of affected cells by N, x
INTEGER(KIND=DBL)::upNXdim			!Upper limit of affected cells by N, x
INTEGER(KIND=DBL)::lowNYdim			!Lower limit of affected cells by N, y
INTEGER(KIND=DBL)::upNYdim			!Upper limit of affected cells by N, y
INTEGER(KIND=DBL)::lowNZdim			!Lower limit of affected cells by N, z
INTEGER(KIND=DBL)::upNZdim			!Upper limit of affected cells by N, z
INTEGER(KIND=DBL)::dimXCell			!Cell that dimer has influence on, x
INTEGER(KIND=DBL)::dimYCell			!Cell that dimer has influence on, y
INTEGER(KIND=DBL)::dimZCell			!Cell that dimer has influence on, z
INTEGER(KIND=DBL)::l				!Placeholder for calculation
REAL(KIND=DBL)::IdRsq				!Inverse of the square of dR
REAL(KIND=DBL)::IdR6				!Inverse of (dR)**6
REAL(KIND=DBL)::IdR8				!Inverse of (dR)**8
REAL(KIND=DBL)::IdR12				!Inverse of (dR)**12
REAL(KIND=DBL)::IdR14				!Inverse of (dR)**14
REAL(KIND=DBL)::Fbrack				!LJ force calc., part with radii
REAL(KIND=DBL)::Ubrack				!LJ pot. calc., part with radii
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::FpairX	!LJ force in x direction
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::FpairY	!LJ force in y direction
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::FpairZ	!LJ force in z direction
REAL(KIND=DBL)::Rn6				!Rn^6
REAL(KIND=DBL)::Rn12				!Rn^12
REAL(KIND=DBL)::r12Rn12				!12*Rn^12
REAL(KIND=DBL)::r6Rn6				!6*Rn^6
REAL(KIND=DBL)::Rc6				!Rc^6
REAL(KIND=DBL)::Rc12				!Rc^12
REAL(KIND=DBL)::r12Rc12				!12*Rc^12
REAL(KIND=DBL)::r6Rc6				!6*Rc^6
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::Potential	!Potential per particle
REAL(KIND=DBL)::TPotential			!Total potential
INTEGER(KIND=DBL)::lowCXdim			!Lower limit of affected cells by C, x
INTEGER(KIND=DBL)::upCXdim			!Upper limit of affected cells by C, x
INTEGER(KIND=DBL)::lowCYdim			!Lower limit of affected cells by C, y
INTEGER(KIND=DBL)::upCYdim			!Upper limit of affected cells by C, y
INTEGER(KIND=DBL)::lowCZdim			!Lower limit of affected cells by C, z
INTEGER(KIND=DBL)::upCZdim			!Upper limit of affected cells by C, z
REAL(KIND=DBL)::RPfg				!Probability of reaction from F to G
REAL(KIND=DBL)::RPgf				!Probability of reaction from G to F
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::affectedCellN	!CellN is within region touching sphere

!Reaction Rates
REAL(KIND=DBL)::k1			! A->F, incorporation of [A] in k1
REAL(KIND=DBL)::k2			! F->A, reverse of k1
REAL(KIND=DBL)::k3			! F+2G->3G
REAL(KIND=DBL)::k4			! 3G->F+2G, reverse of k3
REAL(KIND=DBL)::k5			! G->B
REAL(KIND=DBL)::k6			! B->G, reverse of k5
REAL(KIND=DBL)::a00			! a00 = a1+a2+a3+a4+a5+a6
REAL(KIND=DBL)::a1			!Probability of reaction 1
REAL(KIND=DBL)::a2			!Probability of reaction 2
REAL(KIND=DBL)::a3			!Probability of reaction 3
REAL(KIND=DBL)::a4			!Probability of reaction 4
REAL(KIND=DBL)::a5			!Probability of reaction 5
REAL(KIND=DBL)::a6			!Probability of reaction 6
REAL(KIND=DBL)::a0rand			!Random number between 0 and a0
REAL(KIND=DBL)::ProbMPCt		!Probability for a reaction during MPCt

!File writing
INTEGER(KIND=DBL)::Freq1		!Files for dim. written at nth time step
INTEGER(KIND=DBL)::Freq2		!Files for solv. F/G/I written at nth time step
INTEGER(KIND=DBL)::Freq3		!Files for E and Mom. written at nth time step
INTEGER(KIND=DBL)::Freq4		!Files for vel.on Dcn written at nth time step
INTEGER(KIND=DBL)::Freq5		!Files for tot. F/G written at nth time step
INTEGER(KIND=DBL)::Freq6		!Files for velocity flow field
CHARACTER(LEN=29)::SolventFPosFile	!Used to save dimer F information
CHARACTER(LEN=29)::SolventGPosFile	!Used to save dimer G information
CHARACTER(LEN=29)::SolventIPosFile	!Used to save dimer i information
CHARACTER(LEN=29)::SolventFVelFile	!Used to save dimer F information
CHARACTER(LEN=29)::SolventGVelFile	!Used to save dimer G information
CHARACTER(LEN=29)::SolventIVelFile	!Used to save dimer i information
CHARACTER(LEN=29)::SolventFForFile	!Used to save dimer F information
CHARACTER(LEN=29)::SolventGForFile	!Used to save dimer G information
CHARACTER(LEN=29)::SolventIForFile	!Used to save dimer i information

!Energy and momentum calculations
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::PXsol	!Momentum of solvent, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::PYsol	!Momentum of solvent, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::PZsol	!Momentum of solvent, z
REAL(KIND=DBL)::sumPXsol			!Sum of all momenta of solvent, x
REAL(KIND=DBL)::sumPYsol			!Sum of all momenta of solvent, y
REAL(KIND=DBL)::sumPZsol			!Sum of all momenta of solvent, z
REAL(KIND=DBL)::PXdim			!Momentum of dimer, x
REAL(KIND=DBL)::PYdim			!Momentum of dimer, y
REAL(KIND=DBL)::PZdim			!Momentum of dimer, z
REAL(KIND=DBL)::sumPX			!Total momentum, x
REAL(KIND=DBL)::sumPY			!Total momentum, y
REAL(KIND=DBL)::sumPZ			!Total momentum, z
REAL(KIND=DBL)::KEsol			!Total kinetic energy of solvent
REAL(KIND=DBL)::KEdim			!Kinetic energy of the dimer
REAL(KIND=DBL)::Ekin			!Total kinetic energy of system
REAL(KIND=DBL)::Epot			!Total potential energy of system
REAL(KIND=DBL)::Etot			!Total energy of system
REAL(KIND=DBL)::VXdimav			!Center of mass velocity for dimer, x
REAL(KIND=DBL)::VYdimav			!Center of mass velocity for dimer, y
REAL(KIND=DBL)::VZdimav			!Center of mass velocity for dimer, z
REAL(KIND=DBL)::NucVecX			!Unit vector on internuc. axis, from N to C, x
REAL(KIND=DBL)::NucVecY			!Unit vector on internuc. axis, from N to C, y
REAL(KIND=DBL)::NucVecZ			!Unit vector on internuc. axis, from N to C, z
REAL(KIND=DBL)::DotVZ			!Velocity along internuclear axis, from N to C

!Concentration Calculations
REAL(KIND=DBL)::TotF				!Total number of F
REAL(KIND=DBL)::TotG				!Total number of G
REAL(KIND=DBL)::TotI				!Total number of I
INTEGER(KIND=DBL)::iTotF			!Total number of F
INTEGER(KIND=DBL)::iTotG			!Total number of G
INTEGER(KIND=DBL)::iTotI			!Total number of I
REAL(KIND=DBL)::ConcF				!Concentration of F
REAL(KIND=DBL)::ConcG				!Concentration of G
REAL(KIND=DBL)::ConcI                           !Concentration of I
REAL(KIND=DBL)::tFBulk				!Total F in bulk
REAL(KIND=DBL)::tGBulk				!Total G in bulk
REAL(KIND=DBL)::tIBulk                          !Total I in bulk
REAL(KIND=DBL)::TBulkCells			!Total unaffected cells
REAL(KIND=DBL)::ConcFBulk			!Concentration of F in bulk
REAL(KIND=DBL)::ConcGBulk			!Concentration of G in bulk
REAL(KIND=DBL)::ConcIBulk

!!!!!!Define variables within subroutine
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VXCM		!C.M. pre-collision vel., x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VYCM		!C.M. pre-collision vel., y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::VZCM		!C.M. pre-collision vel., z
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::tCellshift	!Tot.solv.in cellN after shift
REAL(KIND=DBL)::ItCellshift				!Inverse of tCellshift
REAL(KIND=DBL)::Xshift					!Shift in x direction
REAL(KIND=DBL)::Yshift					!Shift in y direction
REAL(KIND=DBL)::Zshift					!Shift in z direction
REAL(KIND=DBL)::Xsolshift				!Shifted position, x
REAL(KIND=DBL)::Ysolshift				!Shifted position, y
REAL(KIND=DBL)::Zsolshift				!Shifted position, z
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::NewCellN	!ith particle's new cell N
REAL(KIND=DBL)::phi					!Azimuthal angle
REAL(KIND=DBL)::theta					!Polar angle
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::UVX		!Unit vector, x
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::UVY		!Unit vector, y
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::UVZ		!Unit vector, z
REAL(KIND=DBL)::dVXsolVXCM		!Diff. betw. CM vel. and ith particle vel., x
REAL(KIND=DBL)::dVYsolVYCM		!Diff. betw. CM vel. and ith particle vel., y
REAL(KIND=DBL)::dVZsolVZCM		!Diff. betw. CM vel. and ith particle vel., z
REAL(KIND=DBL)::ProjVonUV		!Magnitude of projection onto unit vector
REAL(KIND=DBL)::UVXproj			!Projection onto unit vector, x
REAL(KIND=DBL)::UVYproj			!Projection onto unit vector, y
REAL(KIND=DBL)::UVZproj			!Projection onto unit vector, z
REAL(KIND=DBL)::UVXortho		!Vector orthogonal to UVXproj, x
REAL(KIND=DBL)::UVYortho		!Vector orthogonal to UVYproj, y
REAL(KIND=DBL)::UVZortho		!Vector orthogonal to UVZproj, z
REAL(KIND=DBL)::UVXortho2		!Vector orthogonal to UVXproj and UVXortho, x
REAL(KIND=DBL)::UVYortho2		!Vector orthogonal to UVYproj and UVYortho, y
REAL(KIND=DBL)::UVZortho2		!Vector orthogonal to UVZproj and UVZortho, z
REAL(KIND=DBL),PARAMETER::rotang=pi/2	!Rotation angle
REAL(KIND=DBL),PARAMETER::sinrot=DSIN(rotang)		!Sine of rotation angle
REAL(KIND=DBL),PARAMETER::cosrot=DCOS(rotang)		!Cosine of rotation angle
REAL(KIND=DBL)::dVXsolVXCMrot	!Diff. betw. CM vel. and ith particle vel., x, rotated
REAL(KIND=DBL)::dVYsolVYCMrot	!Diff. betw. CM vel. and ith particle vel., y, rotated
REAL(KIND=DBL)::dVZsolVZCMrot	!Diff. betw. CM vel. and ith particle vel., z, rotated

!!!!!!Shell Concentration
REAL(KIND=8)::moveX,moveY,moveZ
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::SolCenterX,SolCenterY,SolCenterZ
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Rx,Ry,Rz,R2
REAL(KIND=8)::Tot3F,Tot3G,Tot4F,Tot4G,Tot6F,Tot6G,Tot8F,Tot8G,Tot10F
REAL(KIND=8)::Tot10G,Tot12F,Tot12G,Tot15F,Tot15G,Tot20F,Tot20G
REAL(KIND=8)::Tot3I,Tot4I,Tot6I,Tot8I,Tot10I,Tot12I,Tot15I,Tot20I
REAL(KIND=8)::Conc3F,Conc3G,Conc4F,Conc4G,Conc6F,Conc6G,Conc8F,Conc8G
REAL(KIND=8)::Conc10F,Conc10G,Conc12F,Conc12G,Conc15F,Conc15G,Conc20F,Conc20G
REAL(KIND=8)::Conc3I,Conc4I,Conc6I,Conc8I,Conc10I,Conc12I,Conc15I,Conc20I
REAL(KIND=8)::Vol2,Vol3,Vol4,Vol6,Vol8,Vol10,Vol12,Vol15,Vol20,ShellVol3
REAL(KIND=8)::ShellVol4,ShellVol6,ShellVol8,ShellVol10,ShellVol12,ShellVol15,ShellVol20

!!!!!!Velocity flow field
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Shift
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Xdimflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Ydimflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Zdimflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Xsolflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Ysolflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Zsolflow
REAL(KIND=8)::XdimNrelCflow
REAL(KIND=8)::YdimNrelCflow
REAL(KIND=8)::ZdimNrelCflow
REAL(KIND=8)::RdimNrelCflow
REAL(KIND=8)::IRdimNrelCflow
REAL(KIND=8)::thetaflow
REAL(KIND=8)::phiflow
REAL(KIND=8)::Xdimrotflow
REAL(KIND=8)::Ydimrotflow
REAL(KIND=8)::Zdimrotflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::XsolrelCflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::YsolrelCflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::ZsolrelCflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Xsolrotflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Ysolrotflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::Zsolrotflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::VXsolrotflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::VYsolrotflow
REAL(KIND=8),DIMENSION(:),ALLOCATABLE::VZsolrotflow
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE::CellVelX
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE::CellVelY
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE::n

!!!!!!Sort
INTEGER(KIND=DBL)::new
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newXsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newYsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newZsol
INTEGER(KIND=DBL),DIMENSION(:),ALLOCATABLE::newFLAG
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newVXsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newVYsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newVZsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newFXsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newFYsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newFZsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newpFXsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newpFYsol
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::newpFZsol

!!!!!!Time
INTEGER(KIND=SGL),DIMENSION(1:3)::now

END MODULE Global_SelkovDimer
