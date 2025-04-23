!--------------------------------------------

MODULE CONSTANT
!Define global constants
IMPLICIT NONE
SAVE

INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND (p=5)
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND (p=15)
INTEGER, PARAMETER :: INTL=SELECTED_INT_KIND (12)

REAL (KIND=DBL), PARAMETER :: PI=4.d0*atan(1.d0)        !Pi
COMPLEX (KIND=DBL), PARAMETER :: IMG_I = CMPLX(0.0,1.0,DBL) !imaginary unit

REAL (KIND=DBL), PARAMETER :: BOLTZ=1.3806485279d-23    !Boltzmann constant [J/K]
REAL (KIND=DBL), PARAMETER :: AFGDL=6.022140857d+23     !Avogadro constant [1/mol]
REAL (KIND=DBL), PARAMETER :: PLANCK=6.62606993489d-34  !Planck constant [J.s]
REAL (KIND=DBL), PARAMETER :: DIRAC=PLANCK/(2.d0*PI)    !reduced Planck constant $\hbar$ [J.s]
REAL (KIND=DBL), PARAMETER :: EV=1.602176634d-19        !1ev [J]
REAL (KIND=DBL), PARAMETER :: RYD=2.179872361103024d-18 !Rydberg unit of energy [J]
REAL (KIND=DBL), PARAMETER :: AMU=1.66053906892d-27     !atomic mass [kg]
REAL (KIND=DBL), PARAMETER :: BOhr=5.2917721054482d-11  !Bohr radius [m]

REAL (KIND=DBL), PARAMETER :: MIN_ERR=1.0d-16
REAL (KIND=DBL), PARAMETER :: MID_ERR=1.0d-10
REAL (KIND=DBL), PARAMETER :: MAX_ERR=1.0d-6

END MODULE CONSTANT

!-----------------------------------------------

MODULE GLOBAL_VARIABLE
!Define and initial global variables
USE CONSTANT
IMPLICIT NONE

REAL (KIND=DBL), SAVE :: TOL=1.0d-5                    !iteration tolerance
INTEGER, SAVE :: TMAX=10000                            !maximum iteration steps
INTEGER, SAVE :: NPOLE=20, NAZIM=40                    !number of discretized polar and azimuthal angles
INTEGER, SAVE :: DEG=1                                 !degree of polynomial of DG discretization
INTEGER, SAVE :: NDOF_TRI                              !number of degree of freedom in element (triangle): NDEG_TRI=(DEG+1)*(DEG+2)/2
INTEGER, SAVE :: NDOF_FC                               !number of degree of freedom on face: NDEG_FC=DEG+1
INTEGER, SAVE :: NP_TRI                                !number of quadrature points in element
INTEGER, SAVE :: NP_FC                                 !number of quadrature points on face

!Material properties 
REAL (KIND=DBL), SAVE :: Ce=2.0d4               !electron heat capacity at $T0$ [J/(m^3.K)]
REAL (KIND=DBL), SAVE :: Ve=1.36d6              !electron group velocity [m/s] 
REAL (KIND=DBL), SAVE :: Cp=2.5d6               !phonon heat capacity at $T0$ [J/(m^3.K)] 
REAL (KIND=DBL), SAVE :: Vp=2142.857            !phonon group velocity [m/s] 
REAL (KIND=DBL), SAVE :: TAU_e_p=2.595d-14      !relaxation time [s]
REAL (KIND=DBL), SAVE :: TAU_p_e=7.187d-13
REAL (KIND=DBL), SAVE :: TAU_p_p=7.187d-13

REAL (KIND=DBL), SAVE :: T0=300.d0              
!Reference quantities used to nondimensionalize!global reference temperature [K]
REAL (KIND=DBL), SAVE :: T_ref = 100.d0
REAL (KIND=DBL), SAVE :: L_ref = 1.d-8
REAL (KIND=DBL), SAVE :: Ve_ref= 1.d0
REAL (KIND=DBL), SAVE :: Vp_ref= 1.d0
REAL (KIND=DBL), SAVE :: Ce_ref= 1.d0
REAL (KIND=DBL), SAVE :: Cp_ref= 1.d0  ! we can use different reference value to nomoralize electron and phonon parameter

!Dimensionless parameters (all variables with '_s' represents corresponding scaled dimensionless parameters)
REAL (KIND=DBL), SAVE :: Kn_e_p=1.d0
REAL (KIND=DBL), SAVE :: Kn_p_e=1.d0
REAL (KIND=DBL), SAVE :: Kn_p_p=1.d0
REAL (KIND=DBL), SAVE :: Kn_p_C=1.d0
REAL (KIND=DBL), SAVE :: Ve_s=1.d0
REAL (KIND=DBL), SAVE :: Vp_s=1.d0
REAL (KIND=DBL), SAVE :: Ce_s=1.d0
REAL (KIND=DBL), SAVE :: Cp_s=1.d0

CHARACTER (LEN=50), SAVE :: FNAME_MSH='Gmsh.msh'                           !File name of spatial mesh
CHARACTER (LEN=50), SAVE :: FNAME_ELECTRON_BAND='ElectronBand.txt'         !File name for elctron band structure
CHARACTER (LEN=50), SAVE :: FNAME_PHONON_DISPERSION='PhononDispersion.txt' !File name for phonon dispersion relation
CHARACTER (LEN=50), SAVE :: FNAME_SCATTERING_MATRIX='ScatteringMatrix.txt' !File name of e-ph scattering matrix

INTEGER, SAVE :: NBC = 1                                    !# of boundary condition
INTEGER, SAVE :: ACCFLAG                                    !0-no acceleration, 1-acceleration
REAL (KIND=DBL), SAVE :: KN_THR=1.d0
REAL (KIND=DBL), SAVE :: gamma=1.d0                         !shape factor: H/L

CONTAINS
   SUBROUTINE Init_Global_Variables (fin)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: fin
   INTEGER :: ios

   NAMELIST /ITERATION/TOL,TMAX
   NAMELIST /GSIS/ACCFLAG,KN_THR
   NAMELIST /VELMSH/NPOLE,NAZIM
   NAMELIST /DG/DEG
   NAMELIST /FILENAME/FNAME_MSH,FNAME_ELECTRON_BAND,FNAME_PHONON_DISPERSION,FNAME_SCATTERING_MATRIX
   NAMELIST /PHYS/T0,gamma
   NAMELIST /FLOW/Ce_s,Cp_s,Ve_s,Vp_s,Kn_e_p,Kn_p_e,Kn_p_p ! If we don't read from materials, the dimensionless parameters can be assigned previously in control.in
   NAMELIST /N_BC/NBC

   WRITE (*,*) '**** Initial global variables ****'
   READ (fin, NML=ITERATION, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read ITERATION parameters: ', ios
   WRITE (*,100) 'Iteration Tolerance: ', TOL
   100 FORMAT (1X, A45, ES10.3)
   WRITE (*,101) 'Maximum iteration step: ', TMAX
   101 FORMAT (1X, A45, I10)
   WRITE (*,*)

   READ (fin, NML=GSIS, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read GSIS parameters: ', ios
   IF (ACCFLAG.EQ.0) WRITE (*,105) 'Iteration Scheme: ','CIS'
   IF (ACCFLAG.EQ.1) WRITE (*,105) 'Iteration Scheme: ','GSIS'
   105 FORMAT (1X, A45, A4)
   WRITE (*,*)

   READ (fin, NML=VELMSH, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read VELMESH parameters: ', ios
   NAZIM = (NAZIM/2)*2
   WRITE (*,102) '# of polar angles: ', NPOLE
   102 FORMAT (1X,A45,I5)
   WRITE (*,102) '# of azimuthal angles: ', NAZIM
   WRITE (*,*)

   READ (fin, NML=DG, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read DG parameters: ', ios
   NDOF_TRI=(DEG+1)*(DEG+2)/2
   NDOF_FC=DEG+1

   IF(DEG.EQ.1) NP_TRI = 3
   IF(DEG.EQ.2) NP_TRI = 6
   IF(DEG.EQ.3) NP_TRI = 12
   IF(DEG.EQ.4) NP_TRI = 16

   NP_FC=15
   WRITE (*,103) 'Degree of polynomial: ', DEG
   103 FORMAT (1X,A45,I5)
   WRITE (*,103) 'Degree of freedom in element: ', NDOF_TRI
   WRITE (*,103) 'Degree of freedom on face: ', NDOF_FC
   WRITE (*,103) 'Quadrature points in element: ', NP_TRI
   WRITE (*,103) 'Quadrature points on face: ', NP_FC
   WRITE (*,*)
   
   READ (fin, NML=FILENAME, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read FILENAME parameters: ', ios
   WRITE (*,109) 'File for spatial mesh: ', FNAME_MSH
   109 FORMAT (1X,A45,A40)
   WRITE (*,109) 'File for elctron band structure: ', FNAME_ELECTRON_BAND
   WRITE (*,109) 'File for phonon dispersion relation: ', FNAME_PHONON_DISPERSION
   WRITE (*,109) 'File for elctron-phonon scattering matrix: ', FNAME_SCATTERING_MATRIX
   WRITE (*,*)

   READ (fin, NML=PHYS, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read TEMP parameters: ', ios
   WRITE (*,108) 'Global reference temperature: ', T0
   WRITE (*,108) 'Geometry size factor H/L: ', gamma
   108 FORMAT (1X,A45,F10.3)
   WRITE (*,*)

   READ (fin, NML=FLOW, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read FLOW parameters:', ios
   Kn_p_C=1.d0/(1.d0/Kn_p_e+1.d0/Kn_p_p)

   READ (fin, NML=N_BC, IOSTAT=ios)
   IF (ios.NE.0) WRITE(*,*) 'Error to read number of boundary condtion: ', ios
   WRITE (*,102) '# of boundary conditions: ', NBC
   WRITE (*,*)


   END SUBROUTINE Init_Global_Variables

END MODULE GLOBAL_VARIABLE
