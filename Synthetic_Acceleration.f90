MODULE ACCELERATION
USE GLOBAL_VARIABLE
USE SPATIAL_GRID
USE NODAL_FUNCTION
USE INTEGRATIONS
USE BOUNDARY_CONDITION
USE VELOCITY_GRID
USE VELOCITY_DISTRIBUTION_FUNCTION
IMPLICIT NONE

!Some variables emerged in dimensionless macroscopic equations
REAL (KIND=DBL), SAVE :: kappa_e, kappa_p, Ge, Gp, Ke_e, Ke_p, Kp_e, Kp_p
REAL (KIND=DBL), SAVE :: weT, wpT

REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:), SAVE :: UQ !DIMENSION(6*NDOF_TRI,N_TRIS), electron and ponon temperature, heatflux (local problem)
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:), SAVE :: U_TRACE !DIMENSION(6*NDOF_FC,N_FCS), traces of macro properties on interior faces (global problem)
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:), SAVE :: ST !DIMENSION(6), stabilization factors
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:,:), SAVE :: inv_AA_SOL !DIMENSION(6*NDOF_TRI, 6*NDOF_TRI, N_TRIS) coefficient matrix of UQ in local problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:,:,:), SAVE :: AA_TRACE !DIMENSION(6*NDOF_TRI, 6*NDOF_FC, 3, N_TRIS) coefficient matrix of traces in local problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:), SAVE :: AA_SRC !DIMENSION(6*NDOF_TRI, N_TRIS) source term

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: KKA ! DIMENSION(N_FCS*NDOF_FC*6, N_FCS*NDOF_FC*6), global matrix in global problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: AA_TMP, D1A_TMP, D2A_TMP
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: BA_SOL !DIMENSION(6*NDOF_FC,6*NDOF_TRI,3,N_TRIS), coefficient matrix for the field variables in global problem
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: KKAcomp, C1A_TMP, C2A_TMP
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: IKKA, JKKA, JOBA

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: FFA, HHA !R.H.S of global problem, DIMENSION(N_FCS*6*NDOF_FC)

!***********************************************************************************************
!pardiso
INTEGER, SAVE :: MAXFCT, MTYPE, PHASE, ERROR, MSGLVL
INTEGER*8, SAVE :: PT(64)
INTEGER, SAVE :: IPARM(64)
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: PERM

CONTAINS

SUBROUTINE Init_Acceleration_New ()
IMPLICIT NONE
EXTERNAL :: DGETRF, DGETRI
INTEGER :: I, M, L, IL, LWORK, INFO, TRID1, TRID2, LFCID1, LFCID2, GFCID1, GFCID2, GFCID,BCID, PAIR, NNZ
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: AA, BB, CC, DD, OO, PP, WW
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: AA_SOL
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:) :: WORK_A, KKAcomp_tmp
INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV_A, JKKA_tmp
REAL (KIND=DBL) :: nx, ny


kappa_e = Ve_s**2.d0*Ce_s*Kn_e_p/3.d0
kappa_p = Vp_s**2.d0*Cp_s*Kn_p_C/3.d0
Ge   = Cp_s / (Kn_p_e*Ce_s + Kn_e_p*Cp_s)
Gp   = Ce_s / (Kn_p_e*Ce_s + Kn_e_p*Cp_s)
weT  = Kn_p_e/(Kn_p_e*Ce_s + Kn_e_p*Cp_s)
wpT  = Kn_e_p/(Kn_p_e*Ce_s + Kn_e_p*Cp_s)
Ke_e = kappa_e*weT
Ke_p = kappa_e*wpT
Kp_e = kappa_p*weT*Kn_p_C/Kn_p_e
Kp_p = kappa_p*wpT*Kn_p_C/Kn_p_e+kappa_p*Kn_p_C/Cp_s/Kn_p_p

ALLOCATE(UQ(6*NDOF_TRI,N_TRIS),U_TRACE(6*NDOF_FC,N_FCS),ST(6))
UQ = 0.d0
U_TRACE = 0.d0
ST(1) = 1.d3 !/Hmin!/DELTA_REF
ST(2) = 1.d3 !/Hmin!/DELTA_REF
ST(3) = 1.d3 !/Hmin!/DELTA_REF
ST(4) = 1.d3 !/Hmin!/DELTA_REF
ST(5) = 1.d3 !/Hmin!/DELTA_REF
ST(6) = 1.d3 !/Hmin!/DELTA_REF

ALLOCATE (inv_AA_SOL(6*NDOF_TRI,6*NDOF_TRI,N_TRIS))
inv_AA_SOL = 0.d0

ALLOCATE(AA(NDOF_TRI,NDOF_TRI),BB(NDOF_TRI,NDOF_TRI),CC(NDOF_TRI,NDOF_TRI))
ALLOCATE(DD(NDOF_TRI,NDOF_TRI),OO(NDOF_TRI,NDOF_TRI),PP(NDOF_TRI,NDOF_TRI))
ALLOCATE(WORK_A(6*NDOF_TRI),IPIV_A(6*NDOF_TRI),AA_SOL(6*NDOF_TRI,6*NDOF_TRI))

LWORK = 6*NDOF_TRI

!$OMP PARALLEL DO PRIVATE(I,M,L,IL,INFO,AA,BB,CC,DD,PP,OO,AA_SOL,WORK_A,IPIV_A)
DO I = 1, N_TRIS
    AA = INT_NODFUNC_TRI_TRI_X(I,:,:)
    BB = INT_NODFUNC_TRI_TRI_Y(I,:,:)
    DO M = 1, NDOF_TRI
        DO L = 1, NDOF_TRI
            CC(M,L) = AA(L,M)
            DD(M,L) = BB(L,M)
        END DO
    END DO
    PP = 0.d0
    DO IL = 1, 3
        PP = PP + INT_NODFUNC_TRI_FC(I,IL,:,:)
    END DO
    OO = INT_NODFUNC_TRI_TRI(:,:,I)

    AA_SOL = 0.d0
    !eE
    AA_SOL(1:NDOF_TRI,            1:  NDOF_TRI) = ST(1)*PP + Ge*OO
    AA_SOL(1:NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = -Gp*OO
    AA_SOL(1:NDOF_TRI, 2*NDOF_TRI+1:3*NDOF_TRI) = CC
    AA_SOL(1:NDOF_TRI, 3*NDOF_TRI+1:4*NDOF_TRI) = DD

    !pE
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI,            1:  NDOF_TRI) = -Ge*OO
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = ST(2)*PP + Gp*OO
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI, 4*NDOF_TRI+1:5*NDOF_TRI) = CC
    AA_SOL(NDOF_TRI+1:2*NDOF_TRI, 5*NDOF_TRI+1:6*NDOF_TRI) = DD

    !eQx
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI,            1:  NDOF_TRI) = -Ke_e*AA
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = -Ke_p*AA
    AA_SOL(2*NDOF_TRI+1:3*NDOF_TRI, 2*NDOF_TRI+1:3*NDOF_TRI) = ST(3)*PP + OO


    !eQy
    AA_SOL(3*NDOF_TRI+1:4*NDOF_TRI,            1:  NDOF_TRI) = -Ke_e*BB
    AA_SOL(3*NDOF_TRI+1:4*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = -Ke_p*BB
    AA_SOL(3*NDOF_TRI+1:4*NDOF_TRI, 3*NDOF_TRI+1:4*NDOF_TRI) = ST(4)*PP + OO

    !pQx
    AA_SOL(4*NDOF_TRI+1:5*NDOF_TRI,            1:  NDOF_TRI) = -Kp_e*AA
    AA_SOL(4*NDOF_TRI+1:5*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = -Kp_p*AA
    AA_SOL(4*NDOF_TRI+1:5*NDOF_TRI, 4*NDOF_TRI+1:5*NDOF_TRI) = ST(5)*PP + OO

    !pQy
    AA_SOL(5*NDOF_TRI+1:6*NDOF_TRI,            1:  NDOF_TRI) = -Kp_e*BB
    AA_SOL(5*NDOF_TRI+1:6*NDOF_TRI,   NDOF_TRI+1:2*NDOF_TRI) = -Kp_p*BB
    AA_SOL(5*NDOF_TRI+1:6*NDOF_TRI, 5*NDOF_TRI+1:6*NDOF_TRI) = ST(6)*PP + OO

    WORK_A = 0.d0
    IPIV_A = 0
    CALL DGETRF (LWORK,LWORK,AA_SOL,LWORK,IPIV_A,INFO)
    IF (INFO.NE.0) WRITE (*,*) 'Warning LU factorization error: ',INFO
    CALL DGETRI (LWORK,AA_SOL,LWORK,IPIV_A,WORK_A,LWORK,INFO)
    IF (INFO.NE.0) WRITE (*,*) 'Warning Inverse Matrix error: ',INFO
    inv_AA_SOL(:,:,I) = AA_SOL

END DO
!$OMP END PARALLEL DO

DEALLOCATE(AA,BB,CC,DD,PP,OO,AA_SOL,WORK_A,IPIV_A)

!*******************************************************************************
ALLOCATE(AA_TRACE(6*NDOF_TRI,6*NDOF_FC,3,N_TRIS))
ALLOCATE(WW(NDOF_TRI,NDOF_FC))
AA_TRACE = 0.d0

!$OMP PARALLEL DO PRIVATE(I,IL,WW,nx,ny)
DO I = 1, N_TRIS
    DO IL = 1, 3
        WW = INT_NODFUNC_TRI_FC_FC(:,:,I,IL)
        nx = TRIANGLES_INF(I,2*IL)
        ny = TRIANGLES_INF(I,2*IL+1)

        !eE
        AA_TRACE(1:NDOF_TRI,           1:  NDOF_FC, IL, I) = ST(1)*WW

        !pE
        AA_TRACE(NDOF_TRI+1:2*NDOF_TRI,   NDOF_FC+1:2*NDOF_FC, IL, I) = ST(2)*WW

        !eQx
        AA_TRACE(2*NDOF_TRI+1:3*NDOF_TRI,           1:  NDOF_FC, IL, I) = -Ke_e*nx*WW
        AA_TRACE(2*NDOF_TRI+1:3*NDOF_TRI,   NDOF_FC+1:2*NDOF_FC, IL, I) = -Ke_p*nx*WW
        AA_TRACE(2*NDOF_TRI+1:3*NDOF_TRI, 2*NDOF_FC+1:3*NDOF_FC, IL, I) = ST(3)*WW
    
        !eQy
        AA_TRACE(3*NDOF_TRI+1:4*NDOF_TRI,           1:  NDOF_FC, IL, I) = -Ke_e*ny*WW
        AA_TRACE(3*NDOF_TRI+1:4*NDOF_TRI,   NDOF_FC+1:2*NDOF_FC, IL, I) = -Ke_p*ny*WW
        AA_TRACE(3*NDOF_TRI+1:4*NDOF_TRI, 3*NDOF_FC+1:4*NDOF_FC, IL, I) = ST(4)*WW

        !pQx
        AA_TRACE(4*NDOF_TRI+1:5*NDOF_TRI,           1:  NDOF_FC, IL, I) = -Kp_e*nx*WW
        AA_TRACE(4*NDOF_TRI+1:5*NDOF_TRI,   NDOF_FC+1:2*NDOF_FC, IL, I) = -Kp_p*nx*WW
        AA_TRACE(4*NDOF_TRI+1:5*NDOF_TRI, 4*NDOF_FC+1:5*NDOF_FC, IL, I) = ST(5)*WW

        !pQy
        AA_TRACE(5*NDOF_TRI+1:6*NDOF_TRI,           1:  NDOF_FC, IL, I) = -Kp_e*ny*WW
        AA_TRACE(5*NDOF_TRI+1:6*NDOF_TRI,   NDOF_FC+1:2*NDOF_FC, IL, I) = -Kp_p*ny*WW
        AA_TRACE(5*NDOF_TRI+1:6*NDOF_TRI, 5*NDOF_FC+1:6*NDOF_FC, IL, I) = ST(6)*WW

        AA_TRACE(:,:,IL,I) = MATMUL(inv_AA_SOL(:,:,I),AA_TRACE(:,:,IL,I))
    END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(WW)

!**********************************************************************************
ALLOCATE(AA_SRC(6*NDOF_TRI, N_TRIS))
AA_SRC = 0.d0

!***********************************************************************************
ALLOCATE(KKA(NDOF_FC*6,N_FCS*NDOF_FC*6))
ALLOCATE(BA_SOL(6*NDOF_FC,6*NDOF_TRI,3,N_TRIS))
ALLOCATE(KKAcomp_tmp(5*N_FCS*6*NDOF_FC*6*NDOF_FC))
ALLOCATE(IKKA(N_FCS*6*NDOF_FC+1))
ALLOCATE(JKKA_tmp(5*N_FCS*6*NDOF_FC*6*NDOF_FC))
ALLOCATE(FFA(N_FCS*6*NDOF_FC),HHA(N_FCS*6*NDOF_FC))

ALLOCATE(BB(NDOF_FC,NDOF_TRI))

KKA = 0.d0
BA_SOL = 0.d0
IKKA = 0
JKKA_tmp = 0
KKAcomp_tmp = 0.d0
FFA = 0.d0
HHA = 0.d0

!*****************************************************************************
!$OMP PARALLEL DO PRIVATE(I,IL,M,L,BB,nx,ny,GFCID,BCID)
DO I = 1, N_TRIS
    DO IL = 1, 3
        DO M = 1, NDOF_FC
            DO L = 1, NDOF_TRI
                BB(M,L) = INT_NODFUNC_TRI_FC_FC(L,M,I,IL)
            END DO
        END DO
        nx = TRIANGLES_INF(I,2*IL)
        ny = TRIANGLES_INF(I,2*IL+1)
        GFCID = TRIANGLES_TAG(I,3+IL)
        BCID = FACES_TAG(GFCID,3)

        IF (BCID.EQ.0) THEN
            !interior 
            BA_SOL(          1:  NDOF_FC,           1:  NDOF_TRI, IL, I) = 0.5d0*BB
            BA_SOL(          1:  NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = nx/2.d0/ST(1)*BB
            BA_SOL(          1:  NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = ny/2.d0/ST(1)*BB

            BA_SOL(  NDOF_FC+1:2*NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = 0.5d0*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = nx/2.d0/ST(2)*BB
            BA_SOL(  NDOF_FC+1:2*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = ny/2.d0/ST(2)*BB

            BA_SOL(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = 0.5d0*BB

            BA_SOL(3*NDOF_FC+1:4*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = 0.5d0*BB

            BA_SOL(4*NDOF_FC+1:5*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = 0.5d0*BB

            BA_SOL(5*NDOF_FC+1:6*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = 0.5d0*BB
        END IF
        IF (BCID.GT.0) THEN
            IF (BC_TYP(BCID).EQ.3) THEN
                !Periodic BC
                BA_SOL(          1:  NDOF_FC,           1:  NDOF_TRI, IL, I) = 0.5d0*BB
                BA_SOL(          1:  NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = nx/2.d0/ST(1)*BB
                BA_SOL(          1:  NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = ny/2.d0/ST(1)*BB

                BA_SOL(  NDOF_FC+1:2*NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = 0.5d0*BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = nx/2.d0/ST(2)*BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = ny/2.d0/ST(2)*BB

                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = 0.5d0*BB

                BA_SOL(3*NDOF_FC+1:4*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = 0.5d0*BB

                BA_SOL(4*NDOF_FC+1:5*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = 0.5d0*BB

                BA_SOL(5*NDOF_FC+1:6*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = 0.5d0*BB
            END IF
            IF (BC_TYP(BCID).EQ.1) THEN
                !eletron thermarlising wall
                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = BB

                BA_SOL(3*NDOF_FC+1:4*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = BB
            END IF
            IF (BC_TYP_P(BCID).EQ.1) THEN
                !phonon thermarlising wall     
                BA_SOL(4*NDOF_FC+1:5*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = BB
    
                BA_SOL(5*NDOF_FC+1:6*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = BB
            END IF

            IF (BC_TYP(BCID).EQ.2) THEN
                !eletron adiabatic boundary
                BA_SOL(          1:  NDOF_FC,           1:  NDOF_TRI, IL, I) = BB
                BA_SOL(          1:  NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = nx/ST(1)*BB
                BA_SOL(          1:  NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = ny/ST(1)*BB

                BA_SOL(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_TRI+1:3*NDOF_TRI, IL, I) = BB

                BA_SOL(3*NDOF_FC+1:4*NDOF_FC,3*NDOF_TRI+1:4*NDOF_TRI, IL, I) = BB
            END IF
            IF (BC_TYP_P(BCID).EQ.2) THEN
                !phonon adiabatic boundary
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,  NDOF_TRI+1:2*NDOF_TRI, IL, I) = BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = nx/ST(2)*BB
                BA_SOL(  NDOF_FC+1:2*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = ny/ST(2)*BB
    
                BA_SOL(4*NDOF_FC+1:5*NDOF_FC,4*NDOF_TRI+1:5*NDOF_TRI, IL, I) = BB
    
                BA_SOL(5*NDOF_FC+1:6*NDOF_FC,5*NDOF_TRI+1:6*NDOF_TRI, IL, I) = BB
            END IF
        END IF
    END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(BB)

ALLOCATE(D1A_TMP(6*NDOF_FC,6*NDOF_FC),D2A_TMP(6*NDOF_FC,6*NDOF_FC))
ALLOCATE(AA(NDOF_FC,NDOF_FC),AA_TMP(6*NDOF_FC,6*NDOF_FC))
ALLOCATE(C1A_TMP(6*NDOF_FC),C2A_TMP(6*NDOF_FC))

NNZ = 0
DO I = 1, N_FCS
    KKA = 0.d0

    AA = INT_NODFUNC_FC_FC(:,:,I)
    AA_TMP(          1:  NDOF_FC,          1:  NDOF_FC) = AA
    AA_TMP(  NDOF_FC+1:2*NDOF_FC,  NDOF_FC+1:2*NDOF_FC) = AA
    AA_TMP(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_FC+1:3*NDOF_FC) = AA
    AA_TMP(3*NDOF_FC+1:4*NDOF_FC,3*NDOF_FC+1:4*NDOF_FC) = AA
    AA_TMP(4*NDOF_FC+1:5*NDOF_FC,4*NDOF_FC+1:5*NDOF_FC) = AA
    AA_TMP(5*NDOF_FC+1:6*NDOF_FC,5*NDOF_FC+1:6*NDOF_FC) = AA

    BCID = FACES_TAG(I,3)

    IF (BCID.EQ.0) THEN
        !interior faces
        TRID1 = FACES_TAG(I,5)
        TRID2 = FACES_TAG(I,6)
        LFCID1 = FACES_TAG(I,7)
        LFCID2 = FACES_TAG(I,8)

        KKA(1:6*NDOF_FC,1+(I-1)*6*NDOF_FC:I*6*NDOF_FC) = AA_TMP
        DO IL = 1,3
            !interior
            GFCID1 = TRIANGLES_TAG(TRID1,3+IL)
            GFCID2 = TRIANGLES_TAG(TRID2,3+IL)

            D1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_TRACE(:,:,IL,TRID1))
            KKA(1:6*NDOF_FC, 1+(GFCID1-1)*6*NDOF_FC:GFCID1*6*NDOF_FC) = &
                KKA(1:6*NDOF_FC, 1+(GFCID1-1)*6*NDOF_FC:GFCID1*6*NDOF_FC) - D1A_TMP

            D2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_TRACE(:,:,IL,TRID2))
            KKA(1:6*NDOF_FC, 1+(GFCID2-1)*6*NDOF_FC:GFCID2*6*NDOF_FC) = &
                KKA(1:6*NDOF_FC, 1+(GFCID2-1)*6*NDOF_FC:GFCID2*6*NDOF_FC) - D2A_TMP
        END DO
    ELSE
        !boundary
        IF (BC_TYP(BCID).EQ.3) THEN
            PAIR = FACES_TAG(I,4)
            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF
            IF (FACES_TAG(PAIR,5).EQ.0) THEN
                TRID2 = FACES_TAG(PAIR,6)
                LFCID2 = FACES_TAG(PAIR,8)
            ELSE
                TRID2 = FACES_TAG(PAIR,5)
                LFCID2 = FACES_TAG(PAIR,7)
            END IF

            KKA(1:6*NDOF_FC,1+(I-1)*6*NDOF_FC:I*6*NDOF_FC) = AA_TMP
            DO IL = 1,3
                !interior
                GFCID1 = TRIANGLES_TAG(TRID1,3+IL)
                GFCID2 = TRIANGLES_TAG(TRID2,3+IL)

                D1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_TRACE(:,:,IL,TRID1))
                KKA(1:6*NDOF_FC, 1+(GFCID1-1)*6*NDOF_FC:GFCID1*6*NDOF_FC) = &
                    KKA(1:6*NDOF_FC, 1+(GFCID1-1)*6*NDOF_FC:GFCID1*6*NDOF_FC) - D1A_TMP

                D2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_TRACE(:,:,IL,TRID2))
                KKA(1:6*NDOF_FC, 1+(GFCID2-1)*6*NDOF_FC:GFCID2*6*NDOF_FC) = &
                    KKA(1:6*NDOF_FC, 1+(GFCID2-1)*6*NDOF_FC:GFCID2*6*NDOF_FC) - D2A_TMP
            END DO

        ELSE

            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF

            KKA(1:6*NDOF_FC,1+(I-1)*6*NDOF_FC:I*6*NDOF_FC) = AA_TMP
            DO IL = 1,3
                !interior
                GFCID1 = TRIANGLES_TAG(TRID1,3+IL)
                D1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_TRACE(:,:,IL,TRID1))
                KKA(1:6*NDOF_FC, 1+(GFCID1-1)*6*NDOF_FC:GFCID1*6*NDOF_FC) = &
                    KKA(1:6*NDOF_FC, 1+(GFCID1-1)*6*NDOF_FC:GFCID1*6*NDOF_FC) - D1A_TMP
            END DO
        END IF
    END IF

    !compression
    DO M = 1, 6*NDOF_FC
        IKKA(M+(I-1)*6*NDOF_FC) = NNZ !row: M+(I-1)*6*NDOF_FC
        DO L = 1, N_FCS*6*NDOF_FC
            IF (ABS(KKA(M,L)).GT.1.d-30) THEN
                NNZ = NNZ + 1
                IF (NNZ.GT.(5*N_FCS*6*NDOF_FC*6*NDOF_FC)) THEN
                    WRITE(*,*)'Range Exceed'
                    STOP
                END IF
                KKAcomp_tmp(NNZ) = KKA(M,L)
                JKKA_tmp(NNZ) = L
            END IF
        END DO
    END DO
END DO

IKKA(N_FCS*6*NDOF_FC+1) = NNZ
IKKA = IKKA + 1
ALLOCATE(KKAcomp(NNZ),JKKA(NNZ))
KKAcomp = KKAcomp_tmp(1:NNZ)
JKKA = JKKA_tmp(1:NNZ)


DEALLOCATE(KKA,AA,AA_TMP,D1A_TMP,D2A_TMP,KKAcomp_tmp,JKKA_tmp)

END SUBROUTINE

SUBROUTINE Calculate_SRC_ACC_HoTfromDVM ()
IMPLICIT NONE
REAL (KIND=DBL) :: eHx, eHy, pHx, pHy, eVDF_tmp, pVDF_tmp, eT_tmp, pT_tmp
INTEGER :: I, M, L, J1, J2

AA_SRC = 0.d0

!$OMP PARALLEL DO PRIVATE(I,M,eHx,eHy,pHx,pHy,eVDF_tmp,pVDF_tmp,eT_tmp,pT_tmp,J1,J2,L)
DO I = 1, N_TRIS
    DO M = 1, NDOF_TRI
        eHx = 0.d0
        eHy = 0.d0
        pHx = 0.d0
        pHy = 0.d0

        !*****************************************************************
        DO J2 = 1, NAZIM
        DO J1 = 1, NPOLE
            eVDF_tmp = 0.d0
            pVDF_tmp = 0.d0
            DO L = 1, NDOF_TRI
                !partial x of VDF
                eVDF_tmp = eVDF_tmp + eVDF(L,I,J1,J2)*INT_NODFUNC_TRI_TRI_X(I,L,M)
                pVDF_tmp = pVDF_tmp + pVDF(L,I,J1,J2)*INT_NODFUNC_TRI_TRI_X(I,L,M)
            END DO

            eHx = eHx + CX(J1,J2)*CX(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
            eHy = eHy + CY(J1,J2)*CX(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
            pHx = pHx + CX(J1,J2)*CX(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)
            pHy = pHy + CY(J1,J2)*CX(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)

            eVDF_tmp = 0.d0
            pVDF_tmp = 0.d0
            DO L = 1, NDOF_TRI
                !partial y of VDF
                eVDF_tmp = eVDF_tmp + eVDF(L,I,J1,J2)*INT_NODFUNC_TRI_TRI_Y(I,L,M)
                pVDF_tmp = pVDF_tmp + pVDF(L,I,J1,J2)*INT_NODFUNC_TRI_TRI_Y(I,L,M)
            END DO

            eHx = eHx + CX(J1,J2)*CY(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
            eHy = eHy + CY(J1,J2)*CY(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
            pHx = pHx + CX(J1,J2)*CY(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)
            pHy = pHy + CY(J1,J2)*CY(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)

        END DO
        END DO

        eT_tmp = 0.d0
        pT_tmp = 0.d0
        DO L = 1, NDOF_TRI
            !partial x of T
            eT_tmp = eT_tmp + eT_s(L,I)*INT_NODFUNC_TRI_TRI_X(I,L,M)
            pT_tmp = pT_tmp + pT_s(L,I)*INT_NODFUNC_TRI_TRI_X(I,L,M)
        END DO
        eHx = (              eHx - Ce_s/3.d0*eT_tmp                                 )*Ve_s**2.d0
        pHx = (        pHx - Cp_s/3.d0*(Kn_p_C/Kn_p_e*eT_tmp + Kn_p_C/Kn_p_p*pT_tmp))*Vp_s**2.d0

        eT_tmp = 0.d0
        pT_tmp = 0.d0
        DO L = 1, NDOF_TRI
            !partial y of T
            eT_tmp = eT_tmp + eT_s(L,I)*INT_NODFUNC_TRI_TRI_Y(I,L,M)
            pT_tmp = pT_tmp + pT_s(L,I)*INT_NODFUNC_TRI_TRI_Y(I,L,M)
        END DO
        eHy = (              eHy - Ce_s/3.d0*eT_tmp                                 )*Ve_s**2.d0
        pHy = (        pHy - Cp_s/3.d0*(Kn_p_C/Kn_p_e*eT_tmp + Kn_p_C/Kn_p_p*pT_tmp))*Vp_s**2.d0

        AA_SRC(M           ,I) = 0.d0
        AA_SRC(M+  NDOF_TRI,I) = 0.d0
        AA_SRC(M+2*NDOF_TRI,I) = -eHx*Kn_e_p
        AA_SRC(M+3*NDOF_TRI,I) = -eHy*Kn_e_p
        AA_SRC(M+4*NDOF_TRI,I) = -pHx*Kn_p_C 
        AA_SRC(M+5*NDOF_TRI,I) = -pHy*Kn_p_C 
    END DO
    AA_SRC(:,I) = MATMUL(inv_AA_SOL(:,:,I),AA_SRC(:,I))
END DO
!$OMP END PARALLEL DO

END SUBROUTINE


SUBROUTINE Global_Problem_Solver_ACC ()
IMPLICIT NONE
EXTERNAL :: PARDISO
INTEGER :: I,TRID1,TRID2,LFCID1,LFCID2,M,L,J1,J2,BCID,PAIR
REAL (KIND=DBL) :: nx, ny
REAL (KIND=DBL) :: eEbc, eqxbc, eqybc, eVDF_tmp
REAL (KIND=DBL) :: pEbc, pqxbc, pqybc, pVDF_tmp


U_TRACE = 0.d0
FFA = 0.d0

!$OMP PARALLEL DO PRIVATE(I,BCID,TRID1,TRID2,LFCID1,LFCID2,C1A_TMP,C2A_TMP,nx,ny,M,eqxbc,eqybc,eEbc, &
!$OMP L,eVDF_tmp,pEbc,pqxbc,pqybc,pVDF_tmp,J1,J2,PAIR)
DO I = 1, N_FCS 
    BCID = FACES_TAG(I,3)
    IF (BCID.EQ.0) THEN
        !interior
        TRID1 = FACES_TAG(I,5)
        TRID2 = FACES_TAG(I,6)
        LFCID1 = FACES_TAG(I,7)
        LFCID2 = FACES_TAG(I,8)

        C1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_SRC(:,TRID1))
        C2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_SRC(:,TRID2))

        FFA(1+(I-1)*6*NDOF_FC:I*6*NDOF_FC) = C1A_TMP + C2A_TMP
    ELSE
        !periodic boundary
        IF (BC_TYP(BCID).EQ.3) THEN
            PAIR = FACES_TAG(I,4)
            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF
            IF (FACES_TAG(PAIR,5).EQ.0) THEN
                TRID2 = FACES_TAG(PAIR,6)
                LFCID2 = FACES_TAG(PAIR,8)
            ELSE
                TRID2 = FACES_TAG(PAIR,5)
                LFCID2 = FACES_TAG(PAIR,7)
            END IF
            C1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_SRC(:,TRID1))
            C2A_TMP = MATMUL(BA_SOL(:,:,LFCID2,TRID2),AA_SRC(:,TRID2))

            FFA(1+(I-1)*6*NDOF_FC:I*6*NDOF_FC) = C1A_TMP + C2A_TMP
        ELSE

            IF (FACES_TAG(I,5).EQ.0) THEN
                TRID1 = FACES_TAG(I,6)
                LFCID1 = FACES_TAG(I,8)
            ELSE
                TRID1 = FACES_TAG(I,5)
                LFCID1 = FACES_TAG(I,7)
            END IF
            nx = TRIANGLES_INF(TRID1,2*LFCID1)
            ny = TRIANGLES_INF(TRID1,2*LFCID1+1)

            C1A_TMP = MATMUL(BA_SOL(:,:,LFCID1,TRID1),AA_SRC(:,TRID1))
            FFA(1+(I-1)*6*NDOF_FC:I*6*NDOF_FC) = C1A_TMP

            !adiabatic boundary
            ! IF (BC_TYP(BCID).EQ.2) THEN
            !     DO M = 1, NDOF_FC
            !         eqxbc = 0.d0
            !         eqybc = 0.d0
            !         eEbc = 0.d0

            !         DO J2 = 1, NAZIM
            !         DO J1 = 1, NPOLE
            !             eVDF_tmp = 0.d0
            !             DO L = 1, NDOF_TRI
            !                 eVDF_tmp = eVDF_tmp + eVDF(L,TRID1,J1,J2)*INT_NODFUNC_TRI_FC_FC(L,M,TRID1,LFCID1)
            !             END DO
            !             eEbc  = eEbc  +                eVDF_tmp*DOMEGA(J1,J2)
            !             eqxbc = eqxbc + Ve_s*CX(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
            !             eqybc = eqybc + Ve_s*CY(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
            !         END DO
            !         END DO
            !         ! FFA(M+          (I-1)*6*NDOF_FC) = eEbc
            !         ! FFA(M+2*NDOF_FC+(I-1)*6*NDOF_FC) = eqxbc
            !         ! FFA(M+3*NDOF_FC+(I-1)*6*NDOF_FC) = eqybc
            !     END DO
            ! END IF
            ! IF (BC_TYP_P(BCID).EQ.2) THEN
            !     DO M = 1, NDOF_FC
            !         pqxbc = 0.d0
            !         pqybc = 0.d0
            !         pEbc = 0.d0

            !         DO J2 = 1, NAZIM
            !         DO J1 = 1, NPOLE
            !             pVDF_tmp = 0.d0
            !             DO L = 1, NDOF_TRI
            !                 pVDF_tmp = pVDF_tmp + pVDF(L,TRID1,J1,J2)*INT_NODFUNC_TRI_FC_FC(L,M,TRID1,LFCID1)
            !             END DO
            !             pEbc  = pEbc  +                pVDF_tmp*DOMEGA(J1,J2)
            !             ! pqxbc = pqxbc + Vp_s*CX(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)*1.d0/A2 - A1/A2*eqxbc
            !             ! pqybc = pqybc + Vp_s*CY(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)*1.d0/A2 - A1/A2*eqybc
            !             pqxbc = pqxbc + Vp_s*CX(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)
            !             pqybc = pqybc + Vp_s*CY(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)
            !         END DO
            !         END DO
            !         ! FFA(M+  NDOF_FC+(I-1)*6*NDOF_FC) = pEbc 
            !         ! FFA(M+4*NDOF_FC+(I-1)*6*NDOF_FC) = pqxbc
            !         ! FFA(M+5*NDOF_FC+(I-1)*6*NDOF_FC) = pqybc
            !     END DO
            ! END IF

            !thermarlising wall
            IF (BC_TYP(BCID).EQ.1) THEN
                DO M = 1, NDOF_FC
                    eqxbc = 0.d0
                    eqybc = 0.d0
                    eEbc = 0.d0

                    DO J2 = 1, NAZIM
                    DO J1 = 1, NPOLE
                        eVDF_tmp = 0.d0
                        DO L = 1, NDOF_TRI
                            eVDF_tmp = eVDF_tmp + eVDF(L,TRID1,J1,J2)*INT_NODFUNC_TRI_FC_FC(L,M,TRID1,LFCID1)
                        END DO
                        eEbc  = eEbc   +                eVDF_tmp*DOMEGA(J1,J2)
                        eqxbc = eqxbc  + Ve_s*CX(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
                        eqybc = eqybc  + Ve_s*CY(J1,J2)*eVDF_tmp*DOMEGA(J1,J2)
                    END DO
                    END DO
                    FFA(M+          (I-1)*6*NDOF_FC) = eEbc
                END DO
            END IF
            IF (BC_TYP_P(BCID).EQ.1) THEN
                DO M = 1, NDOF_FC
                    pqxbc = 0.d0
                    pqybc = 0.d0
                    pEbc = 0.d0

                    DO J2 = 1, NAZIM
                    DO J1 = 1, NPOLE
                        pVDF_tmp = 0.d0
                        DO L = 1, NDOF_TRI
                            pVDF_tmp = pVDF_tmp + pVDF(L,TRID1,J1,J2)*INT_NODFUNC_TRI_FC_FC(L,M,TRID1,LFCID1)
                        END DO
                        pEbc  = pEbc   +                pVDF_tmp*DOMEGA(J1,J2)
                        pqxbc = pqxbc + Vp_s*CX(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)
                        pqybc = pqybc + Vp_s*CY(J1,J2)*pVDF_tmp*DOMEGA(J1,J2)
                    END DO
                    END DO
                    FFA(M+  NDOF_FC+(I-1)*6*NDOF_FC) = pEbc 
                END DO
            END IF
        END IF
    END IF
END DO
!$OMP END PARALLEL DO

!******************************************************************************************************
M = N_FCS*6*NDOF_FC
PHASE = 33
CALL PARDISO(PT,MAXFCT,1,MTYPE,PHASE,M,KKAcomp,IKKA,JKKA,PERM,1,IPARM,MSGLVL,FFA,HHA,ERROR)
IF (ERROR.NE.0) WRITE(*,*) 'PARDISO solution error', ERROR

!$OMP PARALLEL DO PRIVATE(I)
DO I = 1, N_FCS
    U_TRACE(:,I) = HHA(1+(I-1)*6*NDOF_FC:I*6*NDOF_FC)
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE Local_Problem_Solver_ACC ()
IMPLICIT NONE
INTEGER :: I, IL

!$OMP PARALLEL DO PRIVATE (I,IL)
DO I = 1, N_TRIS
    UQ(:,I) = AA_SRC(:,I)
    DO IL = 1,3
        UQ(:,I) = UQ(:,I) + MATMUL(AA_TRACE(:,:,IL,I),U_TRACE(:,TRIANGLES_TAG(I,3+IL)))
    END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE Init_PARDISO ()
IMPLICIT NONE
EXTERNAL :: pardisoinit, PARDISO
INTEGER :: N

WRITE(*,*) 'Initialize PARDISO ...'
MAXFCT = 1
MTYPE = 11
MSGLVL = 0

N = N_FCS*6*NDOF_FC
ALLOCATE(PERM(N))
PERM = 0

CALL pardisoinit(PT,MTYPE,IPARM)

IPARM(27) = 1
PHASE = 12
CALL PARDISO(PT,MAXFCT,1,MTYPE,PHASE,N,KKAcomp,IKKA,JKKA,PERM,1,IPARM,MSGLVL,FFA,HHA,ERROR)
IF (ERROR.NE.0) WRITE(*,*) 'PARDISO analysis & factorization error', ERROR
IF (ERROR.EQ.0) WRITE(*,*) 'PARDISO LU factorization successful'

END SUBROUTINE

SUBROUTINE Correct_VDF_Calculate_Macro_Properties ()
IMPLICIT NONE
INTEGER :: I, M, J1, J2
REAL (KIND=DBL) :: beta, Kn_loc
REAL (KIND=DBL) :: eT_ACC, eE_ACC, eqx_ACC, eqy_ACC
REAL (KIND=DBL) :: eT_VDF, eE_VDF, eqx_VDF, eqy_VDF
REAL (KIND=DBL) :: pT_ACC, pE_ACC, pqx_ACC, pqy_ACC
REAL (KIND=DBL) :: pT_VDF, pE_VDF, pqx_VDF, pqy_VDF
REAL (KIND=DBL) :: l_eT, l_eE, l_eqx, l_eqy
REAL (KIND=DBL) :: l_pT, l_pE, l_pqx, l_pqy

eTemp = 0.d0
eEnergy = 0.d0
eQx = 0.d0
eQy = 0.d0

pTemp = 0.d0
pEnergy = 0.d0
pQx = 0.d0
pQy = 0.d0

!$OMP PARALLEL DO PRIVATE(I,M,J1,J2,beta,Kn_loc, &
!$OMP eT_ACC, eE_ACC, eqx_ACC, eqy_ACC, eT_VDF, eqx_VDF, eqy_VDF,l_eT, l_eE, l_eqx, l_eqy, &
!$OMP pT_ACC, pE_ACC, pqx_ACC, pqy_ACC, pT_VDF, pqx_VDF, pqy_VDF, l_pT, l_pE, l_pqx, l_pqy)
DO I = 1, N_TRIS
        Kn_loc = Kn_e_p/TRIANGLES_Hmin(I)
        beta = 1.d0 !MIN(Kn_loc,KN_THR)/Kn_loc

    DO M = 1, NDOF_TRI
        eT_VDF  = eT_s(M,I)
        eE_VDF  = eE_s(M,I)
        eqx_VDF = eQx_s(M,I)
        eqy_VDF = eQy_s(M,I)
        pT_VDF  = pT_s(M,I)
        pE_VDF  = pE_s(M,I)
        pqx_VDF = pQx_s(M,I)
        pqy_VDF = pQy_s(M,I)

        eE_ACC  = UQ(           M, I)
        pE_ACC  = UQ(  NDOF_TRI+M, I)
        eqx_ACC = UQ(2*NDOF_TRI+M, I)
        eqy_ACC = UQ(3*NDOF_TRI+M, I)
        pqx_ACC = UQ(4*NDOF_TRI+M, I)
        pqy_ACC = UQ(5*NDOF_TRI+M, I)

        eT_ACC =weT*eE_ACC + wpT*pE_ACC
        pT_ACC = pE_ACC/Cp_s


        l_eT  = ( eT_ACC -  eT_VDF)*beta
        l_eE  = ( eE_ACC -  eE_VDF)*beta
        l_eqx = (eqx_ACC - eqx_VDF)*beta
        l_eqy = (eqy_ACC - eqy_VDF)*beta
        l_pT  = ( pT_ACC -  pT_VDF)*beta
        l_pE  = ( pE_ACC -  pE_VDF)*beta
        l_pqx = (pqx_ACC - pqx_VDF)*beta
        l_pqy = (pqy_ACC - pqy_VDF)*beta

        DO J2 = 1, NAZIM
        DO J1 = 1, NPOLE
            eVDF(M,I,J1,J2) = eVDF(M,I,J1,J2) + l_eE/4.d0/PI 
            pVDF(M,I,J1,J2) = pVDF(M,I,J1,J2) + l_pE/4.d0/PI
        END DO
        END DO

        eT_s(M,I)  = l_eT  + eT_VDF
        eE_s(M,I)  = l_eE  + eE_VDF
        eQx_s(M,I) = l_eqx + eqx_VDF
        eQy_s(M,I) = l_eqy + eqy_VDF
        pT_s(M,I)  = l_pT  + pT_VDF
        pE_s(M,I)  = l_pE  + pE_VDF
        pQx_s(M,I) = l_pqx + pqx_VDF
        pQy_s(M,I) = l_pqy + pqy_VDF


        eTemp(I)   =   eTemp(I) +  eT_s(M,I)*INT_NODFUNC_TRI(M,I)
        eEnergy(I) = eEnergy(I) +  eE_s(M,I)*INT_NODFUNC_TRI(M,I)
        eQx(I)     =     eQx(I) + eQx_s(M,I)*INT_NODFUNC_TRI(M,I)
        eQy(I)     =     eQy(I) + eQy_s(M,I)*INT_NODFUNC_TRI(M,I)
        pTemp(I)   =   pTemp(I) +  pT_s(M,I)*INT_NODFUNC_TRI(M,I)
        pEnergy(I) = pEnergy(I) +  pE_s(M,I)*INT_NODFUNC_TRI(M,I)
        pQx(I)     =     pQx(I) + pQx_s(M,I)*INT_NODFUNC_TRI(M,I)
        pQy(I)     =     pQy(I) + pQy_s(M,I)*INT_NODFUNC_TRI(M,I)
    END DO
END DO
!$OMP END PARALLEL DO

eMASS = SUM(eEnergy)
pMASS= SUM(pEnergy)

END SUBROUTINE

SUBROUTINE Release_PARDISO ()
IMPLICIT NONE
EXTERNAL :: PARDISO
INTEGER :: N

WRITE(*,*) 'Release PARDISO ...'
N = N_FCS*6*NDOF_FC
PHASE = -1
CALL PARDISO(PT,MAXFCT,1,MTYPE,PHASE,N,KKAcomp,IKKA,JKKA,PERM,1,IPARM,MSGLVL,FFA,HHA,ERROR)

END SUBROUTINE

END MODULE
