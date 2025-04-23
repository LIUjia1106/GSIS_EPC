MODULE VELOCITY_DISTRIBUTION_FUNCTION
USE CONSTANT
USE GLOBAL_VARIABLE
USE SPATIAL_GRID
USE NODAL_FUNCTION
USE VELOCITY_GRID
USE BOUNDARY_CONDITION
USE INTEGRATIONS
IMPLICIT NONE

REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: eVDF, pVDF  !perturbed eVDF, pVDF, DIMENSION (NDOF_TRI,N_TRIS,NPOLE,NAZIM)
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: eTemp, eEnergy, eQx, eQy !cell average value
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: pTemp, pEnergy, pQx, pQy
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: eT_OLD, eQx_OLD, eQy_OLD, eE_OLD
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:), SAVE :: pT_OLD, pQx_OLD, pQy_OLD, pE_OLD
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: eT_s, eE_s, eQx_s, eQy_s
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: pT_s, pE_s, pQx_s, pQy_s
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:), SAVE :: eFLUX_WALL, pFLUX_WALL

REAL (KIND=DBL), SAVE :: eMASS, pMASS

CONTAINS

SUBROUTINE Init_Velocity_Distribution_Function()
IMPLICIT NONE
INTEGER :: ios
CHARACTER(LEN=100) :: CHAR1

ALLOCATE(eVDF(NDOF_TRI,N_TRIS,NPOLE,NAZIM),pVDF(NDOF_TRI,N_TRIS,NPOLE,NAZIM))
ALLOCATE(eTemp(N_TRIS),eEnergy(N_TRIS),eQx(N_TRIS),eQy(N_TRIS))
ALLOCATE(pTemp(N_TRIS),pEnergy(N_TRIS),pQx(N_TRIS),pQy(N_TRIS))
ALLOCATE(eT_OLD(N_TRIS),eQx_OLD(N_TRIS),eQy_OLD(N_TRIS),eE_OLD(N_TRIS))
ALLOCATE(pT_OLD(N_TRIS),pQx_OLD(N_TRIS),pQy_OLD(N_TRIS),pE_OLD(N_TRIS))
ALLOCATE(eT_s(NDOF_TRI,N_TRIS),eE_s(NDOF_TRI,N_TRIS),eQx_s(NDOF_TRI,N_TRIS),eQy_s(NDOF_TRI,N_TRIS))
ALLOCATE(pT_s(NDOF_TRI,N_TRIS),pE_s(NDOF_TRI,N_TRIS),pQx_s(NDOF_TRI,N_TRIS),pQy_s(NDOF_TRI,N_TRIS))
ALLOCATE(eFLUX_WALL(NDOF_TRI,N_FCS_B),pFLUX_WALL(NDOF_TRI,N_FCS_B))

WRITE(*,*)'**** Initial Velocity Distribution Function ****'
! initialize electron-related quantities
eVDF = 0.d0
eTemp = 0.d0
eEnergy = 0.d0
eQx = 0.d0
eQy = 0.d0

eT_OLD  = 0.d0
eQx_OLD = 0.d0
eQy_OLD = 0.d0
eE_OLD = 0.d0

eT_s = 0.d0
eE_s = 0.d0
eQx_s = 0.d0
eQy_s = 0.d0

! initialize ponon-related quantities
pVDF = 0.d0
pTemp = 0.d0
pEnergy = 0.d0
pQx = 0.d0
pQy = 0.d0

pT_OLD  = 0.d0
pQx_OLD = 0.d0
pQy_OLD = 0.d0
pE_OLD = 0.d0

pT_s = 0.d0
pE_s = 0.d0
pQx_s = 0.d0
pQy_s = 0.d0

WRITE(CHAR1,108)'P',DEG,'T',N_TRIS,'NP',NPOLE,'NA',NAZIM
108 FORMAT(A1,I1,A1,I3,A2,I2,A2,I2)
OPEN (UNIT = 12, FILE='eVDF'//TRIM(CHAR1)//'.out', STATUS='OLD', ACCESS='STREAM', FORM='UNFORMATTED',IOSTAT=ios)
IF (ios.EQ.0) THEN
WRITE(*,*)'read from file: ', 'eVDF'//TRIM(CHAR1)//'.out'
READ (12) eVDF
END IF
CLOSE(12)
OPEN (UNIT = 13, FILE='pVDF'//TRIM(CHAR1)//'.out', STATUS='OLD', ACCESS='STREAM', FORM='UNFORMATTED',IOSTAT=ios)
IF (ios.EQ.0) THEN
WRITE(*,*)'read from file: ', 'pVDF'//TRIM(CHAR1)//'.out'
READ (13) pVDF
END IF
CLOSE(13)

END SUBROUTINE Init_Velocity_Distribution_Function

!************************************************************
SUBROUTINE Calculate_FLUX_WALL ()
IMPLICIT NONE
INTEGER :: I,J1,J2,M,L, BCID, TRID, LFCID
REAL (KIND=DBL) :: ss1, ss2, ss3, ss4, n1, n2, u, A1

eFLUX_WALL =0.d0
pFLUX_WALL =0.d0

A1 = Vp_s**2.d0*Cp_s*Kn_p_C**2.d0/(Ve_s**2.d0*Ce_s*Kn_e_p*Kn_p_e)
!$OMP PARALLEL DO PRIVATE (I,J1,J2,M,L,BCID,TRID,LFCID,ss1,ss2,ss3,ss4,n1,n2,u)
DO I = 1, N_FCS_B
    BCID = FACES_TAG(I,3)
    IF(BC_TYP(BCID).EQ.2)THEN
        !adiabatic for electron
        IF (FACES_TAG(I,5).EQ.0) THEN
            TRID = FACES_TAG(I,6)
            LFCID = FACES_TAG(I,8)
        ELSE
            TRID = FACES_TAG(I,5)
            LFCID = FACES_TAG(I,7)
        END IF
        n1 = -TRIANGLES_INF(TRID,2*LFCID)
        n2 = -TRIANGLES_INF(TRID,2*LFCID+1)

        DO M = 1, NDOF_TRI
            ss1 = 0.d0
            ss2 = 0.d0
            DO J2 = 1, NAZIM
            DO J1 = 1, NPOLE
                u = CX(J1,J2)*n1+CY(J1,J2)*n2
                IF (u.LT.0.d0)THEN
                    DO L = 1, NDOF_TRI
                        ss1 = ss1 + u*Ve_s*DOMEGA(J1,J2)*eVDF(L,TRID,J1,J2)*INT_NODFUNC_TRI_FC(TRID,LFCID,M,L)
                    END DO
                ELSE
                    ss2 = ss2 - u*Ve_s*DOMEGA(J1,J2)
                END IF
            END DO
            END DO
            eFLUX_WALL(M,I) = -ss1/ss2
        END DO
    END IF
    IF(BC_TYP_P(BCID).EQ.2)THEN
        !adiabatic for phonon
        IF (FACES_TAG(I,5).EQ.0) THEN
            TRID = FACES_TAG(I,6)
            LFCID = FACES_TAG(I,8)
        ELSE
            TRID = FACES_TAG(I,5)
            LFCID = FACES_TAG(I,7)
        END IF
        n1 = -TRIANGLES_INF(TRID,2*LFCID)
        n2 = -TRIANGLES_INF(TRID,2*LFCID+1)

        DO M = 1, NDOF_TRI
            ss1 = 0.d0
            ss2 = 0.d0
            ss3 = 0.d0
            ss4 = 0.d0
            DO J2 = 1, NAZIM
            DO J1 = 1, NPOLE
                u = CX(J1,J2)*n1+CY(J1,J2)*n2
                IF (u.LT.0.d0)THEN
                    DO L = 1, NDOF_TRI
                        ss3 = ss3 + u*Vp_s*DOMEGA(J1,J2)*pVDF(L,TRID,J1,J2)*INT_NODFUNC_TRI_FC(TRID,LFCID,M,L)
                    END DO
                ELSE
                    ss4 = ss4 - u*Vp_s*DOMEGA(J1,J2)
                END IF
            END DO
            END DO
            pFLUX_WALL(M,I) = -ss3/ss4
        END DO
    END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE Calculate_Macro_Properties ()
IMPLICIT NONE
INTEGER :: I,L
REAL (KIND=DBL) :: ss,ss_x,ss_y,ssp,ssp_x,ssp_y

eTemp = 0.d0
eEnergy = 0.d0
eQx = 0.d0
eQy = 0.d0

eT_s = 0.d0
eE_s = 0.d0
eQx_s = 0.d0
eQy_s = 0.d0

pTemp = 0.d0
pEnergy = 0.d0
pQx = 0.d0
pQy = 0.d0

pT_s = 0.d0
pE_s = 0.d0
pQx_s = 0.d0
pQy_s = 0.d0

!$OMP PARALLEL DO PRIVATE (I,L,ss,ss_x,ss_y,ssp,ssp_x,ssp_y)
DO I = 1, N_TRIS
    DO L = 1, NDOF_TRI
        ss   = SUM(        eVDF(L,I,:,:)*DOMEGA)
        ss_x = SUM(Ve_s*CX*eVDF(L,I,:,:)*DOMEGA)
        ss_y = SUM(Ve_s*CY*eVDF(L,I,:,:)*DOMEGA)

        ssp   = SUM(        pVDF(L,I,:,:)*DOMEGA)
        ssp_x = SUM(Vp_s*CX*pVDF(L,I,:,:)*DOMEGA)
        ssp_y = SUM(Vp_s*CY*pVDF(L,I,:,:)*DOMEGA)

        eT_s(L,I)  = (Kn_p_e*ss+Kn_e_p*ssp)/(Kn_p_e*Ce_s+Kn_e_p*Cp_s)
        eE_s(L,I)  = ss
        eQx_s(L,I) = ss_x
        eQy_s(L,I) = ss_y

        pT_s(L,I)  = ssp/Cp_s
        pE_s(L,I)  = ssp
        pQx_s(L,I) = ssp_x
        pQy_s(L,I) = ssp_y

        eTemp(I)   = eTemp(I)   + eT_s(L,I)*INT_NODFUNC_TRI(L,I)  
        eEnergy(I) = eEnergy(I) + ss       *INT_NODFUNC_TRI(L,I)
        eQx(I)     = eQx(I)     + ss_x     *INT_NODFUNC_TRI(L,I)
        eQy(I)     = eQy(I)     + ss_y     *INT_NODFUNC_TRI(L,I)

        pTemp(I)   = pTemp(I)   + pT_s(L,I) *INT_NODFUNC_TRI(L,I)
        pEnergy(I) = pEnergy(I) + ssp       *INT_NODFUNC_TRI(L,I)
        pQx(I)     = pQx(I)     + pQx_s(L,I)*INT_NODFUNC_TRI(L,I)
        pQy(I)     = pQy(I)     + pQy_s(L,I)*INT_NODFUNC_TRI(L,I)
    END DO
END DO
!$OMP END PARALLEL DO

eMASS = SUM(eEnergy)
pMASS = SUM(pEnergy)

END SUBROUTINE Calculate_Macro_Properties

SUBROUTINE Calculate_Residual_T (RESIDUALe,RESIDUALp)
IMPLICIT NONE
REAL (KIND=DBL), INTENT(OUT) :: RESIDUALe,RESIDUALp
REAL (KIND=DBL) :: ss1,ss2,ss3,ss4
INTEGER :: I

ss1 = 0.d0
ss2 = 0.d0
ss3 = 0.d0
ss4 = 0.d0
!$OMP PARALLEL DO PRIVATE (I) REDUCTION(+:ss1,ss2,ss3,ss4)
DO I = 1, N_TRIS
    ss1 = ss1 + (eTemp(I)-eT_OLD(I))*(eTemp(I)-eT_OLD(I))
    ss2 = ss2 + eTemp(I)*eTemp(I)
    ss3 = ss3 + (pTemp(I)-pT_OLD(I))*(pTemp(I)-pT_OLD(I))
    ss4 = ss4 + pTemp(I)*pTemp(I)
END DO
!$OMP END PARALLEL DO

RESIDUALe = sqrt(ss1/ss2)
RESIDUALp = sqrt(ss3/ss4)

eT_OLD = eTemp
pT_OLD = pTemp

END SUBROUTINE Calculate_Residual_T

SUBROUTINE Calculate_Residual_ALL (RESIDUALe,RESIDUALp)
IMPLICIT NONE
REAL (KIND=DBL), INTENT(OUT) :: RESIDUALe,RESIDUALp
REAL (KIND=DBL), DIMENSION(4) :: ss1,ss2,ss3,ss4
INTEGER :: I
REAL (KIND=DBL) :: q, q_old

ss1 = 0.d0
ss2 = 0.d0
!$OMP PARALLEL DO PRIVATE (I,q,q_old) REDUCTION(+:ss1,ss2,ss3,ss4)
DO I = 1, N_TRIS
    q = SQRT(eQx(I)*eQx(I)+eQy(I)*eQy(I))
    q_old = SQRT(eQx_OLD(I)*eQx_OLD(I)+eQy_OLD(I)*eQy_OLD(I))

    ss1(1) = ss1(1) + (eTemp(I)-eT_OLD(I))*(eTemp(I)-eT_OLD(I))
    ss1(2) = ss1(2) + (eQx(I)-eQx_OLD(I))*(eQx(I)-eQx_OLD(I))
    ss1(3) = ss1(3) + (eQy(I)-eQy_OLD(I))*(eQy(I)-eQy_OLD(I))
    ss1(4) = ss1(4) + (eEnergy(I)-eE_OLD(I))*(eEnergy(I)-eE_OLD(I))

    ss2(1) = ss2(1) + eTemp(I)*eTemp(I)
    ss2(2) = ss2(2) + eQx(I)*eQx(I)
    ss2(3) = ss2(3) + eQy(I)*eQy(I)
    ss2(4) = ss2(4) + eEnergy(I)*eEnergy(I)

    ss3(1) = ss3(1) + (pTemp(I)-pT_OLD(I))*(pTemp(I)-pT_OLD(I))
    ss3(2) = ss3(2) + (pQx(I)-pQx_OLD(I))*(pQx(I)-pQx_OLD(I))
    ss3(3) = ss3(3) + (pQy(I)-pQy_OLD(I))*(pQy(I)-pQy_OLD(I))
    ss3(4) = ss3(4) + (pEnergy(I)-pE_OLD(I))*(pEnergy(I)-pE_OLD(I))

    ss4(1) = ss4(1) + pTemp(I)*pTemp(I)
    ss4(2) = ss4(2) + pQx(I)*pQx(I)
    ss4(3) = ss4(3) + pQy(I)*pQy(I)
    ss4(4) = ss4(4) + pEnergy(I)*pEnergy(I)
END DO
!$OMP END PARALLEL DO

RESIDUALe = MAXVAL(SQRT(ss1/ss2))
RESIDUALp = MAXVAL(SQRT(ss3/ss4))

eT_OLD  = eTemp
eQx_OLD = eQx
eQy_OLD = eQy
eE_OLD  = eEnergy

pT_OLD  = pTemp
pQx_OLD = pQx
pQy_OLD = pQy
pE_OLD  = pEnergy

END SUBROUTINE Calculate_Residual_ALL

SUBROUTINE Compress_VDF()
IMPLICIT NONE
CHARACTER(LEN=100) :: CHAR1

WRITE(CHAR1,108)'P',DEG,'T',N_TRIS,'NP',NPOLE,'NA',NAZIM
108 FORMAT(A1,I1,A1,I3,A2,I2,A2,I2)
OPEN (UNIT = 12, FILE='eVDF'//TRIM(CHAR1)//'.out', STATUS='UNKNOWN', ACCESS='STREAM', FORM='UNFORMATTED')
WRITE(12) eVDF
CLOSE(12)
WRITE(*,*)'Compress eVDF'
OPEN (UNIT = 13, FILE='pVDF'//TRIM(CHAR1)//'.out', STATUS='UNKNOWN', ACCESS='STREAM', FORM='UNFORMATTED')
WRITE(13) pVDF
CLOSE(13)
WRITE(*,*)'Compress pVDF'

END SUBROUTINE Compress_VDF

END MODULE VELOCITY_DISTRIBUTION_FUNCTION
