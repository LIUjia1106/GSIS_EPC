MODULE OUT_PUT
USE CONSTANT
USE GLOBAL_VARIABLE
USE SPATIAL_GRID
USE NODAL_FUNCTION
USE VELOCITY_GRID
USE VELOCITY_DISTRIBUTION_FUNCTION
USE ACCELERATION
IMPLICIT NONE

INTEGER :: Nx = 109
INTEGER :: Ny = 109
INTEGER :: I, IL, J, K,M,L,KL,JL,J1,J2
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:) :: px,py,pm
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:) :: peT, peE, peqx, peqy
REAL (KIND=DBL), ALLOCATABLE, DIMENSION (:,:) :: ppT, ppE, ppqx, ppqy
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TRIID

CONTAINS


SUBROUTINE Out_Put_Conduction()
IMPLICIT NONE
REAL (KIND=DBL), DIMENSION(4) :: x,y
LOGICAL :: INTRI
REAL (KIND=DBL) :: AREA, A, B, C, D, E, F, G, H, xi, eta
REAL (KIND=DBL) ::  ss,ss_x,ss_y,ss_xx,ss_xy,ss_yy
REAL (KIND=DBL) ::  ssp,ssp_x,ssp_y,ssp_xx,ssp_xy,ssp_yy

CHARACTER(LEN=100) :: FILENAME,INTG3
CHARACTER(LEN=10) :: INTG1, INTG2, INTG4, INTG5, INTG6

WRITE(*,*)'Output Conduction'

ALLOCATE(px(Nx),py(Ny),pm(NDOF_TRI),TRIID(Nx,Ny))
ALLOCATE(peT(Nx,Ny),peE(Nx,Ny),peqx(Nx,Ny),peqy(Nx,Ny))
ALLOCATE(ppT(Nx,Ny),ppE(Nx,Ny),ppqx(Nx,Ny),ppqy(Nx,Ny))
! ALLOCATE(peNxx(Nx,Ny),peNxy(Nx,Ny),peNyy(Nx,Ny),ppNxx(Nx,Ny),ppNxy(Nx,Ny),ppNyy(Nx,Ny))

px = 0.d0
py = 0.d0
pm = 0.d0
peT = 0.d0
peE = 0.d0
peqx = 0.d0
peqy = 0.d0
ppT = 0.d0
ppE = 0.d0
ppqx = 0.d0
ppqy = 0.d0

! peNxx = 0.d0
! peNxy = 0.d0
! peNyy = 0.d0
! ppNxx = 0.d0
! ppNxy = 0.d0
! ppNyy = 0.d0
TRIID = 0



DO J = 1, 5
A = (J-1.d0)/8.d0
px(J) = A*A*A*(10.d0-15.d0*A+6.d0*A*A)*0.02d0
px(Nx-J+1) = 1.d0 - px(J)
END DO

DO J = 6, Nx-5
px(J) = (J-5.d0)*(1.d0-0.02d0)/100.d0 + 0.01d0
END DO
py = gamma*px

DO J = 1, Ny
    DO I = 1, Nx
        DO K = 1, N_TRIS
            x(1) = NODES(TRIANGLES_TAG(K,1),1)
            y(1) = NODES(TRIANGLES_TAG(K,1),2)
            x(2) = NODES(TRIANGLES_TAG(K,2),1)
            y(2) = NODES(TRIANGLES_TAG(K,2),2)
            x(3) = NODES(TRIANGLES_TAG(K,3),1)
            y(3) = NODES(TRIANGLES_TAG(K,3),2)
            x(4) = x(1)
            y(4) = y(1)

            INTRI = .TRUE.
            DO IL = 1,3
                AREA = 0.5*(px(I)*(y(IL)-y(IL+1))-py(J)*(x(IL)-x(IL+1))+(x(IL)*y(IL+1)-x(IL+1)*y(IL)))
                IF (AREA<-1.d-12) THEN
                    INTRI = .FALSE.
                    EXIT
                END IF
            END DO
            IF (INTRI) THEN
                TRIID(I,J) = K
                EXIT
            END IF
        END DO
    END DO
 END DO

DO J = 1, Ny
    DO I = 1, Nx
        IF (TRIID(I,J).EQ.0) THEN
            WRITE(*,*)'Warning: the point (',I,',',J,') is not in any triangle'
        ELSE
            x(1) = NODES(TRIANGLES_TAG(TRIID(I,J),1),1)
            y(1) = NODES(TRIANGLES_TAG(TRIID(I,J),1),2)
            x(2) = NODES(TRIANGLES_TAG(TRIID(I,J),2),1)
            y(2) = NODES(TRIANGLES_TAG(TRIID(I,J),2),2)
            x(3) = NODES(TRIANGLES_TAG(TRIID(I,J),3),1)
            y(3) = NODES(TRIANGLES_TAG(TRIID(I,J),3),2)
            x(4) = x(1)
            y(4) = y(1)

            A=(x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1))
            B=y(3)-y(1)
            C=x(1)-x(3)
            D=(x(3)-x(1))*y(1)-(y(3)-y(1))*x(1)

            E=(x(3)-x(1))*(y(2)-y(1))-(x(2)-x(1))*(y(3)-y(1))
            F=y(2)-y(1)
            G=x(1)-x(2)
            H=(x(2)-x(1))*y(1)-(y(2)-y(1))*x(1)

            xi = B/A*px(I)+C/A*py(J)+D/A
            eta = F/E*px(I)+G/E*py(J)+H/E

            pm = 0.d0
            DO M = 1, NDOF_TRI
                L = 1
                DO KL = 0,DEG
                    DO JL = 0, KL
                        IL = KL - JL
                        pm(M) = pm(M) + NODFUN_TRI_REF(M,L)*(xi**REAL(IL,DBL))*(eta**REAL(JL,DBL))
                        L = L + 1
                    END DO
                END DO
            END DO

            ss = 0.d0
            ss_x = 0.d0
            ss_y = 0.d0

            ssp = 0.d0
            ssp_x = 0.d0
            ssp_y = 0.d0

            DO J2 = 1, NAZIM
                DO J1 = 1, NPOLE
                    DO L = 1, NDOF_TRI
                        ss   = ss   +           DOMEGA(J1,J2)*eVDF(L,TRIID(I,J),J1,J2)*pm(L)
                        ss_x = ss_x + Ve_s*CX(J1,J2)*DOMEGA(J1,J2)*eVDF(L,TRIID(I,J),J1,J2)*pm(L)
                        ss_y = ss_y + Ve_s*CY(J1,J2)*DOMEGA(J1,J2)*eVDF(L,TRIID(I,J),J1,J2)*pm(L)

                        ssp   = ssp   +           DOMEGA(J1,J2)*pVDF(L,TRIID(I,J),J1,J2)*pm(L)
                        ssp_x = ssp_x + Vp_s*CX(J1,J2)*DOMEGA(J1,J2)*pVDF(L,TRIID(I,J),J1,J2)*pm(L)
                        ssp_y = ssp_y + Vp_s*CY(J1,J2)*DOMEGA(J1,J2)*pVDF(L,TRIID(I,J),J1,J2)*pm(L)
                    END DO
                END DO
            END DO

            peT(I,J)  = (Kn_p_e*ss+Kn_e_p*ssp)/(Kn_p_e*Ce_s+Kn_e_p*Cp_s)
            peE(I,J)  = ss
            peqx(I,J) = ss_x
            peqy(I,J) = ss_y

            ppT(I,J)  = ssp/Cp_s
            ppE(I,J)  = ssp
            ppqx(I,J) = ssp_x
            ppqy(I,J) = ssp_y

            ! ss_xx=0.d0
            ! ss_xy=0.d0
            ! ss_yy=0.d0

            ! ssp_xx=0.d0
            ! ssp_xy=0.d0
            ! ssp_yy=0.d0
            ! IF (ACCFLAG.EQ.1) THEN
            ! DO L =1, NDOF_TRI
            !     ss_xx = ss_xx +(-4.d0/3.d0*UQ(3*NDOF_TRI+L,TRIID(I,J))+2.d0/3.d0*UQ(6*NDOF_TRI+L,TRIID(I,J)))*pm(L)
            !     ss_xy = ss_xy +(-          UQ(4*NDOF_TRI+L,TRIID(I,J))-          UQ(5*NDOF_TRI+L,TRIID(I,J)))*pm(L)
            !     ss_yy = ss_yy +(+2.d0/3.d0*UQ(3*NDOF_TRI+L,TRIID(I,J))-4.d0/3.d0*UQ(6*NDOF_TRI+L,TRIID(I,J)))*pm(L)

            !     ssp_xx = ssp_xx +(-4.d0/3.d0*UQ(3*NDOF_TRI+L,TRIID(I,J))+2.d0/3.d0*UQ(6*NDOF_TRI+L,TRIID(I,J)))*pm(L)
            !     ssp_xy = ssp_xy +(-          UQ(4*NDOF_TRI+L,TRIID(I,J))-          UQ(5*NDOF_TRI+L,TRIID(I,J)))*pm(L)
            !     ssp_yy = ssp_yy +(+2.d0/3.d0*UQ(3*NDOF_TRI+L,TRIID(I,J))-4.d0/3.d0*UQ(6*NDOF_TRI+L,TRIID(I,J)))*pm(L)
            ! END DO
            ! END IF
            ! peNxx(I,J)=ss_xx
            ! peNxy(I,J)=ss_xy
            ! peNyy(I,J)=ss_yy
            ! ppNxx(I,J)=ssp_xx
            ! ppNxy(I,J)=ssp_xy
            ! ppNyy(I,J)=ssp_yy
        END IF
    END DO
END DO

! CALL Cal_Dimensional_Property ()

!**********************************************************
WRITE(INTG1,101)'P',DEG
101 FORMAT (A1,I1)
WRITE(INTG2,102)'T',N_TRIS
102 FORMAT (A1,I5)
WRITE(INTG3,103)'_Kn-ep',Kn_e_p,'_Kn-pe',Kn_p_e,'_Kn-pp',Kn_p_p
103 FORMAT (A6,ES9.2,A6,ES9.2,A6,ES9.2)
IF (ACCFLAG.EQ.0) WRITE(INTG4,104) '_ CIS'
IF (ACCFLAG.EQ.1) WRITE(INTG4,104) '_GSIS'
104 FORMAT (A5)
WRITE(INTG5,105)'THE',NPOLE
WRITE(INTG6,105)'PHI',NAZIM
105 FORMAT (A3,I3)

!Field
FILENAME = './results/2D_'//TRIM(INTG1)//'_'//TRIM(INTG2)//'_'//TRIM(INTG5)//'_'//TRIM(INTG6)//TRIM(INTG3)//TRIM(INTG4)//'.dat'
OPEN(13,FILE=FILENAME,STATUS='UNKNOWN',ACTION='WRITE')
WRITE(13,*) 'VARIABLES="x","y","eT","eE","eqx","eqy","pT","pE","pqx","pqy"'
WRITE(13,*) 'ZONE I = 109 J = 109'
DO J = 1, Ny
    DO I = 1, Nx
        WRITE(13,100) px(I),py(J),peT(I,J),peE(I,J),peqx(I,J),peqy(I,J), &
                                  ppT(I,J),ppE(I,J),ppqx(I,J),ppqy(I,J)
        100 FORMAT(1X,10ES16.6)
    END DO
END DO
CLOSE(13)

!analytical solution
! OPEN(13,FILE='Conduction_A.dat',STATUS='UNKNOWN',ACTION='WRITE')
! WRITE(13,*)'VARIABLES="x","y","T","qx","qy"'
! WRITE(13,*)'ZONE I = 109 J = 109'
! DO J = 1, Ny
!     DO I = 1, Nx
!         ss = 0.d0
!         ss_x = 0.d0
!         ss_y = 0.d0
!         DO M = 1, 200
!             ss = ss + ((-1.d0)**REAL(M+1,DBL)+1.d0)/REAL(M,DBL)*sin(REAL(M,DBL)*PI*px(I)) &
!                     *sinh(REAL(M,DBL)*PI*py(J))/sinh(REAL(M,DBL)*PI)
!             ss_x = ss_x + ((-1.d0)**REAL(M+1,DBL)+1.d0)*cos(REAL(M,DBL)*PI*px(I)) &
!                     *sinh(REAL(M,DBL)*PI*py(J))/sinh(REAL(M,DBL)*PI)
!             ss_y = ss_y + ((-1.d0)**REAL(M+1,DBL)+1.d0)*sin(REAL(M,DBL)*PI*px(I)) &
!                     *cosh(REAL(M,DBL)*PI*py(J))/sinh(REAL(M,DBL)*PI)
!         END DO
!         ss = ss*2.d0/PI
!         ss_x = -ss_x*Cv/3.d0*TAU_R
!         ss_y = -ss_y*Cv/3.d0*TAU_R

!         WRITE(13,105)px(I),py(J),ss,ss_x,ss_y
!         105 FORMAT(1X,5ES16.6)
!     END DO
! END DO
! CLOSE(13)

!FILENAME = '1D_'//TRIM(INTG1)//'_'//TRIM(INTG2)//'_'//TRIM(INTG3)//TRIM(INTG4)//'.dat'
!OPEN(13,FILE=FILENAME,STATUS='UNKNOWN',ACTION='WRITE')
!WRITE(13,*) 'VARIABLES="x","T","qx","qy"'
!WRITE(13,*) 'ZONE I = 109'
!J=Ny/2+1
!DO I = 1, Nx
!    WRITE(13,105) px(I),pT(I,J),pqx(I,J),pqy(I,J)
!        105 FORMAT(1X,4ES16.6)
!END DO
!CLOSE(13)

DEALLOCATE(px,py,pm,peT,peE,peqx,peqy,ppT,ppE,ppqx,ppqy,TRIID)
END SUBROUTINE Out_Put_Conduction

SUBROUTINE Cal_Dimensional_Property ()
    IMPLICIT NONE 
    peT = peT*T_ref + T0
    ppT = ppT*T_ref + T0
    peE = peE*Ce_ref*T_ref
    ppE = ppE*Cp_ref*T_ref
    peqx = peqx*Ce_ref*T_ref*Ve_ref
    peqy = peqy*Ce_ref*T_ref*Ve_ref
    ppqx = ppqx*Cp_ref*T_ref*Vp_ref
    ppqy = ppqy*Cp_ref*T_ref*Vp_ref
END SUBROUTINE

END MODULE OUT_PUT
