MODULE SOLVER
USE CONSTANT
USE GLOBAL_VARIABLE
USE MATRICES
USE VELOCITY_DISTRIBUTION_FUNCTION
USE VELOCITY_GRID
USE SPATIAL_GRID
USE BOUNDARY_CONDITION
USE INTEGRATIONS
IMPLICIT NONE

CONTAINS

SUBROUTINE DG_Solver_VDF ()
IMPLICIT NONE
EXTERNAL:: DGETRF, DGETRS
INTEGER :: J1,J2,I,M,IL,L
INTEGER :: TRID, TRIDext, FCID
INTEGER :: BCID, JS1, JS2
INTEGER :: PAIR,INFO
REAL (KIND = DBL) :: speed

!$OMP PARALLEL DO PRIVATE (J1,J2,I,TRID,TRIDext,M,L,IL,INFO,BCID,FCID,JS1,JS2, &
!$OMP PAIR,speed,IPIV,Ae_SOL,Ap_SOL,Ae_SRC,Ap_SRC,FWe,FWp)
DO J2 = 1, NAZIM
	DO J1 = 1, NPOLE
!**************************************************************
		DO I = 1, N_TRIS
			TRID = TRI_ORDER(I,J1,J2)
			!calculate Ae_SOL, Ap_SOL
			Ae_SOL = 0.d0; Ap_SOL = 0.d0
			DO L = 1, NDOF_TRI
    			DO M = 1, NDOF_TRI
    				Ae_SOL(M,L) = 1.d0/Kn_e_p*INT_NODFUNC_TRI_TRI(M,L,TRID) &
    						- Ve_s*CX(J1,J2)*INT_NODFUNC_TRI_TRI_X(TRID,M,L)-Ve_s*CY(J1,J2)*INT_NODFUNC_TRI_TRI_Y(TRID,M,L)
					Ap_SOL(M,L) = 1.d0/Kn_p_C*INT_NODFUNC_TRI_TRI(M,L,TRID) &
    						- Vp_s*CX(J1,J2)*INT_NODFUNC_TRI_TRI_X(TRID,M,L)-Vp_s*CY(J1,J2)*INT_NODFUNC_TRI_TRI_Y(TRID,M,L)
    				DO IL = 1,3
    					speed = CX(J1,J2)*TRIANGLES_INF(TRID,2*IL)+CY(J1,J2)*TRIANGLES_INF(TRID,2*IL+1)
    					Ae_SOL(M,L) = Ae_SOL(M,L) + 0.5d0*Ve_s*(speed+ABS(speed))*INT_NODFUNC_TRI_FC(TRID,IL,M,L)
						Ap_SOL(M,L) = Ap_SOL(M,L) + 0.5d0*Vp_s*(speed+ABS(speed))*INT_NODFUNC_TRI_FC(TRID,IL,M,L)
    				END DO
    			END DO
    		END DO

!*****************************************************************************************************
			!calculate Ae_SRC, Ap_SRC
			Ae_SRC = 0.d0; Ap_SRC = 0.d0
			DO M = 1, NDOF_TRI
				DO L = 1, NDOF_TRI
					Ae_SRC(M) = Ae_SRC(M) + eT_s(L,TRID)*Ce_s/4.d0/PI/Kn_e_p*INT_NODFUNC_TRI_TRI(L,M,TRID)
					Ap_SRC(M) = Ap_SRC(M) + eT_s(L,TRID)*Cp_s/4.d0/PI/Kn_p_e*INT_NODFUNC_TRI_TRI(L,M,TRID) &
					                      + pT_s(L,TRID)*Cp_s/4.d0/PI/Kn_p_p*INT_NODFUNC_TRI_TRI(L,M,TRID)
				END DO
			END DO

			DO IL = 1, 3
				speed = CX(J1,J2)*TRIANGLES_INF(TRID,2*IL)+CY(J1,J2)*TRIANGLES_INF(TRID,2*IL+1)
				FCID = TRIANGLES_TAG(TRID,3+IL)
				BCID = FACES_TAG(FCID,3)

				IF (BCID.EQ.0) THEN
					!interior
					IF (FACES_TAG(FCID,5).EQ.TRID) TRIDext = FACES_TAG(FCID,6)  !neighbor triangle
					IF (FACES_TAG(FCID,6).EQ.TRID) TRIDext = FACES_TAG(FCID,5)

					DO L = 1, NDOF_TRI
						Ae_SRC = Ae_SRC - 0.5d0*Ve_s*(speed-ABS(speed))*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,L)*eVDF(L,TRIDext,J1,J2)
						Ap_SRC = Ap_SRC - 0.5d0*Vp_s*(speed-ABS(speed))*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,L)*pVDF(L,TRIDext,J1,J2)
					END DO
				ELSE
					!electron boundary condition
					IF (BC_TYP(BCID).EQ.1) THEN
						!thermalisation
						FWe = Ce_s/4.d0/PI*BC_TEMP(BCID)*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,1)
						Ae_SRC = Ae_SRC - 0.5d0*Ve_s*(speed-ABS(speed))*FWe
					END IF
					IF (BC_TYP(BCID).EQ.2) THEN
						!adiabatic
						FWe = -eFLUX_WALL(:,FCID)
						Ae_SRC = Ae_SRC - 0.5d0*Ve_s*(speed-ABS(speed))*FWe
					END IF
					!phonon boundary condition
					IF (BC_TYP_P(BCID).EQ.1) THEN
						!thermalisation
						FWp = Cp_s/4.d0/PI*BC_TEMP_P(BCID)*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,1)
						Ap_SRC = Ap_SRC - 0.5d0*Vp_s*(speed-ABS(speed))*FWp
					END IF
					IF (BC_TYP_P(BCID).EQ.2) THEN
						!adiabatic
						FWp = -pFLUX_WALL(:,FCID)
						Ap_SRC = Ap_SRC - 0.5d0*Vp_s*(speed-ABS(speed))*FWp
					END IF
					
					IF (BC_TYP(BCID).EQ.3) THEN
						!Period
						PAIR = FACES_TAG(FCID,4)
						IF (FACES_TAG(PAIR,5).EQ.0) TRIDext = FACES_TAG(PAIR,6)
						IF (FACES_TAG(PAIR,6).EQ.0) TRIDext = FACES_TAG(PAIR,5)

						DO L = 1, NDOF_TRI
							Ae_SRC = Ae_SRC - 0.5d0*Ve_s*(speed-ABS(speed))*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,L)*eVDF(L,PAIR,J1,J2)
							Ap_SRC = Ap_SRC - 0.5d0*Vp_s*(speed-ABS(speed))*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,L)*pVDF(L,PAIR,J1,J2)
						END DO
					END IF
				END IF
			END DO
!**************************************************************************************************
			!local solver
			CALL DGETRF (NDOF_TRI,NDOF_TRI,Ae_SOL,NDOF_TRI,IPIV,INFO)
			IF (INFO.NE.0) WRITE (*,*) 'Warning LU factorization error: ',INFO
    		CALL DGETRS ('N',NDOF_TRI,1,Ae_SOL,NDOF_TRI,IPIV,Ae_SRC,NDOF_TRI,INFO)
			eVDF(:,TRID,J1,J2) = Ae_SRC

			CALL DGETRF (NDOF_TRI,NDOF_TRI,Ap_SOL,NDOF_TRI,IPIV,INFO)	
			IF (INFO.NE.0) WRITE (*,*) 'Warning LU factorization error: ',INFO
			CALL DGETRS ('N',NDOF_TRI,1,Ap_SOL,NDOF_TRI,IPIV,Ap_SRC,NDOF_TRI,INFO)
			pVDF(:,TRID,J1,J2) = Ap_SRC
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

!STOP
! CLOSE(file_unit)

END SUBROUTINE DG_Solver_VDF

END MODULE
