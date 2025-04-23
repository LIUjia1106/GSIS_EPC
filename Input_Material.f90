! Define the equilibrium distribution and derivatives of temperature for electron
REAL(KIND=DBL) FUNCTION FERMI_DIRAC_DISTRIBUTION(En,Mu,T)
USE CONSTANT
USE GLOBAL_VARIABLE
IMPLICIT NONE
REAL (KIND=DBL), INTENT(IN) :: En, Mu, T
REAL(KIND=DBL) :: ANS
ANS = 1.d0/(EXP((En-Mu)/(BOLTZ*T))+1.d0)
FERMI_DIRAC_DISTRIBUTION = ANS
RETURN
END FUNCTION FERMI_DIRAC_DISTRIBUTION

REAL(KIND=DBL) FUNCTION PARTIAL_FERMI_DIRAC_T(En,Mu,T)
USE CONSTANT
USE GLOBAL_VARIABLE
IMPLICIT NONE
REAL (KIND=DBL), INTENT(IN) :: En, Mu, T
REAL(KIND=DBL) :: ANS
ANS = EXP((En-Mu)/(BOLTZ*T))/(EXP((En-Mu)/(BOLTZ*T))+1.d0)**2.d0*(En-Mu)/(BOLTZ*T**2.d0)
PARTIAL_FERMI_DIRAC_T = ANS
RETURN
END FUNCTION PARTIAL_FERMI_DIRAC_T

! Define the equilibrium distribution and derivatives of temperature for phonon
REAL(KIND=DBL) FUNCTION BOSE_EINSTEIN_DISTRIBUTION(Freq,T)
USE CONSTANT
IMPLICIT NONE
REAL (KIND=DBL), INTENT(IN) :: Freq, T
REAL(KIND=DBL) :: ANS
ANS = 1.d0/(EXP((DIRAC*Freq)/(BOLTZ*T))-1.d0)
BOSE_EINSTEIN_DISTRIBUTION = ANS
RETURN
END FUNCTION BOSE_EINSTEIN_DISTRIBUTION

REAL(KIND=DBL) FUNCTION PARTIAL_BOSE_EINSTEIN_T(Freq,T)
USE CONSTANT
IMPLICIT NONE
REAL (KIND=DBL), INTENT(IN) :: Freq, T
REAL(KIND=DBL) :: ANS
ANS = EXP((DIRAC*Freq)/(BOLTZ*T))/(EXP((DIRAC*Freq)/(BOLTZ*T))-1.d0)**2.d0*(DIRAC*Freq)/(BOLTZ*T**2.d0)
PARTIAL_BOSE_EINSTEIN_T = ANS
RETURN
END FUNCTION PARTIAL_BOSE_EINSTEIN_T

! Define the function of heat capacity of electron and phonon
REAL(KIND=DBL) FUNCTION ELECTRON_HEAT_CAPACITY(En,Mu,De,T)
USE CONSTANT
USE GLOBAL_VARIABLE
IMPLICIT NONE
REAL (KIND=DBL), INTENT(IN) :: En, Mu, De, T
REAL(KIND=DBL) :: PARTIAL_FERMI_DIRAC_T
REAL(KIND=DBL) :: ANS
ANS = ABS(En-Mu)*De*PARTIAL_FERMI_DIRAC_T(En,Mu,T)
ELECTRON_HEAT_CAPACITY = ANS
RETURN
END FUNCTION ELECTRON_HEAT_CAPACITY

REAL(KIND=DBL) FUNCTION PHONON_HEAT_CAPACITY(Freq,Dph,T)
USE CONSTANT
USE GLOBAL_VARIABLE
IMPLICIT NONE
REAL (KIND=DBL), INTENT(IN) :: Freq, Dph, T
REAL(KIND=DBL) :: PARTIAL_BOSE_EINSTEIN_T
REAL(KIND=DBL) :: ANS
ANS = DIRAC*Freq*Dph*PARTIAL_BOSE_EINSTEIN_T(Freq,T)
PHONON_HEAT_CAPACITY = ANS
RETURN
END FUNCTION PHONON_HEAT_CAPACITY

MODULE MATERIAL_PROPERTY
USE CONSTANT
USE GLOBAL_VARIABLE
IMPLICIT NONE

INTEGER :: nband, pbranch               !# of electron band and phonon branch
INTEGER :: nkpoint, nqpoint             !# of electron k point and phonon q point
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: electron_energy, dEn                                            !electron energy band structure
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: electron_groupvel_x, electron_groupvel_y, electron_groupvel_z   !electron group velocity
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: electron_groupspeed, electron_DOS, fermi_energy                                !electron density of states
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: phonon_frequency, dFreq                                         !phonon dispersion relation
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: phonon_groupvel_x, phonon_groupvel_y, phonon_groupvel_z         !phonon group velocity
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: phonon_groupspeed, phonon_DOS, gruneisen_param                  !phonon density of states
REAL (KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:) :: scattering_matrix                                             !$C_{e-p}(\epsilon,\epsilon',\omega)$
REAL (KIND=DBL) :: mass_ave, temp_Debye    
REAL (KIND=DBL) :: Te = 0.d0, Tph = 0.d0                                   

CONTAINS
    SUBROUTINE Init_Material_Property()
        IMPLICIT NONE
        INTEGER :: ios, i, j, bandID, kpointID, qpoint_x, qpoint_y, qpoint_z, branchID
        CHARACTER (LEN=100) line     

        !=========================================================================================
        ! Get electron band structure and group velocity
        !=========================================================================================
        WRITE (*,*) '**** Initial Electron Property ****'
        ! open the input file
        OPEN (UNIT=12, FILE=TRIM(FNAME_ELECTRON_BAND),STATUS='OLD',ACTION='READ',IOSTAT=ios)
        ! read header lines to extract nband and nkpoint
        READ (12,*) !skip 1st line
        READ (12, '(A)') line
        READ (line(INDEX(line, ':')+1:), *) nband
        READ (12, '(A)') line
        READ (line(INDEX(line, ':')+1:), *) nkpoint
        ! WRITE(*,*) nband, nkpoint

        ! read electron energy and group velocity components into arrays
        ALLOCATE (electron_energy(nkpoint,nband),dEn(nkpoint,nband))
        ALLOCATE (electron_groupvel_x(nkpoint,nband),electron_groupvel_y(nkpoint,nband),electron_groupvel_z(nkpoint,nband))
        ALLOCATE (electron_groupspeed(nkpoint,nband),electron_DOS(nkpoint,nband),fermi_energy(nkpoint,nband))
        DO i = 1,nkpoint
            DO j = 1,nband
                READ(12,*,IOSTAT=ios) kpointID, bandID, electron_energy(i,j), electron_groupvel_x(i,j), & 
                                      electron_groupvel_y(i,j), electron_groupvel_z(i,j), electron_DOS(i,j), fermi_energy(i,j)
                IF (ios.NE.0) WRITE(*,*) 'Error to read electron energy band structure: ', ios, j, i
                electron_groupspeed(i,j) = sqrt(electron_groupvel_x(i,j)**2.d0+electron_groupvel_y(i,j)**2.d0+electron_groupvel_z(i,j)**2.d0)
            END DO
        END DO 
        CLOSE(12)

        DO j = 1,nband
            dEn(1,j)=(electron_energy(2,j)-electron_energy(1,j))/2.d0
            DO i=2,nkpoint-1
                dEn(i,j)=(electron_energy(i,j)-electron_energy(i-1,j))/2.d0 & 
                       + (electron_energy(i+1,j)-electron_energy(i,j))/2.d0
            END DO
            dEn(nkpoint,j)=(electron_energy(nkpoint,j)-electron_energy(nkpoint-1,j))/2.d0
        END DO
        
        CALL Cal_Electron_GroupVelocity()
        CALL Cal_Electron_HeatCapacity()

        !=========================================================================================
        ! Get phonon dispersion relation and group velocity
        !=========================================================================================
        WRITE (*,*) '**** Initial Phonon Property ****'
        ! open the input file
        OPEN (UNIT=13, FILE=TRIM(FNAME_PHONON_DISPERSION),STATUS='OLD',ACTION='READ',IOSTAT=ios)
        ! read header lines to extract pbranch and nqpoint
        READ (13,*) !skip 1st line
        READ (13, '(A)') line
        READ (line(INDEX(line, ':')+1:), *) pbranch
        READ (13, '(A)') line
        READ (line(INDEX(line, ':')+1:), *) nqpoint
        ! WRITE(*,*) pbranch,nqpoint

        ! read phonon frequency and group velocity components into arrays
        ALLOCATE (phonon_frequency(nqpoint,pbranch),dFreq(nqpoint,pbranch))
        ALLOCATE (phonon_groupvel_x(nqpoint,pbranch),phonon_groupvel_y(nqpoint,pbranch),phonon_groupvel_z(nqpoint,pbranch))
        ALLOCATE (phonon_groupspeed(nqpoint,pbranch),phonon_DOS(nqpoint,pbranch),gruneisen_param(nqpoint,pbranch))
        DO i = 1,nqpoint
            DO j = 1,pbranch
                READ(13,*,IOSTAT=ios) bandID, qpoint_x, qpoint_y, qpoint_z, phonon_frequency(i,j), phonon_groupvel_x(i,j), & 
                                      phonon_groupvel_y(i,j), phonon_groupvel_z(i,j), phonon_DOS(i,j), gruneisen_param(i,j)
                IF (ios.NE.0) WRITE(*,*) 'Error to read phonon dispersion relation: ', ios, j, i
            END DO
        END DO
        CLOSE(13)

        DO j = 1,pbranch
            dFreq(1,j)=(phonon_frequency(2,j)-phonon_frequency(1,j))/2.d0
            DO i=2,nqpoint-1
                dFreq(i,j)=(phonon_frequency(i,j)-phonon_frequency(i-1,j))/2.d0 & 
                         + (phonon_frequency(i+1,j)-phonon_frequency(i,j))/2.d0
            END DO
            dFreq(nqpoint,j)=(phonon_frequency(nqpoint,j)-phonon_frequency(nqpoint-1,j))/2.d0
        END DO

        CALL Cal_Phonon_GroupVelocity()
        CALL Cal_Phonon_HeatCapacity()

        !=========================================================================================
        ! Get electron-phonon scattering matrix
        !=========================================================================================
        WRITE (*,*) '**** Initial Electron-Phonon Scattering Matrix ****'
        ALLOCATE (scattering_matrix(nkpoint*nband,2,nqpoint*pbranch))
        OPEN (UNIT=14, FILE=TRIM(FNAME_SCATTERING_MATRIX),STATUS='OLD',ACTION='READ',IOSTAT=ios)
        DO i = 1,nkpoint
            DO bandID = 1,nband
                DO j = 1,nqpoint
                    DO branchID = 1,pbranch
                        scattering_matrix((i-1)*nband+bandID,1,(j-1)*pbranch+branchID) = 0.d0
                        scattering_matrix((i-1)*nband+bandID,1,(j-1)*pbranch+branchID) = 0.d0
                    END DO
                END DO
            END DO
        END DO

        CALL Cal_TAU_e_p()
        CALL Cal_TAU_p_e()
        CALL Cal_TAU_p_p()

        DEALLOCATE(electron_energy,dEn,STAT=ios)
        DEALLOCATE(electron_groupvel_x,electron_groupvel_y,electron_groupvel_z,STAT=ios)
        DEALLOCATE(electron_groupspeed,electron_DOS,fermi_energy,STAT=ios)
        DEALLOCATE(phonon_frequency,dFreq,STAT=ios)
        DEALLOCATE(phonon_groupvel_x,phonon_groupvel_y,phonon_groupvel_z,STAT=ios)
        DEALLOCATE(phonon_groupspeed,phonon_DOS,gruneisen_param,STAT=ios)
        DEALLOCATE(scattering_matrix,STAT=ios)

        CALL Cal_Nondimensional_Parameters()

        WRITE (*,104) 'Magnitude of electron group velocity: ', Ve
        104 FORMAT (1X, A60, ES14.2)
        WRITE (*,104) 'Volumetric specific heat of electron: ', Ce        
        WRITE (*,*)

        WRITE (*,104) 'Magnitude of phonon group velocity: ', Vp
        WRITE (*,104) 'Volumetric specific heat of phonon: ', Cp
        WRITE (*,*)
      
        WRITE (*,104) 'tau of e-p scattering at T0: ', TAU_e_p
        WRITE (*,104) 'tau of p-e scattering at T0: ', TAU_p_e
        WRITE (*,104) 'tau of p-p scattering at T0: ', TAU_p_p
        WRITE (*,*)

        WRITE (*,104) 'Kn of e-p scattering:', Kn_e_p
        WRITE (*,104) 'Kn of p-e scattering:', Kn_p_e
        WRITE (*,104) 'Kn of p-p scattering:', Kn_p_p

    END SUBROUTINE

    SUBROUTINE Cal_Nondimensional_Parameters()
        IMPLICIT NONE
        ! first determine reference values
        Ce_ref = Ce
        Cp_ref = Ce
        Ve_ref = Ve
        Vp_ref = Ve ! here we adopt electron property as reference
        ! then calculate dimensionless parameter
        Ve_s = Ve / Ve_ref
        Vp_s = Vp / Vp_ref
        Ce_s = Ce / Ce_ref
        Cp_s = Cp / Cp_ref
        Kn_e_p = TAU_e_p*Ve_ref/L_ref
        Kn_p_e = TAU_p_e*Vp_ref/L_ref
        Kn_p_p = TAU_p_p*Vp_ref/L_ref
        Kn_p_C = 1.d0/(1.d0/Kn_p_e+1.d0/Kn_p_p)
    END SUBROUTINE

    SUBROUTINE Cal_Electron_GroupVelocity()
        IMPLICIT NONE    
        INTEGER :: i, j    
        REAL (KIND=DBL) :: s, ss
        s=0.d0; ss=0.d0
        ! calculate the magnitude of electron group velocity
        DO j = 1,nband
            DO i = 1,nkpoint
                s = s + dEn(i,j)
                ss = ss + electron_groupspeed(i,j)*dEn(i,j)
            END DO
        END DO
        Ve = ss / s
    END SUBROUTINE

    SUBROUTINE Cal_Phonon_GroupVelocity()
        IMPLICIT NONE
        INTEGER :: i,j
        REAL (KIND=DBL) :: s, ss
        s=0.d0; ss=0.d0
        ! calculate the magnitude of phonon group velocity
        DO j = 1,pbranch
            DO i = 1,nqpoint
                s = s + dFreq(i,j)
                ss = ss + phonon_groupspeed(i,j)*dFreq(i,j)
            END DO
        END DO
        Vp = ss / s
    END SUBROUTINE

    SUBROUTINE Cal_Electron_HeatCapacity()
        IMPLICIT NONE
        INTEGER :: i, j
        REAL(KIND=DBL) :: ELECTRON_HEAT_CAPACITY
        Ce = 0.d0
        DO j = 1,nband
            DO i = 1,nkpoint
                Ce = Ce + ELECTRON_HEAT_CAPACITY(electron_energy(i,j),fermi_energy(i,j),electron_DOS(i,j),T0)*dEn(i,j)
            END DO           
        END DO
    END SUBROUTINE

    SUBROUTINE Cal_Phonon_HeatCapacity()
        IMPLICIT NONE
        INTEGER :: i, j
        REAL(KIND=DBL) :: PHONON_HEAT_CAPACITY , BOSE_EINSTEIN_DISTRIBUTION  
        Cp = 0.d0; 
        DO j = 1,pbranch
            DO i = 1,nqpoint
                Cp = Cp + PHONON_HEAT_CAPACITY(phonon_frequency(i,j),phonon_DOS(i,j),T0)*dFreq(i,j)    
            END DO
        END DO
    END SUBROUTINE

    SUBROUTINE Cal_TAU_e_p()
        IMPLICIT NONE
        INTEGER :: ie, je, iph, jph
        REAL (KIND=DBL) :: tau, ss
        REAL(KIND=DBL) :: FERMI_DIRAC_DISTRIBUTION,BOSE_EINSTEIN_DISTRIBUTION,ELECTRON_HEAT_CAPACITY
        DO je = 1,nband
            DO ie=1,nkpoint
                tau = 0.d0
                DO jph = 1,pbranch
                    DO iph = 1,nqpoint
                        tau = tau + 2.d0*PI*((BOSE_EINSTEIN_DISTRIBUTION(phonon_frequency(iph,jph),T0)+1-FERMI_DIRAC_DISTRIBUTION(electron_energy(ie,je)-DIRAC*phonon_frequency(iph,jph),fermi_energy(ie,je),T0))*scattering_matrix(ie,1,iph) &
                                   + (BOSE_EINSTEIN_DISTRIBUTION(phonon_frequency(iph,jph),T0)+FERMI_DIRAC_DISTRIBUTION(electron_energy(ie,je)+DIRAC*phonon_frequency(iph,jph),fermi_energy(ie,je),T0))*scattering_matrix(ie,2,iph))*dFreq(iph,jph)
                    END DO
                END DO
                tau = 1.d0/tau
                ss = ss + tau*ELECTRON_HEAT_CAPACITY(electron_energy(ie,je),fermi_energy(ie,je),electron_DOS(ie,je),T0)*dEn(ie,je)
            END DO
        END DO
        TAU_e_p = ss / Ce
    END SUBROUTINE

    SUBROUTINE Cal_TAU_p_e()
        IMPLICIT NONE
        INTEGER :: ie, je, iph, jph
        REAL (KIND=DBL) :: tau, s, ss
        REAL(KIND=DBL) :: FERMI_DIRAC_DISTRIBUTION,PHONON_HEAT_CAPACITY
        DO jph = 1,pbranch
            DO iph=1,nqpoint
                tau = 0.d0
                DO je = 1,nband
                    DO ie = 1,nkpoint
                        tau = tau + 2.d0*PI*electron_DOS(ie,je)/phonon_DOS(iph,jph) &
                            * (FERMI_DIRAC_DISTRIBUTION(electron_energy(ie,je),fermi_energy(ie,je),T0)-FERMI_DIRAC_DISTRIBUTION(electron_energy(ie,je)+DIRAC*phonon_frequency(iph,jph),fermi_energy(ie,je),T0)) &
                            * (scattering_matrix(ie,2,iph))*dEn(ie,je)
                    END DO
                END DO
                tau = 1.d0/tau
                ss = ss + tau*PHONON_HEAT_CAPACITY(phonon_frequency(iph,jph),phonon_DOS(iph,jph),T0)*dFreq(iph,jph)
            END DO
        END DO
        TAU_p_e = ss / Cp
    END SUBROUTINE

    SUBROUTINE Cal_TAU_p_p()
        IMPLICIT NONE
        INTEGER :: i, j
        REAL (KIND=DBL) :: tau, s, ss, vg2
        REAL(KIND=DBL) :: PHONON_HEAT_CAPACITY
        s = 0.d0; ss = 0.d0
        DO j = 1,pbranch
            DO i = 1,nqpoint
                vg2 = phonon_groupspeed(i,j)**2.d0
                tau = DIRAC*gruneisen_param(i,j)**2.d0/(mass_ave*temp_Debye*vg2) &
                    * phonon_frequency(i,j)**2.d0*T0*EXP(-temp_Debye/(3.d0*T0))
                tau = 1.d0/tau
                ss = ss + tau*PHONON_HEAT_CAPACITY(phonon_frequency(i,j),phonon_DOS(i,j),T0)*dFreq(i,j)
            END DO
        END DO
        TAU_p_p = ss / Cp        
    END SUBROUTINE

END MODULE MATERIAL_PROPERTY