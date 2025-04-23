PROGRAM GenerateGmsh
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: xmax = 1.d0, ymax = 0.2d0
INTEGER :: Nx, Ny
INTEGER :: I,J,NE,IDX
CHARACTER (LEN=50) :: FILENAME
CHARACTER (LEN=10) :: INTG1, INTG2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x, y
DOUBLE PRECISION :: s,ss3

Nx = 21
Ny = 5

WRITE(*,*)Nx,Ny

ALLOCATE(x(Nx),y(Ny))

DO I = 1, Nx
	!x(I) = (I-1.d0)*xmax/(Nx-1.d0)
	s = (I-1.d0)/(Nx-1.d0)
	x(I) = s*s*s*(10.d0-15.d0*s+6.d0*s*s)*xmax
END DO
DO J = 1, Ny
	!y(J) = (J-1.d0)*ymax/(Ny-1.d0)
	s = (J-1.d0)/(Ny-1.d0)
	y(J) = s*s*s*(10.d0-15.d0*s+6.d0*s*s)*ymax
END DO

IF (Nx.LT.10) WRITE (INTG1, 101) 'Nx',Nx
IF (Ny.LT.10) WRITE (INTG2, 101) 'Ny',Ny
101 FORMAT (A2,I1)
IF (Nx.GE.10) WRITE (INTG1, 102) 'Nx',Nx
IF (Ny.GE.10) WRITE (INTG2, 102) 'Ny',Ny
102 FORMAT (A2, I2) 
IF (Nx.GE.100) WRITE(INTG1,180) 'Nx', Nx
IF (Ny.GE.100) WRITE(INTG2,180) 'Ny', Ny
180 FORMAT (A2,I3)


FILENAME = './A2_'//TRIM(INTG1)//'_'//TRIM(INTG2)//'.msh'
OPEN (11, FILE = FILENAME, STATUS = 'UNKNOWN', ACTION = 'WRITE')
WRITE(11,103)'$MeshFormat'
103 FORMAT (A11)
WRITE(11,104)'2.2 0 8'
104 FORMAT (A7)
WRITE(11,105)'$EndMeshFormat'
105 FORMAT (A14)
WRITE(11,105)'$PhysicalNames'
WRITE(11,106)'5'
106 FORMAT (A1)
WRITE(11,107)'1 11 "SWall (physical id 11)"'
107 FORMAT (A30)
WRITE(11,107)'1 12 "NWall (physical id 12)"'
WRITE(11,107)'1 13 "WWall (physical id 13)"'
WRITE(11,107)'1 14 "EWall (physical id 14)"'
WRITE(11,107)'2 15 "Fluid (physical id 15)"'
WRITE(11,108)'$EndPhysicalNames'
108 FORMAT (A17)

!nodes
WRITE(11,109)'$Nodes'
109 FORMAT(A6)

WRITE(11,110) Nx*Ny
110 FORMAT (I6)

DO I = 1, Nx
	DO J = 1, Ny
		WRITE(11,113) J+(I-1)*Ny, x(I), y(J), 0
		113 FORMAT (I6,2F19.16,I2)
	END DO
END DO
WRITE(11,116)'$EndNodes'
116 FORMAT(A9)

!Element
WRITE(11,116)'$Elements'
NE = 2*(Nx-1) + 2*(Ny-1) + 2*(Nx-1)*(Ny-1)
WRITE(11,117) NE
117 FORMAT (I6)

IDX = 0

!SWall
J = 1
DO I = 1, Nx-1
IDX = IDX + 1
WRITE(11,118) IDX, '1', '2', '11', '5', J+(I-1)*Ny, J+I*Ny
118 FORMAT(I6,2A2,A3,A2,2I6)
END DO

!NWall
J = Ny
DO I = 1, Nx-1
IDX = IDX + 1
WRITE(11,118) IDX, '1', '2', '12', '6', J+(I-1)*Ny, J+I*Ny
END DO

!WWall
I = 1
DO J = 1, Ny-1
IDX = IDX + 1
WRITE(11, 118) IDX, '1', '2', '13', '7', J+(I-1)*Ny, J+1+(I-1)*Ny
END DO

!EWall
I = Nx
DO J = 1, Ny-1
IDX = IDX + 1
WRITE(11, 118) IDX, '1', '2', '14', '8', J+(I-1)*Ny, J+1+(I-1)*Ny
END DO



!TRIANGLE
DO I = 1, Nx/2
	DO J = 1, Ny/2
		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+(I-1)*Ny, J+I*Ny, J+1+I*Ny
		119 FORMAT (I6,2A2,2A3,3I6)

		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+(I-1)*Ny, J+1+I*Ny, J+1+(I-1)*Ny

	END DO
	DO J = Ny/2+1, Ny-1
		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+(I-1)*Ny, J+I*Ny, J+1+(I-1)*Ny

		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+1+(I-1)*Ny, J+I*Ny, J+1+I*Ny
	END DO
END DO

DO I = Nx/2+1, Nx-1
	DO J = 1, Ny/2

		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+(I-1)*Ny, J+I*Ny, J+1+(I-1)*Ny

		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+1+(I-1)*Ny, J+I*Ny, J+1+I*Ny

	END DO
	DO J = Ny/2+1, Ny-1
		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+(I-1)*Ny, J+I*Ny, J+1+I*Ny

		IDX = IDX + 1
		WRITE(11,119) IDX, '2', '2', '15', '10', J+(I-1)*Ny, J+1+I*Ny, J+1+(I-1)*Ny
	END DO
END DO

WRITE(11,120) '$EndElements'
120 FORMAT(A12)
WRITE(11,*)

CLOSE(11)


STOP
END PROGRAM
