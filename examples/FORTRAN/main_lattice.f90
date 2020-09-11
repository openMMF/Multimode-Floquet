
PROGRAM MULTIMODEFLOQUET

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK  
  USE FLOQUETINITINTERFACE
  USE ARRAYS 


  IMPLICIT NONE
  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_BARE
  INTEGER                                          INFO,m,INDEX0,r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET,E_BARE
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,H_BARE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  DOUBLE PRECISION                              :: T1,T2

  INTEGER index





  OPEN(UNIT=3,FILE="ShakenLattice_spectrum.dat",ACTION="WRITE")
  !OPEN(UNIT=4,FILE="ShakenLattice_eigenvectors.dat",ACTION="WRITE")



  INFO = 0
  D_BARE = 16
  CALL FLOQUETINIT(ID,'lattice',0.1d1*D_BARE,INFO)

  ALLOCATE(E_BARE(D_BARE))
  ALLOCATE(H_BARE(D_BARE,D_BARE))
  ALLOCATE(U_AUX(D_BARE,D_BARE))
  ALLOCATE(MODES_NUM(2))
  
  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 1 !(DRIVING BY ONE HARMONIC)
  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(D_BARE,D_BARE))
     FIELDS(m)%V = 0.0
  END DO
  

  ! STATIC COMPONENT OF THE HAMILTONIAN
  ! E.G. A 1D LATTICE WITH TWO SITES IN A UNIT CELL
  DO m=2,D_BARE,2 ! EVEN SITES
     FIELDS(1)%V(m,m)  = 0.25
  END DO
  DO m=1,D_BARE,2 ! ODD SITES
     FIELDS(1)%V(m,m)  = -0.25
  END DO
  ! AND NEAR NEIGHBOURGH COUPLING
  DO m=1,D_BARE-1
     FIELDS(1)%V(m,m+1)  = 0.125*EXP(DCMPLX(0.0, 1.0)*2*PI*m/D_BARE)!EXP(DCMPLX(0.0, 1.0)*COS(2*PI*m/D_BARE))
     FIELDS(1)%V(m+1,m)  = 0.125*EXP(DCMPLX(0.0,-1.0)*2*PI*m/D_BARE)!EXP(DCMPLX(0.0,-1.0)*COS(2*PI*m/D_BARE))
  END DO
  FIELDS(1)%V(1,D_BARE) = 0.125*EXP(DCMPLX(0.0, 1.0)*2*PI/D_BARE) ! WITH PERIODIC BOUNDARY CONDITIONS
  FIELDS(1)%V(D_BARE,1) = 0.125*EXP(DCMPLX(0.0,-1.0)*2*PI/D_BARE) ! WITH PERIODIC BOUNDARY CONDITIONS
  FIELDS(1)%omega     = 0.0
  FIELDS(1)%N_Floquet = 0
  H_BARE = FIELDS(1)%V

  CALL LAPACK_FULLEIGENVALUES(FIELDS(1)%V,D_BARE,E_BARE,INFO)
  

  U_aux = 0.0
  index = 1
  DO m=1,D_BARE/2,2
!     WRITE(*,*) 2.0*pi*m/D_bare,E_bare(m)!m*4*pi/D_bare - pi,E_BARE(m)
     U_AUX(:,index)  = FIELDS(1)%V(:,m)
     index = index + 1
  END DO
!  WRITE(*,*)
  DO m=D_BARE/2,2,-2
!     WRITE(*,*) 2.0*(pi-pi*(m-1)/D_bare),E_bare(m)!m*4*pi/D_bare - pi,E_BARE(m)
     U_AUX(:,index)  = FIELDS(1)%V(:,m)
     index = index + 1
  END DO
!  write(*,*)
  DO m=D_BARE/2+1,D_BARE,2
!     WRITE(*,*) 2.0*(pi*m/D_bare-pi/2),E_bare(m)!m*4*pi/D_bare - pi,E_BARE(m)
     U_AUX(:,index)  = FIELDS(1)%V(:,m)
     index = index + 1
  END DO
!  WRITE(*,*)
  DO m=D_BARE,D_BARE/2+2,-2!D_BARE/2+2,D_BARE,2
!     WRITE(*,*) 2.0*(2*pi-pi*(m-1)/D_bare - pi/2.0),E_bare(m)!m*4*pi/D_bare - pi,E_BARE(m)
     U_AUX(:,index)  = FIELDS(1)%V(:,m)
     index = index + 1
  END DO
!  write(*,*)
!  write(*,*)
!  CALL WRITE_MATRIX(abs(FIELDS(1)%V))
!  CALL WRITE_MATRIX(abs(u_aux))
!  U_aux = fields(1)%V
  ! BUILD THE STATIC HAMILTONIANN IN THE DIAGONAL BASIS
  FIELDS(1)%V = MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(H_BARE,U_AUX))

!  CALL WRITE_MATRIX(real(FIELDS(1)%V))
  

!  DO m=1,D_BARE
!     write(*,*) real(FIELDS(1)%V(m,m))
!  END DO

!  CALL LAPACK_FULLEIGENVALUES(FIELDS(1)%V,D_BARE,E_BARE,INFO)
!  write(*,*)
!  write(*,*)
!  DO m=1,D_BARE
!     write(*,*) E_BARE(m)
!  END DO

  ! DEFINE HARMONICS COUPLINGS IN THE ORIGINAL BASIS: FIELD ONE
  DO m=1,D_BARE-1
     FIELDS(2)%V(m,m+1)  = 0.125/2
     FIELDS(2)%V(m+1,m)  = 0.125/2
  END DO
!  FIELDS(2)%V(1,D_BARE) = 0.125/2 ! WITH PERIODIC BOUNDARY CONDITIONS
!  FIELDS(2)%V(D_BARE,1) = 0.125/2 ! WITH PERIODIC BOUNDARY CONDITIONS
  FIELDS(2)%omega     = 0.5
  FIELDS(2)%N_Floquet = 6

  ! TRANSFORM TO THE BASIS OF EIGNSTATES OF THE STATIC SYSTEM
  FIELDS(2)%V = MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(FIELDS(2)%V,U_AUX))

!  write(*,*)
!  write(*,*)
  
!  CALL WRITE_MATRIX(REAL(FIELDS(2)%V))

  DO m=1,512,16
     ! --- SET DRIVING PARAMETERS 
     FIELDS(2)%omega = 0.4 + (m-1.0)*0.4/512.0
     !write(*,*) total_frequencies,ID%ID_SYSTEM
     !--- FIND THE MULTIMODE FLOQUET SPECTRUM 
     CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
     !write(3,*)
     !call write_matrix(real(h_floquet))
     ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
     ALLOCATE(U_F(SIZE(H_FLOQUET,1),SIZE(H_FLOQUET,1)))
     E_FLOQUET = 0.0   
     !     DO r=1,SIZE(H_FLOQUET,1)
     !        write(*,*) real(H_FLOQUET(r,r))
     !     END DO
     !     write(*,*)
     !     write(*,*)
     
     CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_FLOQUET,INFO)
     U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
     DEALLOCATE(H_FLOQUET)
     DO r=1,SIZE(U_F,1)
        WRITE(3,*) FIELDS(2)%omega,r,E_FLOQUET(r)
     END DO
!     CALL WRITE_MATRIX(ABS(U_F))
!     write(*,*)
!     write(*,*)
     
!     CALL WRITE_MATRIX(ABS(U_F))
     DEALLOCATE(E_FLOQUET)
     DEALLOCATE(U_F)
  END DO
  
  
END PROGRAM MULTIMODEFLOQUET

