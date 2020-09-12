
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
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,dh2_dt,dh3_dt,U_dt
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG,W2,W3
  DOUBLE PRECISION                              :: T1,T2

  DOUBLE PRECISION m_, eta,  phi1,  phi2 ,  gamma, dt
  INTEGER N_






  OPEN(UNIT=3,FILE="qubit_lattice_DRIVER_.dat",ACTION="WRITE")


  INFO = 0
  CALL FLOQUETINIT(ID,'qubit',INFO)  
  
  D_BARE = ID%D_BARE
  ALLOCATE(P_AVG(D_BARE,D_BARE))
  ALLOCATE(U_AUX(D_BARE,D_BARE))
  ALLOCATE(U_dt(D_BARE,D_BARE))
  ALLOCATE(dh2_dt(D_BARE,D_BARE))
  ALLOCATE(dh3_dt(D_BARE,D_BARE))
  ALLOCATE(W2(D_BARE,D_BARE))
  ALLOCATE(W3(D_BARE,D_BARE))

  ALLOCATE(MODES_NUM(3))
  
  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 1 !(DRIVING BY ONE HARMONIC)
  MODES_NUM(3) = 1 !(DRIVING BY ONE HARMONIC)
  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(ID%D_BARE,ID%D_BARE))
  END DO
  
  !eq. (20), arxiv 1612.02143
  m_    = 2.2
  eta   = 2.0
  phi1  = pi/10.0
  phi2  = 0.0
  gamma = 0.5*(1.0 + sqrt(5.0))

  FIELDS(1)%X     = 0.0
  FIELDS(1)%Y     = 0.0
  FIELDS(1)%Z     = 2.0*m_*eta
  FIELDS(1)%phi_x = 0.0
  FIELDS(1)%phi_y = 0.0
  FIELDS(1)%phi_z = 0.0
  FIELDS(1)%omega = 0.0
  FIELDS(1)%N_Floquet = 0

  FIELDS(2)%X     =  2.0*eta
  FIELDS(2)%Y     =  0.0
  FIELDS(2)%Z     = -2.0*eta
  FIELDS(2)%phi_x = phi1 - pi/2.0
  FIELDS(2)%phi_y = 0.0
  FIELDS(2)%phi_z = phi1 
  FIELDS(2)%omega = 0.1
  FIELDS(2)%N_Floquet = 2

  FIELDS(3)%X     =  0.0
  FIELDS(3)%Y     =  2.0*eta
  FIELDS(3)%Z     = -2.0*eta
  FIELDS(3)%phi_x = 0.0
  FIELDS(3)%phi_y = phi2 - pi/2
  FIELDS(3)%phi_z = phi2
  FIELDS(3)%omega = gamma*FIELDS(2)%OMEGA
  FIELDS(3)%N_Floquet = 2
 
  DO m=1,TOTAL_FREQUENCIES    
     FIELDS(m)%X = FIELDS(m)%X*exp(DCMPLX(0.0,1.0)*FIELDS(m)%phi_x)
     FIELDS(m)%Y = FIELDS(m)%Y*exp(DCMPLX(0.0,1.0)*FIELDS(m)%phi_y)
     FIELDS(m)%Z = FIELDS(m)%Z*exp(DCMPLX(0.0,1.0)*FIELDS(m)%phi_z)
  END DO

  !--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
  U_aux = 0.0
  U_AUX(1,1) = 1.0
  U_AUX(2,2) = 1.0
  
  T1 = 0.0
  T2 = 1.0E-4
  dt = T2
  CALL TIMEEVOLUTIONOPERATOR(ID,D_BARE,SIZE(MODES_NUM,1),SIZE(MODES_NUM,1),MODES_NUM,FIELDS,T1,T2,U_dt,INFO) 

  dt = T2
  N_ = 1024
  DO r=1,N_
     T2 = (r-1.0)*10.0E4/N_
     CALL TIMEEVOLUTIONOPERATOR(ID,D_BARE,SIZE(MODES_NUM,1),TOTAL_FREQUENCIES,MODES_NUM,FIELDS,T1,T2,U_AUX,INFO) 
     dh2_dt = 2.0*eta*FIELDS(2)%OMEGA*(                 J_x*cos(fields(2)%omega*T2+phi1) + J_z*sin(fields(2)%omega*T2+phi1))
     dh3_dt = 2.0*eta*FIELDS(3)%OMEGA*(DCMPLX(0.0,-1.0)*J_y*cos(fields(3)%omega*T2+phi2) + J_z*sin(fields(3)%omega*T2+phi2))
     W2 = W2 + dt*MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(dh2_dt,U_AUX))
     W3 = W3 + dt*MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(dh3_dt,U_AUX))     
     WRITE(3,*) t2,FIELDS(2)%OMEGA,ABS(U_AUX)**2,W2,W3
  END DO
  WRITE(3,*)
  
  
END PROGRAM MULTIMODEFLOQUET


!!$SUBROUTINE MATRIXEXPONENTIATION(D,H,dt,U,info)
!!$
!!$  IMPLICIT NONE
!!$  INTEGER, INTENT(IN) :: D
!!$  INTEGER, INTENT(INOUT) :: INFO
!!$  COMPLEX*16,DIMENSION(D,D), INTENT(IN) :: H
!!$  DOUBLE PRECISION, INTENT(IN):: dt
!!$  COMPLEX*16,DIMENSION(D,D), INTENT(OUT) :: U
!!$
!!$  COMPLEX*16,DIMENSION(D,D) :: U_AUX  
!!$  INTEGER N,i,j
!!$  DOUBLE PRECISION :: i_factorial
!!$  
!!$  N = 6
!!$  
!!$  U = 0
!!$  U_AUX = 0
!!$  DO i=1,D
!!$     U(i,i) = 1.0
!!$     U_AUX(i,i) = 1.0
!!$  END DO
!!$
!!$  DO i=1,N
!!$     U_AUX = MATMUL(U_AUX,DCMPLX(0.0,-1.0)*H*dt)
!!$     i_factorial = 1.0
!!$     DO j=1,i
!!$        i_factorial = i_factorial*j
!!$     END DO
!!$     U = U + U_AUX/i_factorial
!!$  END DO
!!$
!!$END SUBROUTINE MATRIXEXPONENTIATION
