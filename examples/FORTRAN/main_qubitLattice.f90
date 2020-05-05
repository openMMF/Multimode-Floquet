
PROGRAM MULTIMODEFLOQUET

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE ARRAYS 
  USE FLOQUETINITINTERFACE


  IMPLICIT NONE
  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_BARE
  INTEGER                                          INFO,m,INDEX0,r,i_
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET,E_0
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,dh2,dh3,U_dt,Qubit_IDENTITY,H_0
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG,W2,W3
  DOUBLE PRECISION                              :: T1,T2

  DOUBLE PRECISION m_, eta,  phi1,  phi2 ,  gamma, dt,TotalTime,W2_0,W3_0






  OPEN(UNIT=3,FILE="qubit_lattice.dat",ACTION="WRITE")


  INFO = 0
  CALL FLOQUETINIT(ID,'qubit',INFO)  
  
  D_BARE = ID%D_BARE

  ALLOCATE(H_0(D_BARE,D_BARE))
  ALLOCATE(E_0(D_BARE))
  ALLOCATE(P_AVG(D_BARE,D_BARE))
  ALLOCATE(U_AUX(D_BARE,D_BARE))
  ALLOCATE(U_dt(D_BARE,D_BARE))
  ALLOCATE(dh2(D_BARE,D_BARE))
  ALLOCATE(dh3(D_BARE,D_BARE))
  ALLOCATE(W2(D_BARE,D_BARE))
  ALLOCATE(W3(D_BARE,D_BARE))
  ALLOCATE(Qubit_IDENTITY(D_BARE,D_BARE))

  ALLOCATE(MODES_NUM(3))
  
  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 1 !(DRIVING BY ONE HARMONIC)
  MODES_NUM(3) = 1 !(DRIVING BY ONE HARMONIC)
  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(ID%D_BARE,ID%D_BARE))
  END DO

  Qubit_IDENTITY = 0.0
  Qubit_IDENTITY(1,1) = 1.0
  Qubit_IDENTITY(2,2) = 1.0
  
  !eq. (20), arxiv 1612.02143
  DO i_ = 1,1
     m_    = 1.2 + (i_-1)*1.0/5.0
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
     FIELDS(2)%N_Floquet = 13
     
     FIELDS(3)%X     =  0.0
     FIELDS(3)%Y     =  2.0*eta
     FIELDS(3)%Z     = -2.0*eta
     FIELDS(3)%phi_x = 0.0
     FIELDS(3)%phi_y = phi2 - pi/2
     FIELDS(3)%phi_z = phi2
     FIELDS(3)%omega = gamma*FIELDS(2)%OMEGA
     FIELDS(3)%N_Floquet = 13
     
     CALL SETHAMILTONIANCOMPONENTS(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)

     ! ----- INITIAL STATE
     H_0 = (2.0*eta*m_ - 2.0*COS(phi1) - 2.0*COS(phi2))*J_z + eta*2.0*SIN(phi1)*J_x + eta*2.0*SIN(phi2)*DCMPLX(0.0,1.0)*J_y
     CALL LAPACK_FULLEIGENVALUES(H_0,SIZE(H_0,1),E_0,INFO)
     WRITE(*,*) H_0(:,1)
     WRITE(*,*) E_0
    
     !--- FIND THE MULTIMODE FLOQUET SPECTRUM 
     CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
     ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
     ALLOCATE(U_F(SIZE(H_FLOQUET,1),SIZE(H_FLOQUET,1)))
     E_FLOQUET = 0.0  
     CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_FLOQUET,INFO)
     U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
     DEALLOCATE(H_FLOQUET)
     
     !--- EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
     !P_AVG = 0.0
     !CALL MULTIMODETRANSITIONAVG(SIZE(U_F,1),size(MODES_NUM,1),FIELDS,MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,P_AVG,INFO)   
     !WRITE(*,*) FIELDS(2)%omega,P_AVG
     
     !--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
     T1 = 0.0
     dt = 1.0E-3
     T2 = T1 + dt
     CALL MULTIMODETIMEEVOLUTINOPERATOR(SIZE(U_F,1),SIZE(MODES_NUM,1),MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,FIELDS,T1,T2,U_dt,INFO) 
     U_AUX = Qubit_IDENTITY
     TotalTime = 10000.0
     DO r=1,256*4096*8!int(TotalTime/dt)
        dh2 = 2.0*eta*FIELDS(2)%OMEGA*(                 J_x*cos(fields(2)%omega*T2+phi1) + J_z*sin(fields(2)%omega*T2+phi1))!*dt
        dh3 = 2.0*eta*FIELDS(3)%OMEGA*(DCMPLX(0.0,-1.0)*J_y*cos(fields(3)%omega*T2+phi2) + J_z*sin(fields(3)%omega*T2+phi2))!*dt
        U_AUX = MATMUL(U_dt,U_AUX)
        W2 = W2 + MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(dh2,U_AUX))
        W3 = W3 + MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(dh3,U_AUX))     
        T2 = T2 + dt
        !IF(MOD(r,2000).EQ.0) WRITE(3,*) i_,t2,ABS(U_AUX)**2,W2,W3
        W2_0 = DOT_PRODUCT(H_0(:,1),MATMUL(W2,H_0(:,1)))
        W3_0 = DOT_PRODUCT(H_0(:,1),MATMUL(W3,H_0(:,1)))

        IF(MOD(r,256).EQ.0) WRITE(3,*) i_,t2,ABS(U_AUX)**2,W2,W3,W2_0,W3_0
     END DO
     WRITE(3,*)
     WRITE(3,*)
     
!!$  T1 = 0.0
!!$  DO r=1,256
!!$     dt = 50.0/256
!!$     T2 = T2 + dt
!!$     CALL MULTIMODETIMEEVOLUTINOPERATOR(SIZE(U_F,1),SIZE(MODES_NUM,1),MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,FIELDS,T1,T2,U_AUX,INFO) 
!!$     dh2 = 2.0*eta*FIELDS(2)%OMEGA*(                 J_x*cos(fields(2)%omega*T2+phi1) + J_z*sin(fields(2)%omega*T2+phi1))*dt
!!$     dh3 = 2.0*eta*FIELDS(3)%OMEGA*(DCMPLX(0.0,-1.0)*J_y*cos(fields(3)%omega*T2+phi2) + J_z*sin(fields(3)%omega*T2+phi2))*dt
!!$     W2 = W2 + MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(dh2,U_AUX))
!!$     W3 = W3 + MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(dh3,U_AUX))     
!!$     WRITE(3,*) t2,ABS(U_AUX)**2,W2,W3
!!$  END DO
     WRITE(3,*)
     DEALLOCATE(E_FLOQUET)
     DEALLOCATE(U_F)
  END DO
  
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
