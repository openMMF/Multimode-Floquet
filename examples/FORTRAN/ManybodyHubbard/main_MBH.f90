
PROGRAM MANYBODYHUBBARD

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
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_J,H_U,U_F,H_BARE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  INTEGER         , DIMENSION(:,:), ALLOCATABLE :: states_occ
  DOUBLE PRECISION                              :: T1,T2

  INTEGER index,N_BODIES,N_SITES,i_

  DOUBLE PRECISION t1_h,t2_h,mu,d_h
  
  DOUBLE PRECISION t,u,t_driv,omega

  OPEN(UNIT=3,FILE="ManyboydHubbard_Bosons.dat",ACTION="WRITE")
  OPEN(UNIT=3,FILE="ManyboydHubbard_Fermions.dat",ACTION="WRITE")

  INFO = 0
  N_BODIES = 2
  N_SITES  = 7
  D_BARE = D_H(N_SITES,N_BODIES,'F') 
  !D_BARE = INT(D_H(N_SITES,N_BODIES,'B'))
  write(*,*) D_Bare
  CALL FLOQUETINIT(ID,'lattice',0.1d1*D_BARE,INFO)


  ALLOCATE(E_BARE(D_BARE))
  ALLOCATE(H_BARE(D_BARE,D_BARE))
  ALLOCATE(H_J(D_BARE,D_BARE))
  ALLOCATE(H_U(D_BARE,D_BARE))
  ALLOCATE(MODES_NUM(2))
  ALLOCATE(states_occ(D_BARE,N_SITES))

   ! CREATE THE BASIS OF STATES
  CALL Manybody_basis(D_BARE,N_SITES,N_BODIES,'F',states_occ,INFO)
  call write_matrix_int(states_occ)
!  CALL Manybody_basis(D_BARE,N_SITES,N_BODIES,'B',states_occ,INFO)
!  
!  ! EVALUATE THE TUNNELING TERM OF THE HAMILTONIAN
!  CALL Tunneling(D_BARE,N_SITES,N_BODIES,STATES_OCC,H_J,INFO)
!  !CALL WRITE_MATRIX(abs(h_j))
!  !! EVALUATE THE ON-SITE INTERACTION
!  CALL Onsite_twobody(D_BARE,N_SITES,N_BODIES,STATES_OCC,H_U,INFO)
!  !CALL WRITE_MATRIX(abs(h_u))
!
!  MODES_NUM(1) = 1 !(STATIC FIELD)
!  MODES_NUM(2) = 1 !(DRIVING BY ONE HARMONIC)
!  
!  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
!  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
!  DO m=1,TOTAL_FREQUENCIES
!     ALLOCATE(FIELDS(m)%V(D_BARE,D_BARE))
!     FIELDS(m)%V = 0.0
!  END DO
!  
!  
!  t      = 1.0
!  u      = 0.0
!  t_driv = 1.0
!  omega  = 1.0
!  FIELDS(1)%V         = -t*H_J + u*H_U
!  FIELDS(1)%omega     = 0.0
!  FIELDS(1)%N_Floquet = 0
!
!  FIELDS(2)%V         = t_driv*H_J
!  FIELDS(2)%omega     = omega
!  FIELDS(2)%N_Floquet = 2
!
!  H_BARE = FIELDS(1)%V
!  CALL LAPACK_FULLEIGENVALUES(H_BARE,D_BARE,E_BARE,INFO)
!  DO r=1,SIZE(E_BARE,1)
!     WRITE(*,*) FIELDS(2)%omega,r,E_BARE(r)
!  END DO
!  WRITE(*,*)
!  WRITE(*,*)
!  CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
!  !CALL WRITE_MATRIX(ABS(H_FLOUET()))
!
!  ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
!  ALLOCATE(U_F(SIZE(H_FLOQUET,1),SIZE(H_FLOQUET,1)))
!  E_FLOQUET = 0.0   
!  
!  CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_FLOQUET,INFO)
!  U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
!  DEALLOCATE(H_FLOQUET)
!  !WRITE(*,*) SIZE(U_F,1)
!  DO r=1,SIZE(U_F,1)
!     WRITE(*,*) FIELDS(2)%omega,r,E_FLOQUET(r)
!  END DO
   
END PROGRAM MANYBODYHUBBARD

SUBROUTINE Tunneling(D_BARE,N_SITES,N_BODIES,STATE,H_J,INFO)
   
  USE CREATIONDESTRUCTION
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: D_BARE,N_SITES,N_BODIES
  INTEGER,    DIMENSION(D_BARE,N_SITES), INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE,D_BARE),  INTENT(OUT)   :: H_J
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I
    
  INTEGER :: N,I_,J_,k_ 
  N = D_BARE
  
  H_J = 0
  J_=1
  K_=1
  DO k_=1,N_SITES-1 ! loop though all sites
     DO J_=1,N
        DO I_=J_+1,N
           !WRITE(*,*) K_,J_,I_,STATE(J_,:) 
            !WRITE(*,*) STATE(I_,:) 
            STATE_J = STATE(J_,:) 
            STATE_I = STATE(I_,:) 
            NEW_STATE = TUNNELING_(k_,STATE_J)
            !write(*,*) NEW_STATE,dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I))
            IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
                H_J(I_,J_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))
                H_J(J_,I_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))
            END IF           
        END DO
     END DO

     DO J_=1,N
         STATE_J = STATE(J_,:) 
         STATE_I = STATE(J_,:) 
         NEW_STATE = TUNNELING_(k_,STATE_J)
         IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
            H_J(J_,J_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))
         END IF           
     END DO
  END DO

END SUBROUTINE Tunneling

SUBROUTINE Onsite_twobody(D_BARE,N_SITES,N_BODIES,STATE,H_U,INFO)

    USE CREATIONDESTRUCTION
  
    IMPLICIT NONE    
    INTEGER,                               INTENT(IN)    :: D_BARE,N_SITES,N_BODIES
    INTEGER,    DIMENSION(D_BARE,N_SITES), INTENT(IN)    :: STATE
    COMPLEX*16, DIMENSION(D_BARE,D_BARE),  INTENT(OUT)   :: H_U
    INTEGER,                               INTENT(INOUT) :: INFO

    DOUBLE PRECISION, DIMENSION(N_SITES) :: NEW_STATE
    
    
    INTEGER :: N,I_,J_,k_ 
    N = D_BARE

    H_U = 0
    DO k_=1,N_SITES ! loop though all sites
        DO J_=1,N
            H_U(J_,J_) = STATE(J_,k_)
            H_U(J_,J_) = H_U(J_,J_)*(H_U(J_,J_)-1.0)
        END DO
    END DO

END SUBROUTINE Onsite_twobody

