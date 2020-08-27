
PROGRAM MANYBODYHUBBARD

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE FLOQUETINITINTERFACE
  USE ARRAYS 


  IMPLICIT NONE
  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  TYPE(ATOM)                                       ID_F,ID_B
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_BARE
  INTEGER                                          INFO,m,INDEX0,r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET,E_BARE
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,H_BARE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  DOUBLE PRECISION                              :: T1,T2

  INTEGER index

  DOUBLE PRECISION t1_h,t2_h,mu



  OPEN(UNIT=3,FILE="ManyboydHubbard_Bosons.dat",ACTION="WRITE")
  OPEN(UNIT=3,FILE="ManyboydHubbard_Fermions.dat",ACTION="WRITE")



  INFO = 0
  N = 3
  L   = 2
  CALL FLOQUETINIT(ID_F,'ManybodyHubbard',N,L,'F',INFO)
  CALL FLOQUETINIT(ID_B,'ManybodyHubbard',N,L,'B',INFO)

  D_BARE = D_H(N_SITES,N_PARTICLES,'F') 
  ALLOCATE(E_BARE(D_BARE))
  ALLOCATE(H_BARE(D_BARE,D_BARE))
  ALLOCATE(U_AUX(D_BARE,D_BARE))
  ALLOCATE(MODES_NUM(2))
  ALLOCATE(states_F(D_BARE,N_P))
  ALLOCATE(states_B(D_BARE,N_P))

  CALL Manybody_basis(N,L,'F',states_F,INFO)
  CALL Manybody_basis(N,L,'B',states_B,INFO)

  CALL Tunneling(N,L,states_F,H_J,INFO)
  CALL Tunneling(N,L,states_B,H_J,INFO)

  CALL Onsite_twobody(N,L,states_F,H_U,INFO)
  CALL Onsite_twobody(N,L,states_B,H_U,INFO)


  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 1 !(DRIVING BY ONE HARMONIC)
  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(D_BARE,D_BARE))
     FIELDS(m)%V = 0.0
  END DO
  
  FIELDS(1)%V         = J*H_J + U*H_U
  FIELDS(1)%omega     = 0.0
  FIELDS(1)%N_Floquet = 0

  FIELDS(2)%V         = J_d*J
  FIELDS(2)%omega     = 2.0
  FIELDS(2)%N_Floquet = 10

  CALL LAPACK_FULLEIGENVALUES(FIELDS(1)%V,D_BARE,E_BARE,INFO)
  CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)

  ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
  ALLOCATE(U_F(SIZE(H_FLOQUET,1),SIZE(H_FLOQUET,1)))
  E_FLOQUET = 0.0   
  
  CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_FLOQUET,INFO)
  U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
  DEALLOCATE(H_FLOQUET)
  DO r=1,SIZE(U_F,1)
     WRITE(3,*) FIELDS(2)%omega,r,E_FLOQUET(r)
  END DO
  
 
END PROGRAM MANYBODYHUBBARD


SUBROUTINE Tunneling(N,L,state,J,INFO)

  INTEGER :: N,I_,J_,k_ 
  N = SIZE(state,1)
  allocate(new_state(L))

  J = 0
  DO k_=1,L ! loop though all sites
     DO J_=1,N
        DO I_=J_+1,N
           new_state = a_dagger(k_,state(J_))
           new_state = a_ (k_+1,new_state,K)
           J(I_,J_) = dot_product(state(I_,:),new_state)
        END DO
     END DO
     J = J + TRANSPOSE(CONJG(J))
     DO J_=1,N
        new_state = a_dagger(k_,state(J_))
        new_state = a_ (k_+1,new_state)
        J(J_,J_) = dot_product(state(J_,:),new_state)
     END DO
  END DO

END SUBROUTINE Tunneling

SUBROUTINE Onsite_twobody(N,L,state,J,INFO)

  INTEGER :: N,I_,J_,k_ 
  N = SIZE(state,1)
  allocate(new_state(L))

  J = 0
  DO k_=1,L ! loop though all sites
     DO J_=1,N
        DO I_=J_+1,N
           new_state = a_dagger(k_,state(J_))
           new_state = a_ (k_,new_state,K)
           J(I_,J_) = dot_product(state(I_,:),new_state)
        END DO
     END DO
     J = J + TRANSPOSE(CONJG(J))
     DO J_=1,N
        new_state = a_dagger(k_,state(J_))
        new_state = a_ (k_+1,new_state)
        J(J_,J_) = dot_product(state(J_,:),new_state)
     END DO
  END DO

END SUBROUTINE Onsite_twobody


SUBROUTINE Tunneling_F(Nup,Ndown,L,state,T,INFO)

  INTEGER :: N,I_,J_,k_ 
  N = SIZE(state,1)
  allocate(new_state(L))

  J = 0
  DO k_=1,L ! loop though all sites
     DO J_=1,N
        DO I_=J_+1,N
           new_state = a_dagger(k_,state(J_))
           new_state = a_ (k_+1,new_state,K)
           J(I_,J_) = dot_product(state(I_,:),new_state)
        END DO
     END DO
     J = J + TRANSPOSE(CONJG(J))
     DO J_=1,N
        new_state = a_dagger(k_,state(J_))
        new_state = a_ (k_+1,new_state)
        J(J_,J_) = dot_product(state(J_,:),new_state)
     END DO
  END DO

END SUBROUTINE Tunneling_F
