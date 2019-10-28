SUBROUTINE MULTIMODEFLOQUETMATRIX(ATOM_C,NM,NF,MODES_NUM,FIELD_C,INFO)
  !ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO
  !  USE FLOQUET
  !ATOM_ type atom, -> dimension of the bare Hilbert space
  !NM -> number of modes
  !NF -> Number of Fields
  !MODES_NUM -> number of harmonics of each mode
  !FIELD -> Field couplings
  !INFO

  !USE ARRAYS
  !USE ATOMIC_PROPERTIES
  USE TYPES_C
  USE TYPES
  USE MODES_4F
  !USE SUBINTERFACE_LAPACK

  IMPLICIT NONE
  INTEGER,                     INTENT(IN)    :: NM,NF
  INTEGER,                     INTENT(INOUT) :: INFO
  INTEGER,      DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE_C), DIMENSION(NF), INTENT(IN)    :: FIELD_C
  TYPE(ATOM_C),                INTENT(IN)    :: ATOM_C                       

  TYPE(ATOM)                :: ATOM_
  !TYPE(MODE), DIMENSION(NM) :: FIELD

  CALL MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,COUPLING,INFO)

END SUBROUTINE MULTIMODEFLOQUETMATRIX


!!$SUBROUTINE MULTIMODETRANSITIONAVG(D,U_F_MODES,E_MULTIFLOQUET,D_BARE,U,INFO) 
SUBROUTINE MULTIMODETRANSITIONAVG(D,NM,FIELD,MODES_NUM,U_F_MODES,E_MULTIFLOQUET,D_BARE,U,INFO) 

  ! TIME EVOLUTION OPERATOR OF A MULTIMODE DRESSED SYSTEM. THE EVOLUTION OPERATOR IS WRITEN IN THE BASIS USED TO EXPRESS THE 
  ! MULTIMODE FLOQUET HAMILTONIAN
  ! U : MATRIX OF AMPLITUED OF PROBABILITIES FOR TRANSITIONS BETWEEN T1 TO T2
!ATOM_ type atom, -> dimension of the bare Hilbert space
!NM -> number of modes
!NF -> Number of Fields
!U_F_MODES -> TRANSFORMATINO BETWEEN THE BARE AND THE DRESSED BASIS
!MODES_NUM -> number of harmonics of each mode
!FIELD -> Field couplings
!INFO

!  USE FLOQUET
  USE TYPES
  USE SUBINTERFACE_LAPACK

  IMPLICIT NONE
!  INTEGER, INTENT(IN) :: NM,NF
!  INTEGER,                  INTENT(INOUT) :: INFO
!  INTEGER,   DIMENSION(NM), INTENT(IN)     :: MODES_NUM
  TYPE(MODE),DIMENSION(NM), INTENT(IN)     :: FIELD
  INTEGER,   DIMENSION(NM), INTENT(IN)     :: MODES_NUM
!  TYPE(ATOM),               INTENT(IN)    :: ATOM_                       

  INTEGER,                                    INTENT(IN)    :: D,D_BARE,NM ! DIMENSION OF THE MULTIMODE FLOQUET SPACE AND THE BARE BASIS
  INTEGER,                                    INTENT(INOUT) :: INFO
  DOUBLE PRECISION, DIMENSION(D),             INTENT(IN)    :: E_MULTIFLOQUET ! SET OF MULTIMODE FLOQUET ENERGIES, IN Hz, TO AVOID HBAR FACTORS
  COMPLEX*16,       DIMENSION(D,D),           INTENT(IN)    :: U_F_MODES   ! TRANFORMATION MATRIX BETWEEN DRESSED AND BARE BASIS
  DOUBLE PRECISION, DIMENSION(D_BARE,D_BARE), INTENT(OUT)   :: U           ! EVOLUTION OPERATOR U(T2,T1)

  COMPLEX*16, DIMENSION(D,D) :: U_DIAGONAL
  COMPLEX*16, DIMENSION(D_BARE,D) :: U_AUX

  INTEGER :: MULTIMODE_HARMONICS, n,i,j,m,index0,index1,field_index
  INTEGER, DIMENSION(NM) :: N_FLOQUET
  DOUBLE PRECISION, DIMENSION(NM) :: OMEGA_VEC
  TYPE(HARMONIC_FACTORS), DIMENSION(:),ALLOCATABLE:: U_MODES_n
  DOUBLE PRECISION :: PHASE

  MULTIMODE_HARMONICS   = D/D_BARE
  OMEGA_VEC  = 0.0
  U_DIAGONAL = 0.0
  U_AUX      = 0.0
  PHASE      = 0.0
  U          = 0.0

  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))

  DO n=1,MULTIMODE_HARMONICS
     ALLOCATE(U_MODES_n(n)%U(D_BARE,D))
     U_MODES_n(n)%U = 0.0
     ALLOCATE(U_MODES_n(n)%n(NM))
     U_MODES_n(n)%n = 0.0
  END DO

  DO i=0,MULTIMODE_HARMONICS-1
     U_MODES_n(i+1)%n = 0
     index0 = i
     DO j=2,NM
        U_MODES_n(i+1)%n(j)= -N_FLOQUET(j) + MOD(index0,(2*N_FLOQUET(j)+1))
        index0 = INT(index0/(2*N_FLOQUET(j)+1))
     END DO
     !     write(*,*) i+1,U_MODES_n(i+1)%n
  END DO

  i  = 1
  DO index1=1,MULTIMODE_HARMONICS
     DO n=1,D
        j = n
        i = 1 + (index1-1)*D_BARE
        DO m=1,D_BARE
           U_MODES_n(index1)%U(m,n) = U_F_MODES(i,j)
           i = i + 1
        END DO
     END DO
  END DO

  index0 = (MULTIMODE_HARMONICS -1)/2 + 1
  U_AUX  =  U_MODES_n(index0)%U
  U = 0.0
  DO index1=1,MULTIMODE_HARMONICS
     U = U +  MATMUL(ABS(U_AUX)**2,ABS(TRANSPOSE(CONJG(U_MODES_n(index1)%U)))**2)
  END DO

  info = 0

END SUBROUTINE MULTIMODETRANSITIONAVG


!!$SUBROUTINE MULTIMODETRANSITIONPROBABILITY_DRESSEDBASIS(D_F,D,U_F,U_F_MODES,E_DRESSED,E_MULTIFLOQUET, &
!!$     & D_BARE,P,P2,INFO)
!!$
!!$  !D                : dimension of the multimode hilbert space
!!$  !D_F              : dimension of the dressed hilbert space
!!$  !U_F_MODES        : D x D matrix transformation to the multidressed eigenstates
!!$  !U_F              : D_F x D_F matrix transformation to the single-dressed eigenstates
!!$  !E_DRESSED        : single mode eigenvalues == dressed energies
!!$  !E_MULTIFLOQUET   : multimode eigenvalues
!!$  !D_BARE           : dimension of the bare hilbert space
!!$  !P                : matrix of Time and phase average transition probabilities
!!$  !INFO             : for info
!!$
!!$  ! The dressed basis is not always well defined.
!!$  !Time average transition probability
!!$
!!$  USE SUBINTERFACE_LAPACK
!!$  USE FLOQUET
!!$  USE TYPES
!!$
!!$  IMPLICIT NONE
!!$  INTEGER,                                    INTENT(IN)     :: D_F,D,D_BARE
!!$  COMPLEX*16,       DIMENSION(D_F,D_F),       INTENT(INOUT)  :: U_F
!!$  COMPLEX*16,       DIMENSION(D,D),           INTENT(INOUT)  :: U_F_MODES
!!$  DOUBLE PRECISION, DIMENSION(D_F),           INTENT(IN)     :: E_DRESSED
!!$  DOUBLE PRECISION, DIMENSION(D),             INTENT(IN)     :: E_MULTIFLOQUET
!!$  !  COMPLEX*16,       DIMENSION(D_Bare,D_Bare),       INTENT(OUT)   :: P,P2
!!$  COMPLEX*16,       DIMENSION(D_F,D_F),       INTENT(OUT)    :: P,P2
!!$
!!$  INTEGER,                                    INTENT(INOUT)  :: INFO
!!$
!!$
!!$  TYPE(HARMONIC_FACTORS), DIMENSION(:),ALLOCATABLE :: U_n,U_MODES_n,U_nn,U_nn2
!!$  COMPLEX*16, DIMENSION(D,D)       :: U_AUX_F ! TEMPORAL VARIABLE, USED TO CHECK SUM_N U^n \DAGGER U^n = IDENTITY
!!$  !COMPLEX*16, DIMENSION(D_F,D_F)  :: U_AUX
!!$  COMPLEX*16, DIMENSION(D_F,D)     :: U_AUX
!!$  COMPLEX*16, DIMENSION(D,D_F)     :: U_AUX_2
!!$
!!$  DOUBLE PRECISION, DIMENSION(D)   :: PHASE_CONVENTION_MODES
!!$  DOUBLE PRECISION, DIMENSION(D_F) :: PHASE_CONVENTION_SINGLE
!!$
!!$  INTEGER n,m,index1,index2,i,j,index,index_,index3,index0
!!$  INTEGER DRESSED_HARMONICS, MULTIMODE_HARMONICS,D_AUX
!!$
!!$  INTEGER :: a,c,e,g,k,k_,n_,lambda,l,l_
!!$
!!$  DRESSED_HARMONICS   = D_F/D_BARE
!!$  MULTIMODE_HARMONICS =   D/D_BARE
!!$
!!$  ALLOCATE(U_n(DRESSED_HARMONICS))
!!$  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))
!!$  ALLOCATE(U_nn(DRESSED_HARMONICS*MULTIMODE_HARMONICS))
!!$  ALLOCATE(U_nn2(DRESSED_HARMONICS*MULTIMODE_HARMONICS))
!!$
!!$  DO n=1,DRESSED_HARMONICS
!!$     ALLOCATE(U_n(n)%U(D_BARE,D_F))   ! All singlemode dressed states
!!$     !ALLOCATE(U_n(n)%U(D_BARE,D_BARE)) ! Select a manifold of singlemode dressed states
!!$     U_n(n)%U = 0.0
!!$  END DO
!!$
!!$  DO n=1,MULTIMODE_HARMONICS
!!$     ALLOCATE(U_MODES_n(n)%U(D_BARE,D)) !Full set of multimode dressed states
!!$     !ALLOCATE(U_MODES_n(n)%U(D_BARE,D_BARE)) !Single out a multimode dressed manifold
!!$     U_MODES_n(n)%U = 0.0
!!$     ALLOCATE(U_MODES_n(n)%n(NM))
!!$     U_MODES_n(n)%n = 0
!!$  END DO
!!$
!!$  DO n=1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$     !ALLOCATE(U_nn(n)%U(D_F,D_BARE))
!!$     !U_nn(n)%U = 0.0
!!$     !ALLOCATE(U_nn2(n)%U(D_F,D_BARE))
!!$     !U_nn2(n)%U = 0.0
!!$     ALLOCATE(U_nn(n)%U(D_F,D))
!!$     U_nn(n)%U = 0.0
!!$     ALLOCATE(U_nn2(n)%U(D_F,D))
!!$     U_nn2(n)%U = 0.0
!!$  END DO
!!$  i  = 1
!!$  DO index1=1,DRESSED_HARMONICS
!!$     i = 1 + (index1-1)*D_BARE
!!$     U_n(index1)%U = U_F(i:i+D_BARE-1,1:D_F)
!!$     !index2 = D_BARE*N_FLOQUET(2)+1
!!$     !U_n(index1)%U = U_F(i:i+D_BARE-1,index2:index2+D_BARE-1)
!!$     !call write_matrix(ABS(U_n(index1)%U))
!!$     !
!!$  END DO
!!$
!!$
!!$  DO i=0,MULTIMODE_HARMONICS-1
!!$     U_MODES_n(i+1)%n = 0
!!$     index1 = i
!!$     DO j=2,MODES_NUM
!!$        U_MODES_n(i+1)%n(j)= -N_FLOQUET(j) + MOD(index1,(2*N_FLOQUET(j)+1))
!!$        index1 = INT(index1/(2*N_FLOQUET(j)+1))
!!$     END DO
!!$  END DO
!!$
!!$  i  = 1
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     i = 1 + (index1-1)*D_BARE
!!$     U_MODES_n(index1)%U = U_F_MODES(i:i+D_BARE-1,1:D)
!!$     !WRITE(*,*) INDEX1,U_MODES_n(index1)%n
!!$     !call write_matrix(ABS(U_MODES_n(index1)%U))
!!$
!!$     !index3 = MULTIMODE_HARMONICS*D_BARE/2.0 -1
!!$     !     !U_MODES_n(index1)%U = U_F_MODES(i:i+D_BARE-1,index3:index3+1)
!!$  END DO
!!$
!!$  index1 =1
!!$  DO n=1,DRESSED_HARMONICS
!!$     DO m=1,MULTIMODE_HARMONICS
!!$        U_nn(index1)%U  = ABS(MATMUL(TRANSPOSE(CONJG(U_n(n)%U)),U_MODES_n(m)%U))**2
!!$        !      CALL WRITE_MATRIX(ABS(U_n(n)%U))
!!$        !      CALL WRITE_MATRIX(ABS(U_MODES_n(m)%U))
!!$        !      CALL WRITE_MATRIX(ABS(U_nn(index1)%U))
!!$        index1 =  index1 + 1
!!$     END DO
!!$  END DO
!!$
!!$  !P = ABS(MATMUL(TRANSPOSE(CONJG(U_F)),U_F_MODES))**2
!!$  !P = MATMUL(P,TRANSPOSE(CONJG(P)))
!!$  !
!!$  INDEX0 = (DRESSED_HARMONICS -1)/2 + 1
!!$  INDEX1 = (MULTIMODE_HARMONICS -1)/2 + 1
!!$  !write(*,*) '#',index0,index1
!!$
!!$  U_AUX = ABS(MATMUL(TRANSPOSE(CONJG(U_n(INDEX0)%U)),U_MODES_n(INDEX1)%U))**2
!!$
!!$  !as in 10.10.17
!!$  P = 0.0
!!$  DO index1  = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS!,MULTIMODE_HARMONICS!
!!$     !write(*,*) size(U_aux,1),size(U_aux,2),size(u_nn(index1)%U,1),size(u_nn(index1)%U,2)
!!$     !P = P + MATMUL(U_AUX,transpose(conjg(U_nn(index1)%U)))
!!$     P = P + MATMUL(U_nn(index1)%U,TRANSPOSE(U_AUX))
!!$     !P = P + MATMUL(U_AUX,TRANSPOSE(U_AUX))
!!$  END DO
!!$
!!$
!!$  !  write(*,*) real(P(4,4)),real(P(4,5)),real(P(5,5))
!!$  !  P = 0.0
!!$  !  DO index1  = MULTIMODE_HARMONICS+1,MULTIMODE_HARMONICS*2!DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !     P = P + MATMUL(U_nn(index1)%U,TRANSPOSE(U_AUX))
!!$  !  END DO
!!$  !  !CALL WRITE_MATRIX(ABS(P))
!!$  !
!!$  !  P = 0.0
!!$  !  DO index1  = 2*MULTIMODE_HARMONICS+1,MULTIMODE_HARMONICS*3!DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !     P = P + MATMUL(U_nn(index1)%U,TRANSPOSE(U_AUX))
!!$  !  END DO
!!$  !  !CALL WRITE_MATRIX(ABS(P))
!!$  !
!!$  !
!!$  !  P = 0.0
!!$  !  DO index1  = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !     P = P + MATMUL(U_nn(index1)%U,TRANSPOSE(U_AUX))
!!$  !  END DO
!!$  !  !CALL WRITE_MATRIX(ABS(P))
!!$  !
!!$  !  P = 0.0
!!$  !  P = MATMUL(U_AUX,TRANSPOSE(U_AUX))
!!$  !  !CALL WRITE_MATRIX(ABS(P))
!!$
!!$
!!$  ! as in 19.01.17 b.
!!$
!!$  !DO index1 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !   U_nn2(index1)%U = abs(U_nn(index1)%U)**2
!!$  !END DO
!!$  !P = 0.0
!!$  !DO index1 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !   DO index2 = 1,DRESSED_HARMONICS
!!$  !      INDEX3 = INDEX2 + (MULTIMODE_HARMONICS-1)/2
!!$  !      P =  P + MATMUL(U_nn2(index1)%U,TRANSPOSE(CONJG(U_nn2(index3)%U)))
!!$  !   END DO
!!$  !END DO
!!$
!!$
!!$  !! as in 17.01.17 ! selecting a dressed manifold
!!$
!!$  !  ALLOCATE(U_n(DRESSED_HARMONICS))
!!$  !  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))
!!$  !  ALLOCATE(U_nn(DRESSED_HARMONICS*MULTIMODE_HARMONICS))
!!$  !  ALLOCATE(U_nn2(DRESSED_HARMONICS*MULTIMODE_HARMONICS))
!!$  !
!!$  !  DO n=1,DRESSED_HARMONICS
!!$  !     ALLOCATE(U_n(n)%U(D_BARE,D_BARE))
!!$  !     U_n(n)%U = 0.0
!!$  !  END DO
!!$  !
!!$  !  D_AUX = D_BARE
!!$  !  DO n=1,MULTIMODE_HARMONICS
!!$  !     ALLOCATE(U_MODES_n(n)%U(D_BARE,D_aux))
!!$  !     U_MODES_n(n)%U = 0.0
!!$  !     ALLOCATE(U_MODES_n(n)%n(MODES_NUM))
!!$  !     U_MODES_n(n)%n = 0
!!$  !  END DO
!!$  !
!!$  !  DO n=1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !     ALLOCATE(U_nn(n)%U(D_BARE,D_AUX))
!!$  !     U_nn(n)%U = 0.0
!!$  !     ALLOCATE(U_nn2(n)%U(D_BARE,D_AUX))
!!$  !     U_nn2(n)%U = 0.0
!!$  !  END DO
!!$  !
!!$  !  DO i=0,MULTIMODE_HARMONICS-1
!!$  !     U_MODES_n(i+1)%n = 0
!!$  !     index1 = i
!!$  !     DO j=2,MODES_NUM
!!$  !        U_MODES_n(i+1)%n(j)= -N_FLOQUET(j) + MOD(index1,(2*N_FLOQUET(j)+1))
!!$  !        index1 = INT(index1/(2*N_FLOQUET(j)+1))
!!$  !     END DO
!!$  !  END DO
!!$  !
!!$  !
!!$  !  i  = 1
!!$  !  DO index1=1,DRESSED_HARMONICS
!!$  !     i =1 + (index1-1)*D_BARE
!!$  !     index2 = D_BARE*N_FLOQUET(2)+2
!!$  !     U_n(index1)%U = U_F(i:i+D_BARE-1,index2:index2+D_BARE-1)
!!$  !  END DO
!!$  !  i  = 1
!!$  !  DO index1=1,MULTIMODE_HARMONICS
!!$  !     i = 1 + (index1-1)*D_BARE
!!$  !     !U_MODES_n(index1)%U = U_F_MODES(i:i+D_BARE-1,1:D)
!!$  !     index = MULTIMODE_HARMONICS*D_BARE/2.0 + D_BARE+1
!!$  !     U_MODES_n(index1)%U = U_F_MODES(i:i+D_BARE-1,index:index+D_AUX-1)
!!$  !     if(index1.eq.5) then
!!$  !         P = U_MODES_n(index1)%U
!!$  !     end if
!!$  !  END DO
!!$  !
!!$  !  index1 =1
!!$  !  DO n=1,DRESSED_HARMONICS
!!$  !     DO m=1,MULTIMODE_HARMONICS
!!$  !        U_nn(index1)%U  = MATMUL(TRANSPOSE(CONJG(U_n(n)%U)),U_MODES_n(m)%U)
!!$  !        !write(*,*) n,m,index1
!!$  !        index1 =  index1 + 1
!!$  !     END DO
!!$  !  END DO
!!$  !
!!$  !!17.01.17  expressions where the sums are over a single set of multimode dressed states, page 1
!!$  !
!!$  !  DO index1 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !     U_nn2(index1)%U = ABS(U_nn(index1)%U)**2
!!$  !!    CALL WRITE_MATRIX(ABS(U_nn2(index1)%U))
!!$  !  END DO
!!$  !
!!$  !  P = 0.0
!!$  !  DO index1 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !     DO index2 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !              P      = P+MATMUL(U_nn2(INDEX1)%U,TRANSPOSE(U_nn2(INDEX2)%U))
!!$  ! !       CALL WRITE_MATRIX(REAL(P))
!!$  !     END DO
!!$  !  END DO
!!$
!!$  !13.01.17  expressions where the sums are over a single set of multimode dressed states, page 12
!!$
!!$  !  DO index1 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !     U_nn2(index1)%U = MATMUL(U_nn(index1)%U,TRANSPOSE(CONJG(U_nn(index1)%U)))
!!$  !  END DO
!!$  !
!!$  !  P = 0.0
!!$  !  aDO i=1,D_BARE
!!$  !     DO j=1,D_BARE
!!$  !        DO index1 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !           DO index2 = 1,DRESSED_HARMONICS*MULTIMODE_HARMONICS
!!$  !              P(i,j) = P(i,j) + U_nn2(index1)%U(j,j)*U_nn2(index2)%U(i,i)
!!$  !           END DO
!!$  !        END DO
!!$  !     END DO
!!$  !  END DO
!!$  !
!!$
!!$
!!$END SUBROUTINE MULTIMODETRANSITIONPROBABILITY_DRESSEDBASIS
!!$
!!$
!!$
!!$SUBROUTINE MULTIMODEFLOQUETTRANSFORMATION(D,NF,U_F_MODES,E_MULTIFLOQUET,D_BARE,FIELD,T1,U,INFO) 
!!$
!!$
!!$
!!$  ! TIME-DEPENDENT TRANSFORMATION BETWEEN THE BARE BASIS AND THE FLOQUET STATES
!!$  ! U(T1) = sum_ U^n exp(i n omega T1)
!!$  ! 
!!$
!!$  USE TYPES
!!$
!!$  IMPLICIT NONE
!!$  INTEGER,                                    INTENT(IN)    :: D,D_BARE,NF ! DIMENSION OF THE MULTIMODE FLOQUET SPACE AND THE BARE BASIS
!!$  INTEGER,                                    INTENT(INOUT) :: INFO
!!$  TYPE(MODE),       DIMENSION(NF),            INTENT(IN)    :: FIELD  ! FIELDS PROPERTIES: FREQUENCY, AMPLITUDES AND PHASES
!!$  DOUBLE PRECISION,                           INTENT(IN)    :: T1 ! IN SECONDS
!!$  DOUBLE PRECISION, DIMENSION(D),             INTENT(IN)    :: E_MULTIFLOQUET ! SET OF MULTIMODE FLOQUET ENERGIES, IN Hz, TO AVOID HBAR FACTORS
!!$  COMPLEX*16,       DIMENSION(D,D),           INTENT(IN)    :: U_F_MODES   ! TRANFORMATION MATRIX BETWEEN DRESSED FLOQUET AND BARE EXTENDED BASIS
!!$  COMPLEX*16,       DIMENSION(D_BARE,D), INTENT(OUT)   :: U           ! TIME-DEPENDENT TRANSFORMATINO BETWEEN THE DRESSED AND BARE BASIS
!!$
!!$
!!$  COMPLEX*16,       DIMENSION(D,D)      :: U_DIAGONAL
!!$  COMPLEX*16,       DIMENSION(D_BARE,D) :: U_AUX
!!$
!!$  INTEGER                                                :: MULTIMODE_HARMONICS, n,i,j,m,index0,index1
!!$  DOUBLE PRECISION,       DIMENSION(NF)                  :: OMEGA_VEC
!!$  TYPE(HARMONIC_FACTORS), DIMENSION(:),      ALLOCATABLE :: U_MODES_n
!!$  DOUBLE PRECISION                                       :: PHASE
!!$
!!$  MULTIMODE_HARMONICS   = D/D_BARE
!!$  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))
!!$  DO n=1,MULTIMODE_HARMONICS
!!$     ALLOCATE(U_MODES_n(n)%U(D_BARE,D))
!!$     U_MODES_n(n)%U = 0.0
!!$     ALLOCATE(U_MODES_n(n)%n(NF))
!!$     U_MODES_n(n)%n = 0
!!$  END DO
!!$  DO i=1,NF
!!$     OMEGA_VEC(i) = FIELD(i)%OMEGA
!!$  END DO
!!$  i  = 1
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     DO n=1,D
!!$        j = n
!!$        i = 1 + (index1-1)*D_BARE 
!!$        DO m=1,D_BARE
!!$           U_MODES_n(index1)%U(m,n) = U_F_MODES(i,j)
!!$           i = i + 1           
!!$        END DO
!!$     END DO
!!$  END DO
!!$
!!$  DO i=0,MULTIMODE_HARMONICS-1
!!$     U_MODES_n(i+1)%n = 0
!!$     index0 = i
!!$     DO j=2,NF
!!$        U_MODES_n(i+1)%n(j)= -FIELD(j)%N_FLOQUET + MOD(index0,(2*FIELD(j)%N_FLOQUET+1))
!!$        index0 = INT(index0/(2*FIELD(j)%N_FLOQUET+1))
!!$     END DO
!!$  END DO
!!$
!!$
!!$  U = 0.0
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     PHASE = DOT_PRODUCT(U_MODES_n(index1)%n,omega_vec)*T1
!!$     U      = U  + EXP(DCMPLX(0.0,1.0)*PHASE)*U_MODES_n(index1)%U
!!$  END DO
!!$
!!$
!!$END SUBROUTINE MULTIMODEFLOQUETTRANSFORMATION
!!$
!!$SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR(D,NF,U_F_MODES,E_MULTIFLOQUET,D_BARE,FIELD,T1,T2,U,INFO) 
!!$
!!$  ! TIME EVOLUTION OPERATOR OF A MULTIMODE DRESSED SYSTEM. THE EVOLUTION OPERATOR IS WRITEN IN THE BASIS USED TO EXPRESS THE 
!!$  ! MULTIMODE FLOQUET HAMILTONIAN
!!$  ! U : MATRIX OF AMPLITUED OF PROBABILITIES FOR TRANSITIONS BETWEEN T1 TO T2
!!$
!!$  USE TYPES
!!$  USE SUBINTERFACE_LAPACK
!!$
!!$
!!$  IMPLICIT NONE
!!$  INTEGER,                                    INTENT(IN)    :: D,D_BARE,NF ! DIMENSION OF THE MULTIMODE FLOQUET SPACE AND THE BARE BASIS
!!$  INTEGER,                                    INTENT(INOUT) :: INFO
!!$  TYPE(MODE),       DIMENSION(NF),            INTENT(IN)    :: FIELD  ! FIELDS PROPERTIES: FREQUENCY, AMPLITUDES AND PHASES
!!$  DOUBLE PRECISION,                           INTENT(IN)    :: T1,T2  ! IN SECONDS
!!$  DOUBLE PRECISION, DIMENSION(D),             INTENT(IN)    :: E_MULTIFLOQUET ! SET OF MULTIMODE FLOQUET ENERGIES, IN Hz, TO AVOID HBAR FACTORS
!!$  COMPLEX*16,       DIMENSION(D,D),           INTENT(IN)    :: U_F_MODES   ! TRANFORMATION MATRIX BETWEEN DRESSED AND BARE BASIS
!!$  COMPLEX*16,       DIMENSION(D_BARE,D_BARE), INTENT(OUT)   :: U           ! EVOLUTION OPERATOR U(T2,T1)
!!$
!!$  COMPLEX*16, DIMENSION(D,D) :: U_DIAGONAL
!!$  COMPLEX*16, DIMENSION(D_BARE,D) :: U_AUX
!!$
!!$  INTEGER :: MULTIMODE_HARMONICS, n,i,j,m,index0,index1
!!$  DOUBLE PRECISION, DIMENSION(NF) :: OMEGA_VEC
!!$  TYPE(HARMONIC_FACTORS), DIMENSION(:),ALLOCATABLE:: U_MODES_n
!!$  DOUBLE PRECISION :: PHASE
!!$
!!$  MULTIMODE_HARMONICS   = D/D_BARE
!!$
!!$  DO i=1,NF
!!$     OMEGA_VEC(i) = FIELD(i)%OMEGA
!!$  END DO
!!$
!!$
!!$  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))
!!$
!!$  DO n=1,MULTIMODE_HARMONICS
!!$     ALLOCATE(U_MODES_n(n)%U(D_BARE,D))
!!$     U_MODES_n(n)%U = 0.0
!!$     ALLOCATE(U_MODES_n(n)%n(NF))
!!$     U_MODES_n(n)%n = 0.0
!!$  END DO
!!$
!!$  DO i=0,MULTIMODE_HARMONICS-1
!!$     U_MODES_n(i+1)%n = 0
!!$     index0 = i
!!$     DO j=2,NF
!!$        U_MODES_n(i+1)%n(j)= -FIELD(j)%N_FLOQUET + MOD(index0,(2*FIELD(j)%N_FLOQUET+1))
!!$        index0 = INT(index0/(2*FIELD(j)%N_FLOQUET+1))        
!!$     END DO
!!$  END DO
!!$
!!$  i  = 1
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     DO n=1,D
!!$        j = n
!!$        i = 1 + (index1-1)*D_BARE 
!!$        DO m=1,D_BARE
!!$           U_MODES_n(index1)%U(m,n) = U_F_MODES(i,j)
!!$           i = i + 1           
!!$        END DO
!!$     END DO
!!$  END DO
!!$
!!$  index0 = (MULTIMODE_HARMONICS -1)/2 + 1
!!$  U_AUX  =  U_MODES_n(index0)%U
!!$
!!$
!!$  U_DIAGONAL = 0.0
!!$  DO i=1,D
!!$     U_DIAGONAL(i,i) = EXP(-DCMPLX(0.0,1.0)*E_MULTIFLOQUET(i)*(T2-T1))
!!$  END DO
!!$
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     U_MODES_n(index1)%U = MATMUL(U_MODES_n(index1)%U,U_DIAGONAL)
!!$  END DO
!!$
!!$  U = 0.0
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     PHASE = DOT_PRODUCT(U_MODES_n(index1)%n,omega_vec)*T2
!!$     U = U +  MATMUL(U_MODES_n(index1)%U,TRANSPOSE(CONJG(U_AUX)))*EXP(DCMPLX(0.0,1.0)*PHASE)
!!$  END DO
!!$
!!$END SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR
!!$
!!$SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR_RESTRICTED(D,NF,U_F_MODES,E_MULTIFLOQUET,D_BARE,FIELD,T1,T2,U,INFO) 
!!$
!!$  ! TIME EVOLUTION OPERATOR OF A MULTIMODE DRESSED SYSTEM. THE EVOLUTION OPERATOR IS WRITEN IN THE BASIS USED TO EXPRESS THE 
!!$  ! MULTIMODE FLOQUET HAMILTONIAN
!!$  ! U : MATRIX OF AMPLITUED OF PROBABILITIES FOR TRANSITIONS BETWEEN T1 TO T2
!!$  !
!!$  ! BUT THE SUM OVER FLOQUET EIGENVALUES IS RESTRICTED TO THOSE ACCURATELY DEFINED, FAR FROM THE EDGES OF THE MATRICES  
!!$  ! THE SUBROUTINE HAS BEEN DESIGNED FOR 87Rb driven by RF+MW fields, AND IT IS NOT EXPECTED TO WORK FOR OTHER CONFIGURATIONS 
!!$  ! D: THE DIMENSOIN OF THE MULTIMODE FLOQUET IS NOW DEFINED AS THE NUMBER OF FLOQUET MODES USED TO CALCULATE THE EVOLUTION OPERATOR
!!$  ! the point of this is to make a comparison with MULTIMODETIMEEVOLUTIONOPERATOR AND CHECK WHETHER.
!!$  ! MULTIMODETIMEEVOLUTIONOPERATOR SHOULD WORK FOR ANY ATOM/FIELD CONFIGURATION.
!!$  ! MULTIMODETIMEEVOLUTIONOPEATOR_RESTRICTED IS DIFFICULT TO RESTRICT IN GENERAL SINCE WE NEED TO IDENTIFY THE "CENTRAL" MANIFOLD OF FLOQUET 
!!$  ! STATES
!!$
!!$  USE TYPES
!!$  USE SUBINTERFACE_LAPACK
!!$
!!$
!!$  IMPLICIT NONE
!!$  INTEGER,                                    INTENT(IN)    :: D,D_BARE,NF ! DIMENSION OF THE MULTIMODE FLOQUET SPACE AND THE BARE BASIS
!!$  INTEGER,                                    INTENT(INOUT) :: INFO
!!$  TYPE(MODE),       DIMENSION(NF),     INTENT(IN)    :: FIELD  ! FIELDS PROPERTIES: FREQUENCY, AMPLITUDES AND PHASES
!!$  DOUBLE PRECISION,                           INTENT(IN)    :: T1,T2  ! IN SECONDS
!!$  DOUBLE PRECISION, DIMENSION(D),             INTENT(IN)    :: E_MULTIFLOQUET ! SET OF MULTIMODE FLOQUET ENERGIES, IN Hz, TO AVOID HBAR FACTORS
!!$  COMPLEX*16,       DIMENSION(D,D),           INTENT(IN)    :: U_F_MODES   ! TRANFORMATION MATRIX BETWEEN DRESSED AND BARE BASIS
!!$  COMPLEX*16,       DIMENSION(D_BARE,D_BARE), INTENT(OUT)   :: U           ! EVOLUTION OPERATOR U(T2,T1)
!!$
!!$  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: U_DIAGONAL
!!$  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: U_AUX
!!$
!!$  INTEGER :: MULTIMODE_HARMONICS, n,i,j,m,index0,index1,D_r,i_0,i_r,Fup,Fdown
!!$  DOUBLE PRECISION, DIMENSION(NF) :: OMEGA_VEC
!!$  TYPE(HARMONIC_FACTORS), DIMENSION(:),ALLOCATABLE:: U_MODES_n
!!$  DOUBLE PRECISION :: PHASE
!!$  DOUBLE PRECISION :: A,HBAR
!!$
!!$  MULTIMODE_HARMONICS   = D/D_BARE
!!$  IF(NF.GT.2) THEN
!!$     D_r                   = (2*FIELD(3)%N_FLOQUET - 2)*(2*(FIELD(2)%N_FLOQUET-2)-1)*D_bare
!!$  END IF
!!$  Fup                   = 2
!!$  Fdown                 = 1
!!$  hbar                  = 1.054E-34 
!!$  A                     = 2.0*4.0*ATAN(1.0)*hbar*3.417341E9
!!$
!!$  DO i=1,NF
!!$     OMEGA_VEC(i) = FIELD(i)%OMEGA
!!$  END DO
!!$
!!$  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))
!!$
!!$  DO n=1,MULTIMODE_HARMONICS
!!$     ALLOCATE(U_MODES_n(n)%U(D_BARE,D))
!!$     U_MODES_n(n)%U = 0.0
!!$     ALLOCATE(U_MODES_n(n)%U_r(D_BARE,D_r))
!!$     U_MODES_n(n)%U_r = 0.0
!!$     ALLOCATE(U_MODES_n(n)%n(NF))
!!$     U_MODES_n(n)%n = 0.0
!!$  END DO
!!$
!!$  DO i=0,MULTIMODE_HARMONICS-1
!!$     U_MODES_n(i+1)%n = 0
!!$     index0 = i
!!$     DO j=2,NF
!!$        U_MODES_n(i+1)%n(j)= -FIELD(j)%N_FLOQUET + MOD(index0,(2*FIELD(j)%N_FLOQUET+1))
!!$        index0 = INT(index0/(2*FIELD(j)%N_FLOQUET+1))
!!$     END DO
!!$  END DO
!!$
!!$  i  = 1
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     DO n=1,D
!!$        j = n
!!$        i = 1 + (index1-1)*D_BARE 
!!$        DO m=1,D_BARE
!!$           U_MODES_n(index1)%U(m,n) = U_F_MODES(i,j)
!!$           i = i + 1           
!!$        END DO
!!$     END DO
!!$  END DO
!!$
!!$  index0 = (MULTIMODE_HARMONICS -1)/2 + 1
!!$
!!$
!!$  !HOW MANY FLOQUET MODES AND WHICH ONES ARE WE GOING TO USE?
!!$  !
!!$  ! D_r  = (2*FIELD(3)%N_FLOQUET - 2)*(2*(FIELD(2)%N_FLOQUET-2)-1)*RB87%D_bare
!!$
!!$  ! (2*FIELD(3)%N_FLOQUET - 2)       : number of mw manifolds, we neglect one up and one down 
!!$  ! (2*(FIELD(2)%N_FLOQUET - 2) - 1) : number of complete (i.e. having Rb87%D_bare) RF manifolds (not in general, it depends on Bdc and omega_rf)
!!$  ! Rb87%D_bare                : number of bare states
!!$
!!$  !WITHIN EACH MW MANIFOLD, WE KEEP FLOQUETS MODES WITH INDICES 
!!$  ! I \IN  I_0+1, I_0 +  (2*(FIELD(2)%N_FLOQUET - 2) - 1)*Rb87%D_bare    
!!$  ! WITH 
!!$  ! I_0 = 3*(2*Fdown+2*Fup+2)+1
!!$  ! and the starting index for each manifold is
!!$  ! (2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+1)+(m)*(2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+2*Fup+2)
!!$  ! with m the manifold number, but we keep manifolds m=1,(2*FIELD(3)%N_FLOQUET - 2) 
!!$
!!$
!!$  ALLOCATE(U_DIAGONAL(D_r,D_r))
!!$  ALLOCATE(U_AUX(D_bare,D_r))
!!$
!!$
!!$  n  = 1
!!$  U_DIAGONAL = 0.0
!!$  DO i=1,(2*FIELD(3)%N_FLOQUET - 2)
!!$     I_0 = (2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+1) + i*(2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+2*Fup+2)
!!$     DO j=1,(2*(FIELD(2)%N_FLOQUET-2)-1)*D_bare
!!$
!!$        i_r = I_0 + 3*(2*Fdown+2*Fup+2)+1 + j 
!!$        U_DIAGONAL(n,n) = EXP(-DCMPLX(0.0,1.0)*E_MULTIFLOQUET(i_r)*(T2-T1))
!!$        n = n + 1
!!$     END DO
!!$  END DO
!!$
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     n = 1
!!$     DO i=1,(2*FIELD(3)%N_FLOQUET - 2)
!!$        I_0 = (2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+1) + i*(2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+2*Fup+2)
!!$        DO j=1,(2*(FIELD(2)%N_FLOQUET-2)-1)*D_bare        
!!$           i_r = I_0 + 3*(2*Fdown+2*Fup+2)+1 + j 
!!$           U_MODES_n(index1)%U_r(:,n) = U_MODES_n(index1)%U(:,i_r)
!!$           n = n+1
!!$        END DO
!!$     END DO
!!$  END DO
!!$
!!$  index0 = (MULTIMODE_HARMONICS -1)/2 + 1
!!$  U_AUX  =  U_MODES_n(index0)%U_r
!!$
!!$
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     U_MODES_n(index1)%U_r = MATMUL(U_MODES_n(index1)%U_r,U_DIAGONAL)
!!$  END DO
!!$
!!$  U = 0.0
!!$  DO index1=1,MULTIMODE_HARMONICS
!!$     PHASE = DOT_PRODUCT(U_MODES_n(index1)%n,omega_vec)*T2
!!$     U = U +  MATMUL(U_MODES_n(index1)%U_r,TRANSPOSE(CONJG(U_AUX)))*EXP(DCMPLX(0.0,1.0)*PHASE)
!!$  END DO
!!$
!!$
!!$END SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR_RESTRICTED
!!$



!!$SUBROUTINE COORDINATEPACKING(D,A,V,R,C,index,INFO)
!!$  IMPLICIT NONE
!!$  INTEGER,INTENT(IN):: D
!!$  COMPLEX*16,DIMENSION(D,D),INTENT(IN)  :: A
!!$  COMPLEX*16,DIMENSION(D*D),INTENT(OUT) :: V
!!$  INTEGER, DIMENSION(D*D),  INTENT(OUT) :: R,C
!!$  INTEGER, INTENT(OUT)   :: index
!!$  INTEGER, INTENT(INOUT) :: INFO
!!$  
!!$  INTEGER I,J
!!$  V=0
!!$  R=0
!!$  C=0
!!$  
!!$  index = 1
!!$  DO I=1,D
!!$     DO J=1,D
!!$        IF(ABS(A(I,J)).GT.0) THEN
!!$           V(index) = A(I,J)
!!$           R(index) = I
!!$           C(index) = J
!!$           index = index+1
!!$        END IF
!!$     END DO
!!$  END DO
!!$  index = index-1
!!$END SUBROUTINE COORDINATEPACKING
!!$
!!$MODULE MERGINGARRAYS
!!$  INTERFACE
!!$     SUBROUTINE APPENDARRAYS(V,B,INFO)
!!$       COMPLEX*16, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
!!$       COMPLEX*16, DIMENSION(:),INTENT(IN)    :: B
!!$       INTEGER,                 INTENT(INOUT) :: INFO
!!$     END SUBROUTINE APPENDARRAYS
!!$     SUBROUTINE APPENDARRAYSI(V,B,INFO)
!!$       INTEGER, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
!!$       INTEGER, DIMENSION(:),INTENT(IN)    :: B
!!$       INTEGER,                 INTENT(INOUT) :: INFO
!!$     END SUBROUTINE APPENDARRAYSI
!!$  END INTERFACE
!!$END MODULE MERGINGARRAYS
!!$
!!$SUBROUTINE APPENDARRAYS(V,B,INFO)
!!$  COMPLEX*16, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
!!$  COMPLEX*16, DIMENSION(:),INTENT(IN)    :: B
!!$  INTEGER,                 INTENT(INOUT) :: INFO
!!$  
!!$  COMPLEX*16,DIMENSION(:),ALLOCATABLE :: tmp_arr
!!$!  write(*,*) V
!!$!  write(*,*) B
!!$  ALLOCATE(tmp_arr(SIZE(V,1)+SIZE(B,1)))
!!$  tmp_arr(1:SIZE(V,1))=V
!!$  tmp_arr(SIZE(V,1)+1:SIZE(tmp_arr))=B
!!$  DEALLOCATE(V)
!!$  ALLOCATE(V(SIZE(tmp_arr)))
!!$  V=tmp_arr
!!$END SUBROUTINE APPENDARRAYS
!!$
!!$SUBROUTINE APPENDARRAYSI(V,B,INFO)
!!$  INTEGER, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
!!$  INTEGER, DIMENSION(:),INTENT(IN)    :: B
!!$  INTEGER,                 INTENT(INOUT) :: INFO
!!$  
!!$  COMPLEX*16,DIMENSION(:),ALLOCATABLE :: tmp_arr
!!$  
!!$  ALLOCATE(tmp_arr(SIZE(V,1)+SIZE(B,1)))
!!$  tmp_arr(1:SIZE(V,1))=V
!!$  tmp_arr(SIZE(V,1)+1:SIZE(tmp_arr))=B
!!$  DEALLOCATE(V)
!!$  ALLOCATE(V(SIZE(tmp_arr)))
!!$  V=tmp_arr
!!$END SUBROUTINE APPENDARRAYSI

!!$SUBROUTINE MULTIMODEFLOQUETMATRIX_SP(ATOM_,NM,NF,MODES_NUM,FIELDS,VALUES_,ROW_INDEX_,COLUMN_,INFO)
!!$  
!!$  USE TYPES !(modes.f90)
!!$  USE MERGINGARRAYS !(utils.f90)
!!$  
!!$  IMPLICIT NONE
!!$  INTEGER                  ,INTENT(IN)    :: NM,NF
!!$  TYPE(MODE), DIMENSION(NF),INTENT(INOUT) :: FIELDs
!!$  TYPE(ATOM),               INTENT(IN)    :: ATOM_
!!$  INTEGER, DIMENSION(NM),   INTENT(IN)    :: MODES_NUM
!!$  INTEGER,                  INTENT(INOUT) :: INFO
!!$  COMPLEX*16, DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: VALUES_
!!$  INTEGER,    DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: COLUMN_
!!$  INTEGER,    DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: ROW_INDEX_
!!$  
!!$  
!!$  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: V_AUX
!!$  COMPLEX*16, DIMENSION(:),   ALLOCATABLE :: VALUES,ARRAY_AUX,VALUES_OLD,VALUES_OLD_
!!$  INTEGER,    DIMENSION(:),   ALLOCATABLE ::ROW,COLUMN,ARRAYI_AUX,ROW_OLD,COLUMN_OLD,ROW_INDEX,INDEX_ORDERROW
!!$  INTEGER,    DIMENSION(NF)               :: N_FLOQUET
!!$  
!!$  INTEGER :: D_MULTIFLOQUET
!!$  INTEGER m,r,c,index,D,counter, N_MODES,values_dim,D_bare,t
!!$  CHARACTER*1 UPLO
!!$  
!!$  !PARAMETERS REQUIRED TO TEST THE DIAGONALIZATION ROUTINE
!  COMPLEX*16,      DIMENSION(:,:),ALLOCATABLE :: U_F
!  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: E_FLOQUET
!  DOUBLE PRECISION :: E_L,E_R
!!$
!!$  N_MODES = NF
!!$  D_bare = ATOM_%D_BARE
!!$  info   = 0
!!$
!!$  ALLOCATE(VALUES(D_bare*D_bare))
!!$  ALLOCATE(ROW(D_bare*D_bare))
!!$  ALLOCATE(COLUMN(D_bare*D_bare))
!!$
!!$  !DEFINITION OF THE MATRICES: BY HAND  
!!$  DO r=1,NF
!!$     N_FLOQUET(r) = FIELDS(r)%N_FLOQUET
!!$  END DO
!!$    
!!$  ! COORDINATE PACKING OF EACH FIELD
!!$  
!!$  D_MULTIFLOQUET = SIZE(FIELDS(1)%V,1)
!!$  DO m=1,NF
!!$     CALL COORDINATEPACKING(D_bare,FIELDS(m)%V,VALUES,ROW,COLUMN,index,INFO)
!!$     ALLOCATE(FIELDS(m)%VALUES(index))
!!$     IF(m.eq.1) FIELDS(m)%VALUES = VALUES(1:index)
!!$     IF(m.ne.1) FIELDS(m)%VALUES = VALUES(1:index)/2.0
!!$!     write(*,*) real(fields(m)%values)
!!$!     write(*,*) fields(m)%omega
!!$     ALLOCATE(FIELDS(m)%ROW(index))
!!$     FIELDS(m)%ROW = ROW(1:index)
!!$!     WRITE(*,*) FIELDS(m)%ROW
!!$     ALLOCATE(FIELDS(m)%COLUMN(index))
!!$     FIELDS(m)%COLUMN = COLUMN(1:index)
!!$!     WRITE(*,*) FIELDS(m)%COLUMN
!!$     D_MULTIFLOQUET = D_MULTIFLOQUET*(2*N_FLOQUET(m)+1)
!!$!     write(*,*) d_multifloquet
!!$  END DO
!!$!  WRITE(*,*)
!!$
!!$  ! WE HAVE TO REDEFINE VALUES,ROW,COLUMN FOR m > 2 as follows:
!!$  D = D_bare
!!$  index = 1
!!$  DO m=3,NF
!!$     DEALLOCATE(VALUES)
!!$     DEALLOCATE(ROW)
!!$     DEALLOCATE(COLUMN)
!!$     ALLOCATE(VALUES(SIZE(FIELDS(m)%VALUES,1)))
!!$     ALLOCATE(ROW(SIZE(FIELDS(m)%VALUES,1)))
!!$     ALLOCATE(COLUMN(SIZE(FIELDS(m)%VALUES,1)))
!!$     ALLOCATE(ARRAYI_AUX(SIZE(FIELDS(m)%VALUES,1)))
!!$     VALUES = FIELDS(m)%VALUES
!!$     ROW    = FIELDS(m)%ROW
!!$     COLUMN = FIELDS(m)%COLUMN
!!$     DO r=2,2*N_FLOQUET(m-1)+1
!!$        CALL APPENDARRAYS(VALUES,FIELDS(m)%VALUES,INFO)
!!$        ARRAYI_AUX = FIELDS(m)%ROW + (r-1)*D
!!$        CALL APPENDARRAYSI(ROW,ARRAYI_AUX,INFO)
!!$        ARRAYI_AUX =  FIELDS(m)%COLUMN + (r-1)*D
!!$        CALL APPENDARRAYSI(COLUMN,ARRAYI_AUX,INFO)
!!$     END DO
!!$     DEALLOCATE(FIELDS(m)%VALUES)
!!$     DEALLOCATE(FIELDS(m)%ROW)
!!$     DEALLOCATE(FIELDS(m)%COLUMN)
!!$     DEALLOCATE(ARRAYI_AUX)
!!$     ALLOCATE(FIELDS(m)%VALUES(SIZE(VALUES,1)))
!!$     ALLOCATE(FIELDS(m)%ROW(SIZE(VALUES,1)))
!!$     ALLOCATE(FIELDS(m)%COLUMN(SIZE(VALUES,1)))
!!$     FIELDS(m)%VALUES = VALUES
!!$     FIELDS(m)%ROW    = ROW
!!$     FIELDS(m)%COLUMN = COLUMN
!!$!     write(*,*) D
!!$!     write(*,*) abs(values)
!!$!     write(*,*) row
!!$!     write(*,*) column
!!$     D = D*(2*N_FLOQUET(m-1)+1)     
!!$  END DO
!!$  DEALLOCATE(VALUES)
!!$  DEALLOCATE(ROW)
!!$  DEALLOCATE(COLUMN)
!!$  
!!$  
!!$  ! BUILDING THE COORDINATE PACKING OF THE MULTIMODE HAMILTONIAN MATRIX
!!$  D = D_bare
!!$  values_dim = 0  
!!$  ALLOCATE(VALUES_OLD(SIZE(FIELDS(1)%VALUES,1)))
!!$  ALLOCATE(ROW_OLD(SIZE(FIELDS(1)%VALUES,1)))
!!$  ALLOCATE(COLUMN_OLD(SIZE(FIELDS(1)%VALUES,1)))
!!$  ALLOCATE(VALUES(SIZE(FIELDS(1)%VALUES)))
!!$  ALLOCATE(ROW(SIZE(FIELDS(1)%VALUES)))
!!$  ALLOCATE(COLUMN(SIZE(FIELDS(1)%VALUES)))
!!$  VALUES = FIELDS(1)%VALUES
!!$  ROW    = FIELDS(1)%ROW
!!$  COLUMN = FIELDS(1)%COLUMN
!!$  VALUES_OLD = FIELDS(1)%VALUES
!!$  ROW_OLD    = FIELDS(1)%ROW
!!$  COLUMN_OLD = FIELDS(1)%COLUMN  
!!$!  write(*,*)values_old
!!$  DO m=2,NF
!!$     ALLOCATE(ARRAYI_AUX(SIZE(FIELDS(m)%ROW,1)))
!!$     ALLOCATE(VALUES_OLD_(SIZE(VALUES_OLD,1)))
!!$     VALUES_OLD_ = VALUES_OLD
!!$     DO r=1,2*N_FLOQUET(m) + 1
!!$        DO c = 1,2*N_FLOQUET(m) + 1
!!$           index = c-r
!!$           
!!$           IF(index.EQ.0) THEN
!!$              
!!$              DO t=1,SIZE(VALUES_OLD,1)
!!$                 IF(COLUMN_OLD(t).EQ.ROW_OLD(t)) VALUES_OLD_(t) = VALUES_OLD(t) + 1.0*(-N_FLOQUET(m)+c-1)*FIELDS(m)%OMEGA
!!$              END DO
!!$              IF(c.NE.1) THEN
!!$                 CALL APPENDARRAYS(VALUES,VALUES_OLD_,INFO)
!!$                 CALL APPENDARRAYSI(COLUMN,COLUMN_OLD+(c-1)*D,INFO)
!!$                 CALL APPENDARRAYSI(ROW,ROW_OLD+(c-1)*D,INFO)
!!$              ELSE
!!$                 VALUES = VALUES_OLD_
!!$              END IF
!!$           ELSE IF(ABS(index).EQ.1) THEN            
!!$              CALL APPENDARRAYS(VALUES,FIELDS(m)%VALUES,INFO)
!!$              ARRAYI_AUX =  FIELDS(m)%COLUMN + (c-1)*D
!!$              CALL APPENDARRAYSI(COLUMN,ARRAYI_AUX,INFO)
!!$              
!!$              ARRAYI_AUX =  FIELDS(m)%ROW + (r-1)*D
!!$              CALL APPENDARRAYSI(ROW,ARRAYI_AUX,INFO)
!!$           END IF
!!$        END DO
!!$     END DO
!!$     
!!$     D = D*(2*N_FLOQUET(m)+1)
!!$     DEALLOCATE(VALUES_OLD_)
!!$     DEALLOCATE(VALUES_OLD)
!!$     DEALLOCATE(ROW_OLD)
!!$     DEALLOCATE(COLUMN_OLD)
!!$     DEALLOCATE(ARRAYI_AUX)
!!$     
!!$     ALLOCATE(VALUES_OLD(SIZE(VALUES,1)))
!!$     ALLOCATE(ROW_OLD(SIZE(VALUES,1)))
!!$     ALLOCATE(COLUMN_OLD(SIZE(VALUES,1)))
!!$     VALUES_OLD = VALUES
!!$     ROW_OLD    = ROW
!!$     COLUMN_OLD = COLUMN
!!$  END DO
!!$  
!!$  DEALLOCATE(VALUES_OLD)
!!$  DEALLOCATE(ROW_OLD)
!!$  DEALLOCATE(COLUMN_OLD)
!!$  
!!$  ! BUILDING THE VARIATION CRS PACKING OF THE MULTIMODE HAMILTONIAN MATRIX, USING THE COORDINATE PACKING
!!$
!!$  !  DO r=1,size(VALUES,1)
!!$  !     WRITE(*,*) ROW(r),COLUMN(r),REAL(VALUES(r))
!!$  !  END DO
!!$  !    WRITE(*,*) SIZE(VALUES,1)
!!$  !    WRITE(*,*) REAL(VALUES)
!!$  !    WRITE(*,*) ROW
!!$  !    WRITE(*,*) COLUMN
!!$  !WRITE(*,*) 
!!$  ALLOCATE(INDEX_ORDERROW(SIZE(ROW,1)))
!!$  CALL QUICK_SORT_INTEGERS(ROW,INDEX_ORDERROW,SIZE(ROW,1))
!!$  ROW    = ROW(INDEX_ORDERROW)
!!$  COLUMN = COLUMN(INDEX_ORDERROW)
!!$  VALUES = VALUES(INDEX_ORDERROW)
!!$  
!!$  ! write(*,*) row
!!$  ! WRITE(*,*)
!!$  ! write(*,*) column
!!$  ! WRITE(*,*)
!!$  ! WRITE(*,*) D_MULTIFLOQUET
!!$  
!!$  
!!$  ALLOCATE(ROW_INDEX(D_MULTIFLOQUET+1))
!!$  ROW_INDEX = -1
!!$  ROW_INDEX(D_MULTIFLOQUET+1) = SIZE(VALUES,1)+1
!!$  
!!$  counter = 1
!!$  
!!$  D = 1
!!$  ROW_INDEX(1)=1
!!$  DO r = 2,SIZE(ROW,1)
!!$     IF(ROW(r).EQ.ROW(r-1)) THEN
!!$        counter = counter +1        
!!$     ELSE
!!$        D = D + counter
!!$        ROW_INDEX(ROW(r)) =D
!!$        counter = 1
!!$     END IF
!!$  END DO
!!$  !  WRITE(*,*) D_MULTIFLOQUET,nf
!  E_L = -6.0
!  E_R =  6.0
!  ALLOCATE(E_FLOQUET(D_MULTIFLOQUET))
!  ALLOCATE(U_F(D_MULTIFLOQUET,D_MULTIFLOQUET))
!  CALL MKLSPARSE_FULLEIGENVALUES(D_MULTIFLOQUET,SIZE(VALUES,1),VALUES,ROW_INDEX,COLUMN,E_L,E_R,E_FLOQUET,U_F,INFO)
!  write(*,*) E_FLOQUET
!!$  ALLOCATE(VALUES_(SIZE(VALUES,1)))
!!$  ALLOCATE(ROW_INDEX_(SIZE(ROW_INDEX,1)))
!!$  ALLOCATE(COLUMN_(SIZE(COLUMN,1)))
!!$  
!!$  VALUES_    = VALUES
!!$  ROW_INDEX_ = ROW_INDEX
!!$  COLUMN_    = COLUMN
!!$  
!!$  DEALLOCATE(VALUES)
!!$  DEALLOCATE(ROW)
!!$  DEALLOCATE(COLUMN)
!!$  DEALLOCATE(ROW_INDEX)
!!$
!!$END SUBROUTINE MULTIMODEFLOQUETMATRIX_SP ! _SP  sparse packing

!!$SUBROUTINE MKLSPARSE_FULLEIGENVALUES(D,DV,VALUES,ROW_INDEX,COLUMN,E_L,E_R,E_FLOQUET,U_F,INFO)
!!$
!!$!CALCULATES THE ENERGY SPECTRUM OF THE MATRIX REPRESENTED BY VALUES, ROW_INDEX AND COLUMN
!!$! D (IN), MATRIX DIMENSION == NUMBER OF EIGENVALUES
!!$! DV (IN), NUMBER OF VALUES != 0
!!$! VALUES (IN) ARRAY OF VALUES
!!$! ROW_INDEX (IN), ARRAY OF INDICES
!!$! COLUMN (IN),    ARRAY OF COLUMN NUMBERS
!!$! E_L (IN),       LEFT BOUNDARY OF THE SEARCH INTERVAL
!!$! E_R (IN),       RIGHT BOUNDARY OF THE SEARCH INTERVAL
!!$! E_FLOQUET (OUT), ARRAY OF EIGENVALUES
!!$! INFO     (INOUT)  ERROR FLAG
!!$
!!$  USE FEAST
!!$  IMPLICIT NONE
!!$  INTEGER,                          INTENT(IN)    :: D,DV
!!$  COMPLEX*16,       DIMENSION(DV),  INTENT(INOUT) :: VALUES
!!$  INTEGER,          DIMENSION(DV),  INTENT(INOUT) :: COLUMN
!!$  INTEGER,          DIMENSION(D+1), INTENT(INOUT) :: ROW_INDEX
!!$  DOUBLE PRECISION,                 INTENT(IN)    :: E_L,E_R
!!$  DOUBLE PRECISION, DIMENSION(D),   INTENT(OUT)   :: E_FLOQUET
!!$  COMPLEX*16,       DIMENSION(D,D), INTENT(OUT)   :: U_F
!!$  INTEGER,                          INTENT(INOUT) :: INFO
!!$
!!$  CHARACTER*1 UPLO
!!$!  DOUBLE PRECISION :: Emin,Emax
!!$  
!!$
!!$  ! ----- 2. RUN FEASTINIT
!!$  
!!$  CALL feastinit(fpm)
!!$  
!!$  ! ----- 3. SOLVE THE STANDARD EIGENVALUE PROBLEM
!!$
!!$  M0  = D  ! number of eignevalues requested
!!$  
!!$  ALLOCATE(E(M0)) ! array of eigenvalues
!!$  ALLOCATE(RES(M0)) ! array of residuals
!!$  ALLOCATE(X(D,M0)) ! matrix with eigenvectors
!!$  E   = 0
!!$  RES = 0
!!$  X   = 0
!!$
!!$
!!$  info_FEAST = 0
!!$  Emin = E_L!-15.0 ! SEARCH INTERVAL: LOWER BOUND
!!$  Emax = E_R! 15.0 ! SEARCH INTERVAL: UPPER BOUND
!!$  UPLO = 'F'
!!$!  WRITE(*,*) d,SIZE(VALUES,1),SIZE(ROW_INDEX,1),SIZE(COLUMN,1),SIZE(FPM,1)
!!$  CALL zfeast_hcsrev(UPLO,D,VALUES,ROW_INDEX,COLUMN,fpm,epsout,loop, &
!!$       &   Emin,Emax,M0,E,X,M1,res,info_FEAST)
!!$  !PRINT  *,'FEAST OUTPUT INFO ',info
!!$  !WRITE(*,*) 'GUESSED NUMBER OF EIGENVALUES:', m0
!!$  !WRITE(*,*) 'NUMBER OF EIGENVALUES FOUND:', m1
!!$  !WRITE(*,*) 'EIGENVALUES:', E, D
!!$  !WRITE(*,*) 'EIGENVECTORS:', ABS(X)
!!$  !CALL WRITE_MATRIX(D_MULTIFLOQUET,D_MULTIFLOQUET,ABS(X))
!!$  E_FLOQUET = E
!!$  U_F       = X
!!$
!!$  DEALLOCATE(RES)
!!$  DEALLOCATE(X)
!!$  DEALLOCATE(E)
!!$  
!!$END SUBROUTINE MKLSPARSE_FULLEIGENVALUES
!!$
