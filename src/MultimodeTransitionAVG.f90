SUBROUTINE MULTIMODETRANSITIONAVG(D,NM,FIELD,MODES_NUM,U_F_MODES,E_MULTIFLOQUET,D_BARE,U,INFO) 
!!$   AVERAGE TIME EVOLUTION OPERATOR OF A MULTIMODE DRESSED SYSTEM. THE AVERAGE EVOLUTION OPERATOR IS WRITEN IN THE BASIS USED TO EXPRESS THE 
!!$   MULTIMODE FLOQUET HAMILTONIAN
!!$   U : MATRIX OF AVERAGE TRANSITION PROBABILITIES
!!$
!!$  D              (IN)   : DIMENSION OF THE EXTENDED HILBERT SPACE (SIZE OF THE MULTIMODE FLOQUET MATRIX)
!!$  NM             (IN)   : NUMBER OF MODES            
!!$  MODES_NUM      (IN)   : VECTOR (NM) INDICATING THE NUMBER OF HARMONICS OF EACH MODE
!!$  U_F_MODES      (IN)   : TRANSFORMATION, DIMENSOON (D,D) 
!!$  E_MULTIFLOQUET (IN)   : MULTIMODE FLOQUET SPECTRUM
!!$  D_BARE         (IN)   : DIMENSION OF THE BARE HILBERT SPACE
!!$  U              (OUT)  :  MATRIX OF AVERAGE TRANSITION PROBABILITIES
!!$  INFO           (INOUT): (POSSIBLE) ERROR FLAG

  USE TYPES

  IMPLICIT NONE
  TYPE(MODE),DIMENSION(NM), INTENT(IN)     :: FIELD
  INTEGER,   DIMENSION(NM), INTENT(IN)     :: MODES_NUM

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

  DO n=2,NM
     FIELD_INDEX = 2+SUM(MODES_NUM(2:n-1))
     N_FLOQUET(n)=FIELD(FIELD_INDEX)%N_Floquet
!     write(*,*) n,N_FLOQUET(n),field_index,modes_num(n)
     IF(modes_num(n).GT.N_FLOQUET(n)+1) THEN
        WRITE(*,*) "TO BUILD THE EXTENDED HAMILTONIAN THE NUMBER OF FLOQUET MODES MUST BE DEFINED"
        WRITE(*,*) "LARGER THAN THE NUMBER OF FIELD MODES"
        INFO = -10
     END IF
  END DO

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
!!$!  USE FLOQUET
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
!!$     ALLOCATE(U_MODES_n(n)%n(MODES_NUM))
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

