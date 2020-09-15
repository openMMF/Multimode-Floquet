DOUBLE PRECISION FUNCTION  D_H(N_SITES,N_PARTICLES,stats) 

  !Calculates the number of states of N_particles Bosons/Ferminos in N_sites lattice sites
  ! Wikipedia:
  ! D_B  = (N_particles + N_sites - 1)!)/(N_particles! * (L-1)!)))
  ! D_F  = N_sites!/(N_particles! (N_sites-N_particles)!)))

    !N_SITES      (IN): INTEGER, number of lattice sites
    !N_PARTICLES  (IN): INTEGER, umber of partcles
    !stats        (IN)! 'F' or 'B' for Fermions or Bosons, respectively

  IMPLICIT NONE
  INTEGER, INTENT(IN):: N_SITES,N_PARTICLES
  CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: stats

  INTEGER          :: N_,K_,i
  DOUBLE PRECISION :: Num,Den,N_STATES

  IF(PRESENT(stats)) THEN
     SELECT CASE (stats)
     CASE("B")

         IF(N_PARTICLES.GT.1) THEN
            N_ = N_SITES+N_PARTICLES-1;
            K_ = N_PARTICLES;
            write(*,*) N_,k_
            Num = 1.0;
            do while (n_.gt.0) 
                Num = Num*n_;
                n_ = n_-1;
            end do
            Den = 1.0;
            do while (k_.gt.0)  
                Den = Den*k_;
                k_ = k_-1;
            end do
         ELSE
             num = n_sites
             den = 1
         END IF

     CASE("F")
        N_STATES = N_SITES
        N_ = N_STATES
        Num = 1.0;
        do while (n_.gt.0) 
           Num = Num*n_
           n_ = n_-1;
        end do
        k_ = N_STATES - N_PARTICLES
        Den = 1.0
        do while (k_.gt.0)
           Den = Den*k_;
           k_ = k_-1;
        end do
        k_ = N_PARTICLES
        do while (k_.gt.0) 
           Den = Den*k_;
           k_ = k_-1;
        end do

     CASE DEFAULT
        D_H = N_SITES
     END SELECT

     D_H = Num/Den;

  ELSE

     D_H = N_SITES

  END IF

END FUNCTION D_H


MODULE CREATIONDESTRUCTION

  implicit none
  private
  public :: A_DAGGER,A_,TUNNELING_,TUNNELING_F_

  interface A_DAGGER
    procedure A_DAGGER_INT,A_DAGGER_REAL
  end interface A_DAGGER

  interface A_
    procedure A_INT,A_REAL
  end interface A_

contains
    !Bosonic tunelling
    FUNCTION TUNNELING_(k,STATE) result(NEW_STATE)

    
        INTEGER, INTENT(IN) :: k
        INTEGER, DIMENSION(:), INTENT(IN) :: STATE
    
        DOUBLE PRECISION, DIMENSION(SIZE(STATE,1)) :: NEW_STATE
            
        NEW_STATE      = STATE 
        IF(STATE(K).GT.0) THEN
            NEW_STATE(k)   = NEW_STATE(k) - 1            
            NEW_STATE(k+1) = NEW_STATE(k+1) + 1            
        ELSE
            NEW_STATE = 0
        END IF  
    
    END FUNCTION TUNNELING_

    !FERMIONIC  tunelling
    FUNCTION TUNNELING_F_(k,STATE) result(NEW_STATE)

    
        INTEGER, INTENT(IN) :: k
        INTEGER, DIMENSION(:), INTENT(IN) :: STATE
    
        DOUBLE PRECISION, DIMENSION(SIZE(STATE,1)) :: NEW_STATE
    
        NEW_STATE      = STATE 
        IF(K.LT.SIZE(STATE,1))THEN
            IF(STATE(K).EQ.1 .AND. STATE(K+1).EQ.0) THEN
                NEW_STATE(k)   = NEW_STATE(k) - 1            
                NEW_STATE(k+1) = NEW_STATE(k+1) + 1
            ELSE
                NEW_STATE = 0
            END IF  
        ELSE
!            write(*,*) 'mw,k',k,state(k),state(k-1)
            IF(STATE(K).EQ.1 .AND. STATE(K-1).EQ.0) THEN
                NEW_STATE(k)   = NEW_STATE(k) - 1            
                NEW_STATE(k-1) = NEW_STATE(k-1) + 1
            ELSE
                NEW_STATE = 0
            END IF  
        END IF
            

    END FUNCTION TUNNELING_F_
   

    !BOSONIC CREATION AND DESTRUCTION OPERATORE
    FUNCTION A_DAGGER_INT(k,STATE) result(NEW_STATE)

    
        INTEGER, INTENT(IN) :: k
        INTEGER, DIMENSION(:), INTENT(IN) :: STATE
    
        DOUBLE PRECISION, DIMENSION(SIZE(STATE,1)) :: NEW_STATE
    
        !WRITE(*,*) k,NEW_STATE
        NEW_STATE    = STATE 
        NEW_STATE(k) = NEW_STATE(k) + 1
        NEW_STATE    = SQRT(1.0*NEW_STATE(K))*NEW_STATE
               !WRITE(*,*) k,NEW_STATE
    END FUNCTION A_DAGGER_INT

    FUNCTION A_DAGGER_REAL(k,STATE) result(NEW_STATE)

    
        INTEGER, INTENT(IN) :: k
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: STATE
    
        DOUBLE PRECISION, DIMENSION(SIZE(STATE,1)) :: NEW_STATE
    
        NEW_STATE    = STATE 
        !WRITE(*,*) k,NEW_STATE
        NEW_STATE(k) = NEW_STATE(k) + 1
        NEW_STATE    = SQRT(1.0*NEW_STATE(K))*NEW_STATE
        !WRITE(*,*) k,NEW_STATE
    END FUNCTION A_DAGGER_REAL

    FUNCTION A_INT(k,STATE) result(NEW_STATE)
    !A_DAGGER |n> = sqrt(n+1) |n+1>
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: K
        INTEGER, DIMENSION(:), INTENT(IN) :: STATE
    
        DOUBLE PRECISION, DIMENSION(SIZE(STATE,1)) :: NEW_STATE

        NEW_STATE    = STATE 
!        WRITE(*,*) k,NEW_STATE
        IF(STATE(k).GT.0) THEN
            NEW_STATE(k) = NEW_STATE(k) - 1
            NEW_STATE    = SQRT(1.0*NEW_STATE(K))*NEW_STATE
        ELSE
            NEW_STATE =0
        END IF
 !       WRITE(*,*) k,NEW_STATE
    END FUNCTION A_INT
    
    FUNCTION A_REAL(k,STATE) result(NEW_STATE)
    !A_DAGGER |n> = sqrt(n+1) |n+1>
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: K
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: STATE
    
        DOUBLE PRECISION, DIMENSION(SIZE(STATE,1)) :: NEW_STATE

        NEW_STATE    = STATE 
        IF(STATE(k).GT.0) THEN
            NEW_STATE(k) = NEW_STATE(k) - 1
            NEW_STATE    = SQRT(1.0*NEW_STATE(K))*NEW_STATE
        ELSE
            NEW_STATE =0
        END IF
    END FUNCTION A_REAL


end module   CREATIONDESTRUCTION


SUBROUTINE Manybody_basis(D_BARE,N_SITES,N_BODIES,STATS,states_occ,INFO)

!   CREATES A MATRIX WITH AS MANY ROWS AS STATES AND AS MANY COLUMNS AS LATTICE SITES
!   EACH ROW OF THIS MATRIX IS A BASIS STATE IN THE OCCUPATION NUMBER
!    
!    D_BARE   (IN) INTEGER  : NUMBER OF STATES
!    N_SITES  (IN) INTEGER  : NUMBER OF LATTICE SITES
!    N_BODIES (IN) INTEGER  : NUMBER OF PARTICLES
!    STATS    (IN) CHAR     : 'F' or '' B
!    states_occ (OUT) (D_BARE,N_SITES))      : STATES LABELLED AS OCCUPATION OF THE LATTICE
!    INFO       (INOUT)     : ERROR FLAG
    
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: D_BARE,N_SITES,N_BODIES
  CHARACTER(LEN=*),INTENT(IN):: stats
  INTEGER, DIMENSION(D_BARE,N_SITES), intent(out) :: states_occ
  INTEGER, INTENT(INOUT) :: INFO

  LOGICAL MORE
  INTEGER i,j,INTEGER_SPACE

  INTEGER, ALLOCATABLE,DIMENSION(:) :: STATE_

  i = 0
  STATES_OCC = 0
  MORE = .FALSE.

  INFO = 0
  SELECT CASE (stats)
  CASE("B")
     ALLOCATE(STATE_(N_SITES))
     STATE_ = 0
     STATE_(1) = N_BODIES  
     DO i=1,D_BARE        
        CALL COMP_NEXT(N_BODIES,N_SITES,STATE_,MORE) ! DEFINED IN subset.f90
        STATES_OCC(i,:) = STATE_
        !write(*,*) state_
     END DO
     
  CASE("F")

     IF(N_BODIES .LE. 2*N_SITES) THEN
        ALLOCATE(STATE_(N_SITES))
        STATE_=0
        j  =1
        INTEGER_SPACE = 2**N_SITES -1
        DO i=1,INTEGER_SPACE
            
            CALL BVEC_NEXT_GRLEX(N_SITES,STATE_) ! DEFINED IN bvec.f90
            IF(SUM(STATE_) .EQ. N_BODIES) THEN
                STATES_OCC(j,:) = STATE_
                !write(*,*) j
                j=j+1
            ELSE
                
            END IF  
        END DO  
     ELSE
        WRITE(*,*) "ERROR"
        WRITE(*,*) "THE NUMBER OF PARTICLES IS LARGER THAN THE NUMBER OF AVAILABLE STATES"
        INFO  = -1
     END IF

  END  SELECT

END SUBROUTINE Manybody_basis


SUBROUTINE TENSORMULT(N,A,B,C,INFO) 
    ! TENSOR MULTIMPLICATION OF MATRICES A AND B TO PRODUCE C 
    !N, INTEGER(3), IN: DIMENSION OF MATRIX A,B,C  N(1,2,3), RESPECTIVELY
    !A  COMPLEX*18, DIMENSION (N(1),N(2)), IN
    !B  COMPLEX*18, DIMENSION (N(1),N(2)), IN
    !C  COMPLEX*18, DIMENSION (N(3),N(3)), OUT
    !INFO: ERROR FLAG
    IMPLICIT NONE
    INTEGER, DIMENSION(3), INTENT(IN) :: N
    INTEGER, INTENT(INOUT) :: INFO
    COMPLEX*16, DIMENSION(N(1),N(1)), INTENT(IN) :: A
    COMPLEX*16, DIMENSION(N(2),N(2)), INTENT(IN) :: B
    COMPLEX*16, DIMENSION(N(3),N(3)), INTENT(OUT) :: C
    
    INTEGER i,j,r
    
    !WRITE(*,*) A
    !WRITE(*,*) B
    !write(*,*) N
    DO i=1,N(2)
        DO j=1,N(2)
            !if(abs(B(i,j)).gt.0) then
            !    do r=1,size(a,1)
            !     write(*,*) abs(A(r,:))
            !    end do
            !    write(*,*)
            !    write(*,*)
            !    write(*,*) (i-1)*N(1)+1,i*N(1),(j-1)*N(1)+1,j*N(1)
            !end if
            C((i-1)*N(1)+1:i*N(1),(j-1)*N(1)+1:j*N(1)) = A*B(i,j)
        END DO
    END DO
        
END SUBROUTINE TENSORMULT
