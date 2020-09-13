SUBROUTINE Tunneling_B(D_BARE,N_SITES,N_BODIES,STATE,H_J,INFO)
   
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

END SUBROUTINE Tunneling_B

SUBROUTINE Onsite_twobody_B(D_BARE,N_SITES,N_BODIES,STATE,H_U,INFO)

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

END SUBROUTINE Onsite_twobody_B

SUBROUTINE Tunneling_F(D_BARE,N_SITES,N_BODIES,STATE,H_J,INFO)
   
  USE CREATIONDESTRUCTION
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: N_SITES
  INTEGER,    DIMENSION(2),              INTENT(IN)    :: D_BARE,N_BODIES
  INTEGER,    DIMENSION(D_BARE(1)+D_BARE(2),N_SITES),      INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE(1)+D_BARE(2),D_BARE(1)+D_BARE(2)), INTENT(OUT)   :: H_J
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I
    
  COMPLEX*16, DIMENSION(D_BARE(1),D_BARE(1)) :: T_UP
  COMPLEX*16, DIMENSION(D_BARE(2),D_BARE(2)) :: T_DOWN
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: T
  
  
  INTEGER :: N,I_,J_,k_,l_
  INTEGER, DIMENSION(3) :: N_
  
  N_(1) = D_BARE(1)
  N_(2) = D_BARE(2)
  N_(3) = D_BARE(1)*D_BARE(2)
  ALLOCATE(T(N_(3),N_(3)))
  
  T_UP   = 0
  T_DOWN = 0
  H_J    = 0
  T      = 0
  DO k_=1,N_SITES ! loop though all sites
    DO l_=1,2 ! loop through spin up and spin down
        N = D_BARE(l_)
        DO J_=1,N ! Nested loop through all spin up/down states
            DO I_=J_,N
                STATE_J = STATE(J_ + (l_-1)*D_BARE(1),:) 
                STATE_I = STATE(I_ + (l_-1)*D_BARE(1),:) 
                NEW_STATE = TUNNELING_F_(k_,STATE_J)
                
                
                
                
                !write(*,*) NEW_STATE,dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I))
                IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
                    IF(l_.EQ.1) THEN
                        T_UP(I_,J_) = 1.0
                        T_UP(J_,I_) = 1.0
                    ELSE
                        T_DOWN(I_,J_) = 1.0
                        T_DOWN(J_,I_) = 1.0
                    END IF  

                END IF           
            END DO
        END DO

        DO J_=1,N
            STATE_J = STATE(J_+(l_-1)*D_BARE(1),:) 
            STATE_I = STATE(J_+(l_-1)*D_BARE(1),:) 
            NEW_STATE = TUNNELING_F_(k_,STATE_J)
            IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
                IF(l_.EQ.1) THEN
                    T_UP(J_,J_) = 1.0
                ELSE
                    T_DOWN(J_,J_) = 1.0
                END IF
            END IF           
        END DO
    END DO  
    CALL TENSORMULT(N_,T_UP,T_DOWN,T,INFO)
    H_J = H_J + T
  END DO

END SUBROUTINE Tunneling_F


SUBROUTINE Onsite_twobody_F(D_BARE,N_SITES,N_BODIES,STATE,H_U,INFO)
   
  USE CREATIONDESTRUCTION
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: N_SITES
  INTEGER,    DIMENSION(2),              INTENT(IN)    :: D_BARE,N_BODIES
  INTEGER,    DIMENSION(D_BARE(1)+D_BARE(2),N_SITES),      INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE(1)+D_BARE(2),D_BARE(1)+D_BARE(2)), INTENT(OUT)   :: H_U
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I    

  
  INTEGER :: N,I_,J_,k_,l_
  
  COMPLEX*16, DIMENSION(D_BARE(1),D_BARE(1)) :: T_UP
  COMPLEX*16, DIMENSION(D_BARE(2),D_BARE(2)) :: T_DOWN
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: T

  INTEGER, DIMENSION(3) :: N_
  
  H_U = 0
  T_UP   = 0
  T_DOWN = 0
  N_(1) = D_BARE(1)
  N_(2) = D_BARE(2)
  N_(3) = D_BARE(1)*D_BARE(2)
  ALLOCATE(T(N_(3),N_(3)))
    
  !write(*,*) N_ 
  T = 0.0
  DO l_=1,N_(2)
      DO k_=1,N_(1)
          !write(*,*) N_,(l_-1)*N_(1)+k_,(l_-1)*N_(1)+k_,DOT_PRODUCT(STATE(k_,:),STATE(N_(1)+l_,:))
        T((l_-1)*N_(1)+k_,(l_-1)*N_(1)+k_) = DOT_PRODUCT(STATE(k_,:),STATE(N_(1)+l_,:))
    END DO
  END DO
  H_U = T
END SUBROUTINE Onsite_twobody_F
