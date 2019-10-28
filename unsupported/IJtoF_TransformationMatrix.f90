!!$! NOTE: this program only works for L=0!
!!$MODULE delta_fun
!!$  IMPLICIT NONE
!!$  
!!$CONTAINS
!!$
!!$  FUNCTION delta_kr(M,N)
!!$    IMPLICIT NONE
!!$    DOUBLE PRECISION :: M,N
!!$    INTEGER DELTA_KR    
!!$    IF(M.EQ.N) DELTA_KR = 1
!!$    IF(M.NE.N) DELTA_KR = 0    
!!$  END FUNCTION delta_kr 
!!$END MODULE delta_fun
!!$
!!$PROGRAM IJtoF_MatrixTest
!!$
!!$  IMPLICIT NONE
!!$  DOUBLE PRECISION L,S,I,Total_states_LSI,F,J
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CLEBSH_GORDAN_JtoF
!!$  INTEGER N,i_
!!$
!!$  L = 0.0
!!$  S = 0.5
!!$  I = 1.5
!!$  
!!$  !----- Counting the number of states  
!!$  J = L+S
!!$  Total_states_LSI = 0    
!!$  DO WHILE(J.GE.ABS(L-S))
!!$     F = I+J
!!$     DO WHILE(F.GE.ABS(J-I))          
!!$        Total_states_LSI = Total_states_LSI + 2*F + 1
!!$        F = F - 1
!!$     END DO
!!$     J= J - 1 ! only two values allowed for J: J = L+1/2 and J=L-1/2     
!!$  END DO
!!$  J = L+S  ! reseting value of J modified before
!!$  
!!$  N = INT(Total_states_LSI)
!!$  write(*,*) Total_states_LSI,J,I
!!$  ALLOCATE(CLEBSH_GORDAN_JtoF(N,N))
!!$  
!!$  CALL IJtoF_Matrix(L,S,I,N,CLEBSH_GORDAN_JtoF)
!!$
!!$  write(*,208) (CLEBSH_GORDAN_JtoF(:,i_), i_=1,N)
!!$  
!!$208 FORMAT(8E15.6)
!!$END PROGRAM IJtoF_MatrixTest
 

SUBROUTINE IJtoF_Matrix(L,S,I,Total_states_LSI,CLEBSH_GORDAN_JtoF)
  
  !USE delta_fun
  USE FUNCIONES
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) ::L,S,I
  INTEGER, INTENT(IN) :: Total_states_LSI
  DOUBLE PRECISION, DIMENSION(Total_states_LSI,Total_states_LSI), INTENT(OUT) :: CLEBSH_GORDAN_JtoF


  DOUBLE PRECISION, DIMENSION(Total_states_LSI,4) :: STATE
  DOUBLE PRECISION J,F,J_,CB,mF
  DOUBLE PRECISION mI,mI_,mJ,mJ_,M,M_
  INTEGER u,r,p,o,v,k,w

  !Build the Clebsh-Gordan Matrix for {|IJmImJ> --> |IJFmF>}

!  write(*,*) Total_states_LSI
  r = 0
  p = 0    
  J = L+S 
  CLEBSH_GORDAN_JtoF  = 0.0
  DO u =0,int(2*(L+S+I))                                           ! POSSIBLE VALUES FOR mF
     mF = L+S+I - u
     DO v=0,0
        J_ = L+S-v                                                 ! SET THE ELECTRONIC ANGULAR MOMENTUM J_
        F = L+S+I - v                                              ! SET THE TOTAL ANGULAR MOMENTUM F
        o = 0
        DO WHILE (ABS(mF).LE.F .AND. F.GE.ABS(J_-I))
           r = r + 1
           p = r - o
           
           STATE(r,1) = I
           STATE(r,2) = J_
           STATE(r,3) = F
           STATE(r,4) = mF
           
           
           mJ = J
           mI = mF-J
           if(abs(mI).GT.I)  then
              mJ = -J
              mI = mF + J
           end if
!           write(*,*) mJ,mI
           DO WHILE ((mJ.GE.-J).AND.(ABS(mI).LE.I))
              !r = r + 1
              
              
              !DO w = 0,int(2*I)
              !mI = I-w                                             ! POSSIBLE VALUES OF mI  
              ! INITIAL VALUE OF J
              !DO WHILE(J.GE.ABS(L-S))
              !DO k = 0,int(2*J)
              !   mJ= J - k                                      ! POSSIBLE VALUES OF mJ
              
              CALL cb_fun(int(2.0*J),int(2.0*I),int(2.0*F), &
                   & int(2.0*mJ),int(2.0*mI),int(2.0*mF),CB) ! CLEBSH-GORDAN COEFFICIENTE <F,mF|J,I,mJ,mi>               
                 
              CLEBSH_GORDAN_JtoF(r,p) = CB*delta_kr(J_,J)
!              write(*,*) r,p,mI,mJ,F,mF,CLEBSH_GORDAN_JtoF(r,p)
              p = p+1
              mJ = mJ - 1
              mI = mI + 1
              !   J= J - 1                                         ! only two values allowed for J: J = L+1/2 and J=L-1/2     
              !END DO
              o=1
           END DO
           F =F-1                                                 ! change the value of F
        END DO
     END DO
  END DO

END SUBROUTINE IJtoF_Matrix
FUNCTION delta_kr(M,N)
  IMPLICIT NONE
  DOUBLE PRECISION :: M,N
  INTEGER DELTA_KR    
  IF(M.EQ.N) DELTA_KR = 1
  IF(M.NE.N) DELTA_KR = 0    
END FUNCTION delta_kr
