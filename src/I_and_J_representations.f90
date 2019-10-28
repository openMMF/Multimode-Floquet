SUBROUTINE I_and_J_representations(j_x,j_y,j_z,I_x,I_y,I_z,L,S,I)

  USE FUNCIONES
  
  IMPLICIT  NONE
  DOUBLE PRECISION, DIMENSION(:,:),INTENT(INOUT) :: j_x,j_y,j_z,I_x,I_y,I_z
  DOUBLE PRECISION, INTENT(IN) :: L,S,I
  
  DOUBLE PRECISION J,F,M,M_,mI,mI_,mJ,mJ_
  INTEGER :: r,p,k,v
  J = L+S
  r = 0
  p = 0    
  j_x = 0
  j_y = 0
  j_z = 0
  I_x = 0
  I_y = 0
  I_z = 0
  
  DO k=1,int(2*(I+S)+1)
     M  = I+S - (k-1)
     mJ = J
     mI = M-J
     if(abs(mI).GT.I)  then
        mJ = -J
        mI = M + J
     end if
     DO WHILE ((mJ.GE.-J).AND.(ABS(mI).LE.I))
        r = r + 1
        p = 0        
        DO v=1,int(2*(I+S)+1)
           M_ = I+S - (v-1) 
           mJ_ = J
           mI_ = M_ - J
           if(abs(mI_).GT.I) THEN
              mJ_ = -J
              mI_ = M_ + J
           END if
           DO WHILE ((mJ_.GE.-J) .AND. (ABS(mI_).LE.I))               
              p = p + 1
              j_x(r,p) = 0.5*(sqrt((J-mJ_)*(J+mJ_+1))*delta_kr(mJ,mJ_+1) + &
                   & sqrt((J+mJ_)*(J-mJ_+1))*delta_kr(mJ,mJ_-1))*delta_kr(mI,mI_)
              j_y(r,p) = 0.5*(sqrt((J-mJ_)*(J+mJ_+1))*delta_kr(mJ,mJ_+1) - &
                   & sqrt((J+mJ_)*(J-mJ_+1))*delta_kr(mJ,mJ_-1))*delta_kr(mI,mI_)
              j_z(r,p) = mJ_*delta_kr(mJ,mJ_)*delta_kr(mI,mI_)
              
              I_x(r,p) = 0.5*(sqrt((I-mI_)*(I+mI_+1))*delta_kr(mI,mI_+1) + &
                   & sqrt((I+mI_)*(I-mI_+1))*delta_kr(mI,mI_-1))*delta_kr(mJ,mJ_)
              I_y(r,p) = 0.5*(sqrt((I-mI_)*(I+mI_+1))*delta_kr(mI,mI_+1) - &
                   & sqrt((I+mI_)*(I-mI_+1))*delta_kr(mI,mI_-1))*delta_kr(mJ,mJ_)                    
              I_z(r,p) = mI_*delta_kr(mJ,mJ_)*delta_kr(mI,mI_)
              !write(*,*) r,p,mI,mJ,mI_,mJ_,j_z(r,p)
              mJ_ = mJ_ - 1
              mI_ = mI_ + 1              
              
           END DO
        END DO
        mJ = mJ - 1
        mI = mI + 1        
     END DO
  END DO
END SUBROUTINE I_and_J_representations
