SUBROUTINE F_representation(Fx,Fy,Fz,Ftotal)

  USE FUNCIONES
  
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT):: Fx,Fy,Fz
  DOUBLE PRECISION, INTENT(IN) :: Ftotal
  INTEGER k,p,N_k
  double precision k_

  Fx = 0.0
  Fy = 0.0 
  Fz = 0.0

  N_k = 2*Ftotal + 1
  k_ = -(N_k - 1.0)/2.0
 
  DO k=1,N_k
     Fz(k,k) = k_
     k_ = k_ + 1
  END DO

  k_ = -(N_k - 1.0)/2.0  
  DO k=1,N_k
     DO p=1,N_k
        Fx(k,p) = 0.5*(SQRT(Ftotal*(Ftotal+1) - k_*(k_+1))*delta_kr_int(k+1,p)  +& 
             &  SQRT(Ftotal*(Ftotal+1) - k_*(k_-1))*delta_kr_int(k-1,p))
        Fy(k,p) = 0.5*(SQRT(Ftotal*(Ftotal+1) - k_*(k_+1))*delta_kr_int(k+1,p)  -&
             & SQRT(Ftotal*(Ftotal+1) - k_*(k_-1))*delta_kr_int(k-1,p))
        
     END DO
     k_ = k_ + 1
  END DO

END SUBROUTINE F_representation
