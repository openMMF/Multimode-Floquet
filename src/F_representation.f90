SUBROUTINE F_representation(Fx,Fy,Fz,Ftotal)

  USE FUNCIONES
  
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT):: Fx,Fy,Fz
  DOUBLE PRECISION, INTENT(IN) :: Ftotal
  !INTEGER, INTENT(IN) :: Ftotal_

  !DOUBLE PRECISION
  INTEGER k,p,N_k
  double precision k_!,Ftotal

  Fx = 0.0
  Fy = 0.0 
  Fz = 0.0
!  Ftotal = 1.0*Ftotal_
!  write(*,*) "with F",ftotal
!!$  N_k = 2*Ftotal + 1
!!$  DO k=-int(Ftotal),int(Ftotal),1
!!$     
!!$     Fz(k+int(Ftotal) + 1,k+int(Ftotal) + 1) = 1.0*k
!!$     
!!$  END DO
!!$
!!$   DO k=-int(Ftotal),int(Ftotal),1
!!$     DO p=-int(Ftotal),int(Ftotal),1
!!$        
!!$        Fx(k+int(Ftotal+1),p+int(Ftotal+1)) = 0.5*(SQRT(Ftotal*(Ftotal+1) - k*(k+1))*delta_kr_int(k+1,p)  +& 
!!$             &  SQRT(Ftotal*(Ftotal+1) - k*(k-1))*delta_kr_int(k-1,p))
!!$        Fy(k+int(Ftotal+1),p+int(Ftotal+1)) = 0.5*(SQRT(Ftotal*(Ftotal+1) - k*(k+1))*delta_kr_int(k+1,p)  -&
!!$             & SQRT(Ftotal*(Ftotal+1) - k*(k-1))*delta_kr_int(k-1,p))
!!$
!!$     END DO
!!$  END DO
  N_k = 2*Ftotal + 1
  k_ = -(N_k - 1.0)/2.0
 
  DO k=1,N_k!-int(Ftotal),int(Ftotal),1
     Fz(k,k) = k_
     k_ = k_ + 1
  END DO

  k_ = -(N_k - 1.0)/2.0  
  DO k=1,N_k!-int(Ftotal),int(Ftotal),1
     DO p=1,N_k!-int(Ftotal),int(Ftotal),1
!        write(*,*) k_,k,p,N_k,Ftotal
        Fx(k,p) = 0.5*(SQRT(Ftotal*(Ftotal+1) - k_*(k_+1))*delta_kr_int(k+1,p)  +& 
             &  SQRT(Ftotal*(Ftotal+1) - k_*(k_-1))*delta_kr_int(k-1,p))
        Fy(k,p) = 0.5*(SQRT(Ftotal*(Ftotal+1) - k_*(k_+1))*delta_kr_int(k+1,p)  -&
             & SQRT(Ftotal*(Ftotal+1) - k_*(k_-1))*delta_kr_int(k-1,p))
        
     END DO
     k_ = k_ + 1
  END DO
!!$  do p=1,int(2*FTotal+1)
!!$     write(*,*) (abs(Fy(p,k)), k=1,int(2*FTotal+1))
!!$  end do
!!$  write(*,*) 
!!$  do p=1,int(2*FTotal+1)
!!$     write(*,*) (abs(Fx(p,k)), k=1,int(2*FTotal+1))
!!$  end do
!!$  write(*,*)

 205 FORMAT(5e15.6)
END SUBROUTINE F_representation
