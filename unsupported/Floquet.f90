SUBROUTINE FLOQUETINIT


! calculate the dimenson of the Hilbert space
! initialize all the matrices required for a full Floquet calcuations
! Calculate the nuclear, electron and total angular momentum operators

  USE physical_constants ! Standard Module with constants
  USE ATOMIC_PROPERTIES  ! gF, F , etc. factors for several species
  USE subinterface       ! To ubroutines for representation of I and J operators
  USE ARRAYS
  USE FLOQUET            ! Number of floquet modes
  USE SUBINTERFACE_LAPACK
  IMPLICIT NONE

  INTEGER  r,D_F2,INFO,P,r_,p_
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: Energy


  J = L+S
  F    = 2
  gF_2 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
!  write(*,*) gF_2
  F    = 1
  gF_1 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
  G_F  = (g_J-g_I)/16.0


 
  !----- Counting the number of states
  J = L+S
  Total_states_LSI = 0
  DO WHILE(J.GE.ABS(L-S))
     F = I+J
     DO WHILE(F.GE.ABS(J-I))
!        write(*,*) Total_states_LSI,F,I,J
        Total_states_LSI = Total_states_LSI + 2*int(F) + 1
        F = F - 1
     END DO
     J= J - 1 ! only two values allowed for J: J = L+1/2 and J=L-1/2
  END DO
  J = L+S  ! reseting value of J modifeid before

  if(s.eq.0.5 .and. i.LT.0.5) total_states_lsi = 2
   !------ ALLOCATE NEEDED ARRAYS: Hamiltonian and Lapack
  ALLOCATE(Energy(TOTAL_STATES_LSI))
!  write(*,*) Total_states_LSI
!  ALLOCATE(H_hyperfine(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_IJ(Total_states_LSI,Total_states_LSI))
  ALLOCATE(U_ZEEMAN(Total_states_LSI,Total_states_LSI))
  ALLOCATE(U_RF(Total_states_LSI,Total_states_LSI))
  ALLOCATE(HAMILTONIAN(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_AUX(Total_states_LSI,Total_states_LSI)) 
  ALLOCATE(H_RF(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_RF_DAGGER(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_MW(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_ALPHA(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_ALPHA_DAGGER(Total_states_LSI,Total_states_LSI))
  !ALLOCATE(Identity(Total_states_LSI,Total_states_LSI))
  ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
  ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
  ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
  ALLOCATE(jz_dash(Total_states_LSI,Total_states_LSI))
  ALLOCATE(I_x(Total_states_LSI,Total_states_LSI))
  ALLOCATE(I_y(Total_states_LSI,Total_states_LSI))
  ALLOCATE(I_z(Total_states_LSI,Total_states_LSI))
  ALLOCATE(CLEBSH_GORDAN_JtoF(Total_states_LSI,Total_states_LSI))
  !ALLOCATE(Z_M((2*N_floquet+1)*Total_states_LSI,(2*N_floquet+1)*Total_states_LSI))
  !ALLOCATE(Fx(int(2*Ftotal+1),int(2*Ftotal+1)))
  !ALLOCATE(Fy(int(2*Ftotal+1),int(2*Ftotal+1)))
  !ALLOCATE(Fz(int(2*Ftotal+1),int(2*Ftotal+1)))
  !ALLOCATE(Hamiltonian_F(int(2*Ftotal+1),int(2*Ftotal+1)))
  !ALLOCATE(Identity_F(int(2*Ftotal+1),int(2*Ftotal+1)))
  ALLOCATE(g_F_matrix(Total_states_LSI,Total_states_LSI))
  ALLOCATE(Fx(Total_states_LSI,Total_states_LSI))
  ALLOCATE(Fy(Total_states_LSI,Total_states_LSI))
  ALLOCATE(Fz(Total_states_LSI,Total_states_LSI))
  ALLOCATE(Fz_DASH(Total_states_LSI,Total_states_LSI))
  ALLOCATE(Hamiltonian_F(Total_states_LSI,Total_states_LSI))
  ALLOCATE(Identity_F(Total_states_LSI,Total_states_LSI))
  !ALLOCATE(index_state(SUBSPACES*int(2*Ftotal+1)))
  ALLOCATE(F_t(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_w(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_M(Total_states_LSI,Total_states_LSI))
  ALLOCATE(H_J(Total_states_LSI,Total_states_LSI))
  !ALLOCATE(H_FLOQUET((2*N_floquet(1)+1)*Total_states_LSI,(2*N_floquet(1)+1)*Total_states_LSI))
!  ALLOCATE(H_FLOQUET_COPY((2*N_floquet(1)+1)*Total_states_LSI,(2*N_floquet(1)+1)*Total_states_LSI))

  !H_hyperfine   = 0.0
  !H_AUX         = 0.0
  !HAMILTONIAN   = 0.0
  !H_RF          = 0.0
  !H_FLOQUET     = 0.0
  !H_OLD         = 0.0
  !Identity      = 0.0
  j_x           = 0.0
  j_y           = 0.0
  j_z           = 0.0
  I_x           = 0.0
  I_y           = 0.0
  I_z           = 0.0
  CLEBSH_GORDAN_JtoF = 0.0
  !Z_M           = 0.0
  Hamiltonian_F = 0.0
  Fx            = 0.0
  Fy            = 0.0
  Fz            = 0.0
  !index_state   = 0
  F_t           = 0
!  H_w           = 0

  !Identity = 0.0
  !DO r=1,Total_states_LSI
  !   Identity(r,r) =1.0
  !END DO

  Identity_F = 0.0
  DO r=1,SIZE(Identity_F,1)
     Identity_F(r,r) = 1.0
  END DO

  !---------------- Build the angular momentum operators in the JImImJ basis,
  !---------------- in decresing order of m=mI+mJ, with mJ=1/2,-1/2
  CALL I_and_J_representations(j_x,j_y,j_z,I_x,I_y,I_z,L,S,I)
!  write(*,*) 'jy',L,S,I
  
!  CALL WRITE_MATRIX(j_x)
!  CALL WRITE_MATRIX(j_y)
!  CALL WRITE_MATRIX(j_z)
!
!  CALL WRITE_MATRIX(I_x)
!  CALL WRITE_MATRIX(I_y)
!  CALL WRITE_MATRIX(I_z)
  !  CALL IJtoF_Matrix(L,S,I,INT(Total_states_LSI),CLEBSH_GORDAN_JtoF)
  !CALL F_representation(Fx,Fy,Fz,Ftotal)

  Fx = J_x+I_x
  Fy = J_y+I_y
  Fz = J_z+I_z

!  CALL WRITE_MATRIX(J_Z+I_Z)

  IF(TOTAL_STATES_LSI.EQ.8) THEN
     F_t = 0
     DO r=1,3
        F_t(r,r) = 1
     END DO
     DO r=4,8
        F_t(r,r) = -1
     END DO
  END IF


  IF(TOTAL_STATES_LSI.EQ.8) THEN
     g_F_matrix = 0
     DO r=1,3
        g_F_matrix(r,r) = gF_1
     END DO
     DO r=4,8
        g_F_matrix(r,r) = gF_2
     END DO
  END IF

  HAMILTONIAN_F = 0.00001*(g_J*J_z+g_I*I_z) + 2*(MATMUL(I_Z,J_Z) +MATMUL(I_x,J_x) - MATMUL(I_y,J_y) )

  CALL LAPACK_FULLEIGENVALUES(HAMILTONIAN_F,TOTAL_STATES_LSI,Energy,INFO)
  
  Fx=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(Fx,HAMILTONIAN_F))
  Fy=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(Fy,HAMILTONIAN_F))
  Fz=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(Fz,HAMILTONIAN_F))
  H_IJ = MATMUL(I_X,J_X)+MATMUL(I_Z,J_Z)-MATMUL(I_Y,J_Y)
  H_IJ=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(H_IJ,HAMILTONIAN_F))
 
  DO r=1,TOTAL_STATES_lsi
     DO P=1,TOTAL_STATES_lsi
        IF(ABS(HAMILTONIAN_F(r,p)).LT.1E-2) HAMILTONIAN_F(r,p) = 0.0
        IF(ABS(H_IJ(r,p)).LT.1E-2) H_IJ(r,p) = 0.0
        IF(ABS(Fx(r,p)).LT.1E-2) Fx(r,p) = 0.0
        IF(ABS(Fy(r,p)).LT.1E-2) Fy(r,p) = 0.0
        IF(ABS(Fz(r,p)).LT.1E-2) Fz(r,p) = 0.0
     END DO
  END DO
  CLEBSH_GORDAN_JtoF = real(HAMILTONIAN_F)
  
 ! CALL WRITE_MATRIX(FZ)

  DO r=1,TOTAL_STATES_LSI
     IF(r.le.3) r_ = -(r-2)
     IF(r.gt.3) r_ =   r-6
     DO p=1,TOTAL_STATES_LSI
        IF(p.le.3) p_ = -(p-2)
        IF(p.gt.3) p_ =   p-6
        H_W(r,p) = r_ - p_        
     END DO
  END DO

  Fz_dash = 0.0
  DO r=1,TOTAL_STATES_LSI
     ! DC DRESSED HAMILTONIAN IN A FRAME OF REFERENCE WHERE THE TOP FIELD IS STATIC
     FZ_DASH(r,r) = int(g_F_matrix(r,r)/abs(g_F_matrix(r,r)))*Fz(r,r)
  END DO
!  CALL WRITE_MATRIX(1.0D0*(Fz))
!  CALL WRITE_MATRIX(1.0D0*(Fz_dash))
!  WRITE(*,*) Fz(3,3),Fz_dash(3,3)

  DO r=1,TOTAL_STATES_LSI
     IF(r.le.3) r_ = (r-2)
     IF(r.gt.3) r_ =  r-6
     DO p=1,TOTAL_STATES_LSI
        IF(p.le.3) p_ = (p-2)
        IF(p.gt.3) p_ =  p-6
        H_M(p,r) = p_ - r_        
!        write(*,*) p,r,Fz_dash(p,p) , Fz_dash(r,r),Fz_dash(p,p) - Fz_dash(r,r)
!        H_M(p,r) = Fz_dash(p,p) - Fz_dash(r,r)
     END DO
  END DO

!  CALL WRITE_MATRIX(1.0D0*H_M)

  Jz_dash = 0.0
  
  DO r=1,TOTAL_STATES_LSI  
     Jz_dash(r,r) = INT(g_F_matrix(r,r)/abs(g_F_matrix(r,r)))
  END DO

!  CALL WRITE_MATRIX(1.0D0*(Jz_DASH))
  DO r=1,TOTAL_STATES_LSI
!     IF(r.le.3) r_ = (r-2)
!     IF(r.gt.3) r_ =  r-6
     DO p=1,TOTAL_STATES_LSI
!        IF(p.le.3) p_ = (p-2)
!        IF(p.gt.3) p_ =  p-6
!        H_M(p,r) = p_ - r_        
        H_J(p,r) = Jz_dash(p,p) - Jz_dash(r,r)  
     END DO
  END DO

!  CALL WRITE_MATRIX(1.0D0*H_J)
!  CALL write_MATRIX(MATMUL(TRANSPOSE(CLEBSH_GORDAN_JtoF),CLEBSH_GORDAN_JtoF))
!  CALL write_MATRIX(MATMUL(CLEBSH_GORDAN_JtoF,TRANSPOSE(CLEBSH_GORDAN_JtoF)))
!!$  CALL WRITE_MATRIX(CLEBSH_GORDAN_JtoF)
!!$  CALL WRITE_MATRIX(Fx)
!!$  CALL WRITE_MATRIX(Fy)
!!$  call WRITE_MATRIX(Fz)
!!$  CALL WRITE_MATRIX(real(H_IJ))
  
END SUBROUTINE FLOQUETINIT


SUBROUTINE COUPLINGMATRICES(FIELD, INFO)

!!$ THE FIRST MODE IS STATIC AND PRODUCES ZEEMAN SHIFTS.
!!$ COUPLINGS OF OTHER MODES ARE TRANSFORMED TO THE BASIS OF ZEEMAN STATES VIA U_ZEEMAN
!!$ FIELD, INPUT, FIELD-TYPE, FIEL PARAMETERS
!!$ INFO, INOUT, INTEGER, ERROR FLAG

  USE FLOQUET
  USE ARRAYS
  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE_LAPACK

  IMPLICIT NONE
  TYPE(MODE), DIMENSION(MODES_NUM),INTENT(INOUT) :: FIELD
  INTEGER,                         INTENT(INOUT) :: INFO

  INTEGER m
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: E_ZEEMAN,mFs
  DOUBLE PRECISION :: RESONANTrfFREQUENCY
  ALLOCATE(E_ZEEMAN(TOTAL_STATES_LSI))
  ALLOCATE(mFs(TOTAL_STATES_LSI))


  DO m=1,total_states_lsi
     if(m.lt.4) mFs(m) = 2-m
     if(m.ge.4) mFs(m) = m-6
  end do
  
  U_ZEEMAN = 0.0
  !write(*,*) "3",info
  DO m=1,MODES_NUM

     IF(m.EQ.1) THEN
        FIELD(m)%V = A*(MATMUL(I_x,J_x) - MATMUL(I_y,J_y) + MATMUL(I_z,J_z)) + &
             &       mu_B*(g_J*(FIELD(m)%Bx*J_x  + DCMPLX(0.0,-1.0)*FIELD(m)%By*J_y + FIELD(m)%Bz*J_z) + &
             &            (g_I*(FIELD(m)%Bx*I_x  + DCMPLX(0.0,-1.0)*FIELD(m)%By*I_y + FIELD(m)%Bz*I_z)))
        FIELD(m)%V = FIELD(m)%V/A
        U_ZEEMAN = FIELD(m)%V 
        !call write_matrix(real(u_zeeman))
        CALL LAPACK_FULLEIGENVALUES(U_ZEEMAN,SIZE(FIELD(m)%V,1),E_ZEEMAN,INFO)
        !write(*,*) "3",info
        !call write_matrix(real(u_zeeman))
        FIELD(m)%V = MATMUL(TRANSPOSE(CONJG(U_ZEEMAN)),MATMUL(FIELD(m)%V,U_ZEEMAN))
!        write(*,*) sqrt(abs(FIELD(m)%Bx)**2+abs(FIELD(m)%By)**2+abs(FIELD(m)%Bz)**2),(real(E_zeeman))-0.0*field(2)%OMEGA*hbar*mFs/A
 !       RESONANTrfFREQUENCY = A*ABS(E_ZEEMAN(3)-E_ZEEMAN(2))/HBAR
 !       call write_matrix(real(FIELD(m)%V))
     END IF

     IF(m.GT.1) THEN
        FIELD(m)%V = mu_B*(g_J*(FIELD(m)%Bx*J_x  + DCMPLX(0.0,-1.0)*FIELD(m)%By*J_y + FIELD(m)%Bz*J_z) + &
             &            (g_I*(FIELD(m)%Bx*I_x  + DCMPLX(0.0,-1.0)*FIELD(m)%By*I_y + FIELD(m)%Bz*I_z)))
        FIELD(m)%V = FIELD(m)%V/A
        FIELD(m)%V = MATMUL(TRANSPOSE(CONJG(U_ZEEMAN)),MATMUL(FIELD(m)%V,U_ZEEMAN))
!        IF(m.EQ.2) FIELD(m)%OMEGA = RESONANTrfFREQUENCY
!       call write_matrix(ABS(FIELD(m)%V))
     END IF
  END DO

  KD = SIZE(U_ZEEMAN,1)
  DO m=2,MODES_NUM-1
     KD = KD*(2*N_floquet(m)+1)
  END DO
  KD = KD + SIZE(U_ZEEMAN,1) - 1

  DEALLOCATE(E_ZEEMAN)
  DEALLOCATE(mFs)

END SUBROUTINE COUPLINGMATRICES

SUBROUTINE SINGLEMODEFLOQUETMATRIX(ATOM_,FIELD,INFO)

! SET UP THE HAMILTONIAN MATRIX PERFROMING A TRANSFORMATION TO A DOUBLE
! ROTATING FRAME
! ATOM_, IN, ATOM-TYPE, PARAMETERS OF AN ATOM
! FIELD, IN, FIELD-TYPE, FIELD PARAMETERS
! INFO, INOUT, INTEGER, ERROR FLAG


  USE FLOQUET
  USE ARRAYS
  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE_LAPACK

  IMPLICIT NONE
  INTEGER,                        INTENT(INOUT) :: INFO
  TYPE(MODE),DIMENSION(MODES_NUM),INTENT(IN)    :: FIELD
  TYPE(ATOM),                     INTENT(IN)    :: ATOM_

  INTEGER m,n,D,r,o
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_TEMP,H_STATIC,COUPLING,Z_M_COPY
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E_DRESSED
  INTEGER index_l_lower,index_l_upper,index_r_lower,index_r_upper


  D        = ATOM_%D_BARE
  ALLOCATE(H_FLOQUET_COPY(D,D))

  !write(*,*) D
  H_FLOQUET_COPY = FIELD(1)%V  ! STATIC HAMILTONIAN

  !DEALLOCATE(COUPLING)
  DO n=2,MODES_NUM_DRESSING!2,2!MODES_NUM  ! RUN OVER EACH MODE

     ! D : UPDATED AT THE ENDO OF THE LOOP. DIMENSION OF THE MULTIMODE FLOQUET MATRIX

     H_STATIC  = H_FLOQUET_COPY
     DEALLOCATE(H_FLOQUET_COPY)

     ALLOCATE(IDENTITY(D,D))
     IDENTITY  = 0.0
     DO m= 1,D
        IDENTITY(m,m) = 1.0
     END DO

     !     WRITE(*,*) n,"COUPLIGN ALLOCATE",D
     ALLOCATE(COUPLING(D,D))
     COUPLING  = 0.0

     DO r=1,(2*N_FLOQUET_DRESSING(n-1)+1)

        index_l_lower = ATOM_%D_BARE*(r - 1) + 1
        index_l_upper = ATOM_%D_BARE*(r - 1) + ATOM_%D_BARE
        index_r_lower = index_l_lower
        index_r_upper = index_l_upper
        COUPLING(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
             &     FIELD(n)%V  ! COUPLING MATRIX OF MODE n
        !        write(*,*) n,r,index_l_lower,index_l_upper
     END DO

     D = D*(2*N_FLOQUET_DRESSING(n)+1)
     ALLOCATE(H_FLOQUET(D,D))
     H_FLOQUET = 0.0


     !    WRITE(*,*) SIZE(COUPLING,1),D
     DO m=-N_FLOQUET_DRESSING(n),N_FLOQUET_DRESSING(n)

        index_l_lower = (m + N_FLOQUET_DRESSING(n)    )*SIZE(COUPLING,1) + 1
        index_l_upper = index_l_lower + SIZE(COUPLING,1) - 1
        index_r_lower =  index_l_lower
        index_r_upper =  index_l_upper
        H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
             &  1.0*H_STATIC + 1.0*m*hbar*FIELD(n)%OMEGA*IDENTITY/A
        !        Observable_extended(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
        !             &  Observable

        IF(m.LT.N_FLOQUET_DRESSING(n)) THEN

           index_l_lower =  (m + N_Floquet_DRESSING(n) + 1)*SIZE(COUPLING,1) + 1
           index_l_upper =  index_l_lower + SIZE(COUPLING,1) - 1
           index_r_lower =  (m + N_Floquet_DRESSING(n)    )*SIZE(COUPLING,1) + 1
           index_r_upper =  index_r_lower + SIZE(COUPLING,1) - 1
           H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
                &     0.5*COUPLING
           !           write(*,*) index_l_lower,index_l_upper,index_r_lower,index_r_upper

           index_l_lower =  (m + N_Floquet_DRESSING(n)    )*SIZE(COUPLING,1) + 1
           index_l_upper =  index_r_lower + SIZE(COUPLING,1)  - 1
           index_r_lower =  (m + N_Floquet_DRESSING(n) + 1)*SIZE(COUPLING,1) + 1
           index_r_upper =  index_l_lower + SIZE(COUPLING,1)  - 1
           H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
                &     0.5*TRANSPOSE(CONJG(COUPLING))
        END IF

     END DO

     DEALLOCATE(IDENTITY)
     DEALLOCATE(COUPLING)

  END DO

END SUBROUTINE SINGLEMODEFLOQUETMATRIX

