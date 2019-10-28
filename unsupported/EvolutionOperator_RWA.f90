!/*
! * EvolutionOperator_RWA.f90
! *
! *  Created on: 8 Oct 2017
! *      Author: german
! */

!!$SUBROUTINE EVOLUTIONOPERATOR_RWA(D,U_MODES,ENERGIES,&
!!$     & FIELD,T1,T2,P,INFO)
!!$
!!$  USE FLOQUET
!!$  USE TYPES
!!$  USE SUBINTERFACE_LAPACK
!!$
!!$
!!$  IMPLICIT NONE
!!$  INTEGER,                                INTENT(IN)    :: D
!!$  COMPLEX*16,       DIMENSION(D,D),       INTENT(IN)    :: U_MODES
!!$  DOUBLE PRECISION, DIMENSION(D),         INTENT(IN)    :: ENERGIES
!!$  TYPE(MODE),       DIMENSION(MODES_NUM), INTENT(IN)    :: FIELD  ! FIELDS PROPERTIES: FREQUENCY, AMPLITUDES AND PHASES
!!$  DOUBLE PRECISION,                       INTENT(IN)    :: T1,T2
!!$  COMPLEX*16,       DIMENSION(D,D),       INTENT(OUT)   :: P
!!$  INTEGER,                                INTENT(INOUT) :: INFO
!!$
!!$  COMPLEX*16,       DIMENSION(D,D) :: IDENTITY
!!$
!!$  INTEGER m
!!$
!!$  P = 0.0
!!$  DO m=1,D
!!$     P(m,m) = EXP(-DCMPLX(0.0,1.0)*ENERGIES(m)*(T2-T1))
!!$  END DO
!!$
!!$  P = MATMUL(U_MODES,MATMUL(P,TRANSPOSE(CONJG(U_MODES))))
!!$
!!$END SUBROUTINE EVOLUTIONOPERATOR_RWA


SUBROUTINE EVOLUTIONOPERATOR_RWA_qubit(D,FIELD,T1,T2,P_BARE,INFO)

  USE FLOQUET
  USE TYPES
  USE SUBINTERFACE_LAPACK


  IMPLICIT NONE
  INTEGER,                         INTENT(IN)     :: D
  INTEGER,                         INTENT(INOUT)  :: INFO
  DOUBLE PRECISION,                INTENT(IN)     :: T1,T2
  COMPLEX*16, DIMENSION(D,D),      INTENT(OUT)     :: P_BARE
  TYPE(MODE),       DIMENSION(MODES_NUM),INTENT(OUT) :: FIELD

  COMPLEX*16,       DIMENSION(D,D) :: IDENTITY
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: U_Rot_1,U_Rot_2,U_RWA_1,U_RWA_2,U_TimeEvol

  INTEGER m
  DOUBLE PRECISION :: delta_1,delta_2,omega_1,omega_2,phi_1,phi_2,cos_2theta_1,cos_2theta_2,E_d_1,E_d_2,OMEGA_R_1,OMEGA_R_2


  IDENTITY = 0.0

  DO m=1,D
   IDENTITY(m,m) = 1.0
  END DO

  omega_1      = FIELD(2)%OMEGA
  phi_1        = FIELD(2)%phi_x
  OMEGA_R_1    = FIELD(2)%BX
  delta_1      = 1.0 - FIELD(2)%OMEGA
  E_d_1        = SQRT(delta_1**2+(OMEGA_R_1/2.0)**2)
  cos_2theta_1 = -delta_1/(E_d_1)


  omega_2      = FIELD(3)%OMEGA
  phi_2        = FIELD(3)%phi_z
  OMEGA_R_2    = FIELD(3)%BZ*sqrt(1-cos_2theta_1)
  delta_2      = E_D_1 - FIELD(3)%OMEGA
  E_d_2        = SQRT(delta_2**2+(OMEGA_R_2/2.0)**2)
  cos_2theta_2 = -delta_2/(E_d_2)

!  write(*,*) omega_2,sqrt((1+cos_2theta_1)/2.0)*sqrt((1+cos_2theta_2)/2.0), sqrt((1+cos_2theta_1)/2.0)*sqrt((1-cos_2theta_2)/2.0), &
!                &    sqrt((1-cos_2theta_1)/2.0)*sqrt((1-cos_2theta_2)/2.0),-sqrt((1-cos_2theta_1)/2.0)*sqrt((1+cos_2theta_2)/2.0), &
!                &    sqrt((1-cos_2theta_1)/2.0)*sqrt((1+cos_2theta_2)/2.0), sqrt((1-cos_2theta_1)/2.0)*sqrt((1-cos_2theta_2)/2.0), &
!                &   -sqrt((1+cos_2theta_1)/2.0)*sqrt((1-cos_2theta_2)/2.0), sqrt((1+cos_2theta_1)/2.0)*sqrt((1+cos_2theta_2)/2.0)


  ALLOCATE(U_ROT_1(D,D))
  ALLOCATE(U_RWA_1(D,D))

  ALLOCATE(U_ROT_2(D,D))
  ALLOCATE(U_RWA_2(D,D))

  ALLOCATE(U_TimeEvol(D,D))


  U_ROT_1 = 0.0
  U_ROT_2 = 0.0
  U_RWA_1 = 0.0
  U_RWA_2 = 0.0

  U_ROT_1(1,1) = exp(-DCMPLX(0.0,1.0)*omega_1*t2/2)
  U_ROT_1(2,2) = exp(DCMPLX(0.0,1.0)*omega_1*t2/2)

  U_ROT_2(1,1) = exp(-DCMPLX(0.0,1.0)*omega_2*t2/2)
  U_ROT_2(2,2) = exp(DCMPLX(0.0,1.0)*omega_2*t2/2)


  U_RWA_1(1,1)  = ( exp(-0.5*DCMPLX(0.0,1.0)*phi_1))*sqrt(0.5*(1-cos_2theta_1))
  U_RWA_1(1,2)  = ( exp(-0.5*DCMPLX(0.0,1.0)*phi_1))*sqrt(0.5*(1+cos_2theta_1))
  U_RWA_1(2,1)  = ( exp( 0.5*DCMPLX(0.0,1.0)*phi_1))*sqrt(0.5*(1+cos_2theta_1))
  U_RWA_1(2,2)  = (-exp( 0.5*DCMPLX(0.0,1.0)*phi_1))*sqrt(0.5*(1-cos_2theta_1))


  U_RWA_2(1,1)  = ( exp(-0.5*DCMPLX(0.0,1.0)*phi_2))*sqrt(0.5*(1-cos_2theta_2))
  U_RWA_2(1,2)  = ( exp(-0.5*DCMPLX(0.0,1.0)*phi_2))*sqrt(0.5*(1+cos_2theta_2))
  U_RWA_2(2,1)  = ( exp( 0.5*DCMPLX(0.0,1.0)*phi_2))*sqrt(0.5*(1+cos_2theta_2))
  U_RWA_2(2,2)  = (-exp( 0.5*DCMPLX(0.0,1.0)*phi_2))*sqrt(0.5*(1-cos_2theta_2))

  !P_BARE = MATMUL(U_ROT_1,U_RWA_1)

  !P_BARE = MATMUL(U_ROT_1,MATMUL(U_RWA_1,MATMUL(U_ROT_2,U_RWA_2)))

!write(*,*) abs(U_RWA_2(1,1)),cos_2theta_2,-delta_2,E_d_2
  ! Time evolution in the dressed basis of the double rotating frame
  U_TimeEvol = 0.0
  U_TimeEvol(1,1) = exp(-DCMPLX(0.0,1.0)*E_d_2*(t2-t1)/2.0)
  U_TimeEvol(2,2) = exp( DCMPLX(0.0,1.0)*E_d_2*(t2-t1)/2.0)

  ! Time evolution in the dressed basis
  U_TimeEvol = MATMUL(MATMUL(U_RWA_2,U_TimeEvol),TRANSPOSE(CONJG(U_RWA_2)))
  U_TimeEvol = MATMUL(U_ROT_2,U_TimeEvol)
  U_ROT_2(1,1) = exp(-DCMPLX(0.0,1.0)*omega_2*t1/2)
  U_ROT_2(2,2) = exp( DCMPLX(0.0,1.0)*omega_2*t1/2)
  U_TimeEvol = MATMUL(U_TimeEvol,transpose(conjg(U_ROT_2)))

  ! Time evolution in the bare basis
  U_TimeEvol = MATMUL(MATMUL(U_RWA_1,U_TimeEvol),TRANSPOSE(CONJG(U_RWA_1)))
  U_TimeEvol = MATMUL(U_ROT_1,U_TimeEvol)
  U_ROT_1(1,1) = exp(-DCMPLX(0.0,1.0)*omega_2*t1/2)
  U_ROT_1(2,2) = exp( DCMPLX(0.0,1.0)*omega_2*t1/2)
  U_TimeEvol = MATMUL(U_TimeEvol,transpose(conjg(U_ROT_1)))

  P_BARE = ABS(U_TimeEvol)**2

END SUBROUTINE EVOLUTIONOPERATOR_RWA_qubit
