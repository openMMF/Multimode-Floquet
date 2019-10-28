MODULE physical_constants
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: pi           = 4.0*ATAN(1.0)
  DOUBLE PRECISION, PARAMETER :: e            = 1.602176462E-19
  DOUBLE PRECISION, PARAMETER :: h_P          = 6.62606957E-34
  DOUBLE PRECISION, PARAMETER :: hbar         = 1.054571800E-34
  DOUBLE PRECISION, PARAMETER :: mu_B         = 9.274009994E-24
  DOUBLE PRECISION, PARAMETER :: k_B          = 1.3806488E-23
  DOUBLE PRECISION, PARAMETER :: mu_cero      = 12.566370614E-7
  DOUBLE PRECISION, PARAMETER :: epsilon_cero = 8.854187817E-12 
  DOUBLE PRECISION, PARAMETER :: amu          = 1.660538921E-27
  DOUBLE PRECISION, PARAMETER :: g_t          = 9.8
  DOUBLE PRECISION, PARAMETER :: SB_ct        = 5.6704E-8
  COMPLEX*16,       PARAMETER :: J_IMAG       = DCMPLX(0.0,1.0)
  DOUBLE PRECISION, PARAMETER :: speedoflight = 299792458.0
  DOUBLE PRECISION            :: TOTAL_TIME
END MODULE physical_constants

MODULE ARRAYS

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Identity,j_x,j_y,j_z,I_x,I_y,I_z 
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: HAMILTONIAN, H_FLOQUET!H_hyperfine
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_IJ,U_ZEEMAN!Z_M,H_OLD,U_ZEEMAN,Z_M_SUBSET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_FLOQUET_COPY
  !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CLEBSH_GORDAN_JtoF,h_hyperfine
  !COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_RF_DAGGER,H_ALPHA_DAGGER,H_ALPHA,H_FLOQUET_COPY,H_MW,H_AUX,U_RF
  !COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: observable, observable_extended, MW_coupling_dressedbasis
  !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: W_SPACE,W_SPACEF,W_SPACEF_0,E_OLD
  !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Fx,Fy,Fz,g_F_matrix
  !COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: Hamiltonian_F,Identity_F,H_AUX,U_RF
  !COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_FLOQUET_INTERACTION,H_FLOQUET_INTERACTION_DAGGER,H_FLOQUET_2D,H_FLOQUETBAND
  !INTEGER,          DIMENSION(:,:), ALLOCATABLE :: F_t,H_w,H_J,H_M,Jz_dash,Fz_dash
  !INTEGER,          DIMENSION(:),   ALLOCATABLE :: index_state
  !INTEGER                                       :: KD
  !DOUBLE PRECISION, DIMENSION(3)                :: POSITION,DELTA_POSITION

END MODULE ARRAYS
 
!!$MODULE DCRFMWFIELDS
!!$ IMPLICIT NONE
!!$  ! DC, RF and MW parameters
!!$  !! DC FIELD
!!$  DOUBLE PRECISION B_DC_X, B_DC_Y, B_DC_Z   
!!$
!!$  !! RF FIELD
!!$  DOUBLE PRECISION  OMEGA_RF, PHI_RF_X,PHI_RF_Y,PHI_RF_Z
!!$  DOUBLE PRECISION  B_RF_X,B_RF_Y,B_RF_Z
!!$
!!$  !! MW FIELD
!!$  DOUBLE PRECISION  OMEGA_MW, PHI_MW_X,PHI_MW_Y,PHI_MW_Z
!!$  DOUBLE PRECISION  B_MW_X,B_MW_Y,B_MW_Z
!!$
!!$END MODULE DCRFMWFIELDS

MODULE subinterface
  IMPLICIT NONE
  INTERFACE
     SUBROUTINE I_and_J_representations(j_x,j_y,j_z,I_x,I_y,I_z,L,S,I)
       
       IMPLICIT  NONE
       DOUBLE PRECISION, DIMENSION(:,:),INTENT(INOUT) :: j_x,j_y,j_z,I_x,I_y,I_z
       DOUBLE PRECISION, INTENT(IN) :: L,S,I
     END SUBROUTINE I_and_J_representations

     SUBROUTINE F_representation(Fx,Fy,Fz,Ftotal_)       
       
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT):: Fx,Fy,Fz
       DOUBLE PRECISION, INTENT(IN) :: Ftotal_
     END SUBROUTINE F_representation
          

     SUBROUTINE SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,JTOTAL,ID,INFO)
     !SUBROUTINE SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,JTOTAL,JTOTAL_HALF,ID,INFO)
       ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
       ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
       ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
       !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES
       USE TYPES
       IMPLICIT NONE
       CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: ATOMICSPECIE
       CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: MANIFOLD  !
       !INTEGER,          OPTIONAL, INTENT(IN) :: JTOTAL
       DOUBLE PRECISION,          OPTIONAL, INTENT(IN) :: JTOTAL
       TYPE(ATOM),OPTIONAL,INTENT(OUT) :: ID       
       INTEGER, INTENT(INOUT) :: INFO
     END SUBROUTINE SET_ATOMIC_PARAMETERS

     SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS(ID,DRESSINGFIELDS_INDICES,MODES_NUM,FIELDS, U_FD,E_DRESSED,INFO)
       ! ID        (in)    :: TYPE(ATOM) system ID
       ! DRESSINGFIELDS_INDICES (in) :: integer array indicating the indices of the dressing modes
       ! MODES_NUM (in)    :: integer array indicating the number of harmonics of all driving modes 
       ! FIELDS    (in)    :: Array of TYPE(MODE) of dimension 
       ! U_FD      (out)   :: complex*16 matrix fourier decomposition of the micromotion operator of the dressed basis
       ! E_DRESSED (out)   :: dressed energies
       ! INFO      (inout) :: error flag
       USE TYPES

       TYPE(ATOM),                     INTENT(IN)  :: ID
       INTEGER,    DIMENSION(:),       INTENT(IN)  :: DRESSINGFIELDS_INDICES
       INTEGER,    DIMENSION(:),       INTENT(IN)  :: MODES_NUM
       TYPE(MODE), DIMENSION(:),       INTENT(IN)  :: FIELDS
       COMPLEX*16, DIMENSION(:,:),     INTENT(OUT) :: U_FD
       DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: E_DRESSED
     END SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS

     SUBROUTINE MICROMOTIONDRESSEDBASIS(ID,MODES_NUM,DRESSINGFIELDS_INDICES,&
          & FIELDS,U_F_MODES,E_MULTIFLOQUET,T1,U,INFO) 
       
       ! ID (in)        :: TYPE(ATOM) system ID
       ! MODES_NUM (in) :: integer array indicating the number of harmonics of each driving mode
       ! DRESSINFIELDS_INDICES :: integer array indicating the indices of the dressing modes
       ! FIELDS         :: Array of TYPE(MODES) with NM components (all driving fields)
       ! U_F_MODES      :: complex*16 matrix of dimension DxD. Fourier decomposition of the micromotion operator of the dressed basis
       ! E_MULTIFLOQUET :: dressed energies
       ! T1             :: double precision, time
       ! U              :: complex*16 matrix of dimension D_BARE x D_BARE. micromotion operator at time T1
       ! INFO           :: error flag
       
       
       USE TYPES
       IMPLICIT NONE
       TYPE(ATOM),                       INTENT(IN)    :: ID
       INTEGER,          DIMENSION(:),   INTENT(IN)    :: MODES_NUM
       INTEGER,          DIMENSION(:),   INTENT(IN)    :: DRESSINGFIELDS_INDICES
       COMPLEX*16,       DIMENSION(:,:), INTENT(IN)    :: U_F_MODES
       DOUBLE PRECISION, DIMENSION(:),   INTENT(IN)    :: E_MULTIFLOQUET
       TYPE(MODE),       DIMENSION(:),   INTENT(IN)    :: FIELDS
       DOUBLE PRECISION ,                INTENT(IN)    :: T1
       COMPLEX*16,       DIMENSION(:,:), INTENT(OUT)   :: U
       INTEGER,                          INTENT(INOUT) :: INFO
     END SUBROUTINE MICROMOTIONDRESSEDBASIS
  END INTERFACE
END MODULE subinterface


MODULE SUBINTERFACE_LAPACK

  IMPLICIT NONE
  PUBLIC 
  INTERFACE
!     SUBROUTINE LAPACK_FULLEIGENVALUES(H,N,W_SPACE,INFO)
!       IMPLICIT NONE
!       INTEGER,                        INTENT(IN)    :: N
!       COMPLEX*16, DIMENSION(:,:),     INTENT(INOUT) :: H  
!       DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)   :: W_SPACE
!       INTEGER,                        INTENT(INOUT) :: INFO
!       
!     END SUBROUTINE LAPACK_FULLEIGENVALUES
!     
!     SUBROUTINE LAPACK_SELECTEIGENVALUES(H,N,W_SPACE,L1,L2,Z_M,INFO)
!       
!       INTEGER,                        INTENT(IN)    :: N,L1,L2
!       COMPLEX*16, DIMENSION(:,:),     INTENT(INOUT) :: H
!       COMPLEX*16, DIMENSION(:,:),     INTENT(OUT)   :: Z_M
!       DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)   :: W_SPACE
!       INTEGER,                        INTENT(OUT)   :: INFO
!     END SUBROUTINE LAPACK_SELECTEIGENVALUES

     SUBROUTINE WRITE_MATRIX(A)
       DOUBLE PRECISION, DIMENSION(:,:) :: A
     END SUBROUTINE  WRITE_MATRIX

     SUBROUTINE WRITE_MATRIX_INT(A)
       INTEGER, DIMENSION(:,:) :: A
     END SUBROUTINE  WRITE_MATRIX_INT
     
  END INTERFACE
END MODULE SUBINTERFACE_LAPACK


