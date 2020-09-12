MODULE FEAST
  integer     fpm(128)
  real*8      Emin,Emax
  real*8      epsout
  integer     loop
  integer     M0 ! initial guess 
  integer     M1 ! total number of eigenvalues found
  integer     info_FEAST
  real*8,     DIMENSION(:),   ALLOCATABLE :: E, RES ! vector of eigenvalues
  complex*16, DIMENSION(:,:), ALLOCATABLE :: X      ! matrix with eigenvectore
END MODULE FEAST

MODULE TYPES

  TYPE :: MODE
     DOUBLE PRECISION :: OMEGA
     COMPLEX*16       :: X,Y,Z
     DOUBLE PRECISION :: phi_x,phi_y,phi_z
     INTEGER          :: N_Floquet
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: V
     COMPLEX*16, DIMENSION(:),   ALLOCATABLE :: VALUES
     INTEGER,    DIMENSION(:),   ALLOCATABLE :: ROW,COLUMN
  END TYPE MODE
  
  TYPE :: ATOM
     INTEGER          :: id_system
     INTEGER          :: D_BARE
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_BARE
  END TYPE ATOM

  TYPE :: HARMONIC_FACTORS
     COMPLEX*16,DIMENSION(:,:), ALLOCATABLE :: U,U_r,U_AVG
     INTEGER,   DIMENSION(:),   ALLOCATABLE :: n
  END type HARMONIC_FACTORS

!!$  TYPE :: MWCOUPLING
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: TOP
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: TOP_DAGGER
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: DC
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: DC_DAGGER
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: MW
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: MW_DAGGER
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: RF
!!$     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: RF_DAGGER
!!$  END type MWCOUPLING
END MODULE TYPES

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


MODULE QUICKSORTINTERFACE
  IMPLICIT NONE
  PUBLIC
  INTERFACE
     recursive subroutine QUICK_SORT_I_T(a,idx_a,na)
       implicit none 
       
       ! DUMMY ARGUMENTS
       integer,                         intent(in)    :: na       ! nr or items to sort
       DOUBLE PRECISION, dimension(nA), intent(inout) :: a     ! vector to be sorted
       integer,          dimension(nA), intent(inout) :: idx_a ! sorted indecies of a
     end subroutine QUICK_SORT_I_T
     
!     recursive subroutine QUICK_SORT_INTEGERS(a,idx_a,na)
!     implicit none 

     ! DUMMY ARGUMENTS
!     integer,                intent(in) :: na       ! nr or items to sort
!     integer, dimension(nA), intent(inout) :: a     ! vector to be sorted
!     integer, dimension(nA), intent(inout) :: idx_a ! sorted indecies of a
     
!     end subroutine quick_sort_integers
 
     subroutine InsertionSort_I(a,idx_a,na)

     !use sangoma_base, only: REALPREC, INTPREC
 
     implicit none

     ! DUMMY ARGUMENTS
     integer,                         intent(in)    :: na
     integer, dimension(nA), intent(inout) :: a
     integer,          dimension(nA), intent(inout) :: idx_a
     end subroutine InsertionSort_I
     
     subroutine InsertionSort_d(a,idx_a,na)

     !use sangoma_base, only: REALPREC, INTPREC
 
     implicit none

     ! DUMMY ARGUMENTS
     integer,                         intent(in)    :: na
     double precision, dimension(nA), intent(inout) :: a
     integer,          dimension(nA), intent(inout) :: idx_a
    
     end subroutine InsertionSort_d
 

  END INTERFACE
END MODULE QUICKSORTINTERFACE
SUBROUTINE WRITE_MATRIX(A)
! it writes a matrix of doubles nxm on the screen
  DOUBLE PRECISION, DIMENSION(:,:) :: A
  CHARACTER(LEN=105) STRING
  CHARACTER(LEN=105) aux_char
  integer :: aux

  aux = int(UBOUND(A,2))
  !write(*,*) aux
  write(aux_char,"(I4)") aux
  aux_char = trim(aux_char)
  write(string,"(A1,I4,A6)") "(",aux,"E15.6)"

  DO I = LBOUND(A,1), UBOUND(A,1)
     WRITE(*,string) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
  END DO
  WRITE(*,*)
  WRITE(*,*)
END SUBROUTINE WRITE_MATRIX

SUBROUTINE WRITE_MATRIX_INT(A)
!it writes a matrix of integer nxm on the screen
  INTEGER, DIMENSION(:,:) :: A
  WRITE(*,*)
  DO I = LBOUND(A,1), UBOUND(A,1)
     WRITE(*,*) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
  END DO
END SUBROUTINE WRITE_MATRIX_INT


SUBROUTINE COORDINATEPACKING(D,A,V,R,C,index,INFO)
  IMPLICIT NONE
  INTEGER,INTENT(IN):: D
  COMPLEX*16,DIMENSION(D,D),INTENT(IN)  :: A
  COMPLEX*16,DIMENSION(D*D),INTENT(OUT) :: V
  INTEGER, DIMENSION(D*D),  INTENT(OUT) :: R,C
  INTEGER, INTENT(OUT)   :: index
  INTEGER, INTENT(INOUT) :: INFO
  
  INTEGER I,J
  V=0
  R=0
  C=0
  
  index = 1
  DO I=1,D
     DO J=1,D
        IF(ABS(A(I,J)).GT.0) THEN
           V(index) = A(I,J)
           R(index) = I
           C(index) = J
           index = index+1
        END IF
     END DO
  END DO
  index = index-1
END SUBROUTINE COORDINATEPACKING

MODULE MERGINGARRAYS
  INTERFACE
     SUBROUTINE APPENDARRAYS(V,B,INFO)
       COMPLEX*16, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
       COMPLEX*16, DIMENSION(:),INTENT(IN)    :: B
       INTEGER,                 INTENT(INOUT) :: INFO
     END SUBROUTINE APPENDARRAYS
     SUBROUTINE APPENDARRAYSI(V,B,INFO)
       INTEGER, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
       INTEGER, DIMENSION(:),INTENT(IN)    :: B
       INTEGER,                 INTENT(INOUT) :: INFO
     END SUBROUTINE APPENDARRAYSI
  END INTERFACE
END MODULE MERGINGARRAYS

SUBROUTINE APPENDARRAYS(V,B,INFO)
  COMPLEX*16, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
  COMPLEX*16, DIMENSION(:),INTENT(IN)    :: B
  INTEGER,                 INTENT(INOUT) :: INFO
  
  COMPLEX*16,DIMENSION(:),ALLOCATABLE :: tmp_arr
!  write(*,*) V
!  write(*,*) B
  ALLOCATE(tmp_arr(SIZE(V,1)+SIZE(B,1)))
  tmp_arr(1:SIZE(V,1))=V
  tmp_arr(SIZE(V,1)+1:SIZE(tmp_arr))=B
  DEALLOCATE(V)
  ALLOCATE(V(SIZE(tmp_arr)))
  V=tmp_arr
END SUBROUTINE APPENDARRAYS

SUBROUTINE APPENDARRAYSI(V,B,INFO)
  INTEGER, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
  INTEGER, DIMENSION(:),INTENT(IN)    :: B
  INTEGER,                 INTENT(INOUT) :: INFO
  
  COMPLEX*16,DIMENSION(:),ALLOCATABLE :: tmp_arr
  
  ALLOCATE(tmp_arr(SIZE(V,1)+SIZE(B,1)))
  tmp_arr(1:SIZE(V,1))=V
  tmp_arr(SIZE(V,1)+1:SIZE(tmp_arr))=B
  DEALLOCATE(V)
  ALLOCATE(V(SIZE(tmp_arr)))
  V=tmp_arr
END SUBROUTINE APPENDARRAYSI


MODULE funciones
  IMPLICIT NONE

CONTAINS

  FUNCTION delta_kr(M,N)
    IMPLICIT NONE
    DOUBLE PRECISION :: M,N
    INTEGER DELTA_KR    
    IF(M.EQ.N) DELTA_KR = 1
    IF(M.NE.N) DELTA_KR = 0    
  END FUNCTION delta_kr

  FUNCTION delta_kr_int(M,N)    
    IMPLICIT NONE
    INTEGER :: M,N
    INTEGER DELTA_KR_INT    
    IF(M.EQ.N) DELTA_KR_INT = 1
    IF(M.NE.N) DELTA_KR_INT = 0    
  END FUNCTION delta_kr_int

END MODULE funciones




MODULE ATOMIC_PROPERTIES
  USE physical_constants
  IMPLICIT NONE
  ! Parameters values by default
  DOUBLE PRECISION :: L = 0.0
  DOUBLE PRECISION :: S = 0.5
  DOUBLE PRECISION :: J,F,gf,mf  
  DOUBLE PRECISION :: gF_2,gF_1,G_F
  INTEGER          :: Total_states_LSI
  INTEGER          :: Fup,Fdown
  DOUBLE PRECISION :: mass_at,Ftotal,I,g_I,g_J,A,a_s,alpha_E
  CHARACTER(LEN=7) :: ID_name
   
  !87Rb
  DOUBLE PRECISION, PARAMETER :: mass_at_87Rb = 87*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_87Rb  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_87Rb       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_87Rb       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_87Rb     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_87Rb     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_87Rb       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_87Rb     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_87Rb = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_87Rb     =  2
  INTEGER,          PARAMETER :: Fdown_87Rb   =  1
  CHARACTER(LEN=7) :: ID_name_87Rb = "87Rb"
 
  !85Rb
  DOUBLE PRECISION, PARAMETER :: mass_at_85Rb = 85*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_85Rb  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_85Rb       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_85Rb       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_85Rb     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_85Rb     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_85Rb       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_85Rb     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_85Rb = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_85Rb     =  2
  INTEGER,          PARAMETER :: Fdown_85Rb   =  1
  CHARACTER(LEN=7) :: ID_name_85Rb = "85Rb"


  !6Li  
  DOUBLE PRECISION, PARAMETER :: mass_at_6Li = 6*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6Li  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_6Li       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6Li       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6Li     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6Li     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6Li       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6Li     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6Li = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6Li     =  2
  INTEGER,          PARAMETER :: Fdown_6Li   =  1
  CHARACTER(LEN=7) :: ID_name_6Li = "6Li"


  !Cs
  DOUBLE PRECISION, PARAMETER :: mass_at_6Cs = 87*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6Cs  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_Cs       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6Cs       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6Cs     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6Cs     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6Cs       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6Cs     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6Cs = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6Cs     =  2
  INTEGER,          PARAMETER :: Fdown_6Cs   =  1
  CHARACTER(LEN=7) :: ID_name_Cs = "Cs"


  !41K
  DOUBLE PRECISION, PARAMETER :: mass_at_6K = 41*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6K  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_K       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6K       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6K     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6K     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6K       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6K     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6K = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6K     =  2
  INTEGER,          PARAMETER :: Fdown_6K   =  1
  CHARACTER(LEN=7) :: ID_name_41K = "41K"

  !Na
  DOUBLE PRECISION, PARAMETER :: mass_at_6Na = 87*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6Na  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_Na        =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6Na       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6Na     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6Na     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6Na       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6Na     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6Na = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6Na     =  2
  INTEGER,          PARAMETER :: Fdown_6Na   =  1
  CHARACTER(LEN=7) :: ID_name_Na = "Na"


  !Qubit
  DOUBLE PRECISION, PARAMETER :: mass_at_qubit =  amu
  DOUBLE PRECISION :: Ftotal_qubit  =  0.0
  DOUBLE PRECISION :: J_qubit       =  0.0
  DOUBLE PRECISION :: I_qubit       =  0.0
  DOUBLE PRECISION, PARAMETER :: g_J_qubit     =  1.0
  DOUBLE PRECISION, PARAMETER :: g_I_qubit     =  0.0
  DOUBLE PRECISION, PARAMETER :: A_qubit       =  1.0
  DOUBLE PRECISION, PARAMETER :: a_s_qubit     =  0.0
  DOUBLE PRECISION, PARAMETER :: alpha_E_qubit =  0.0
  INTEGER,          PARAMETER :: Fup_qubit     =  0
  INTEGER,          PARAMETER :: Fdown_qubit   =  0.5
  CHARACTER(LEN=7) :: ID_name_qubit = "qubit"

  !spin
  DOUBLE PRECISION :: I_spin   =  0.0
  DOUBLE PRECISION :: J_spin   =  0.0  
  DOUBLE PRECISION :: gJ_spin  =  1.0
  DOUBLE PRECISION :: gI_spin  =  0.0
  DOUBLE PRECISION :: A_spin   =  1.0
  DOUBLE PRECISION :: a_s_spin =  0.0
  DOUBLE PRECISION :: alpha_E_spin = 0.0
  INTEGER          :: Fup_spin     =  1
  INTEGER          :: Fdown_spin   =  1
  CHARACTER(LEN=7) :: ID_name_spin = "spin"


  !lattice
  CHARACTER        :: PERIODIC      
  CHARACTER(LEN=7) :: ID_name_lattice = "lattice"

END MODULE ATOMIC_PROPERTIES

