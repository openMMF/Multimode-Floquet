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
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: HAMILTONIAN, H_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_IJ,U_ZEEMAN
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_FLOQUET_COPY

END MODULE ARRAYS
 
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
