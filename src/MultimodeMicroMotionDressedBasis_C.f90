!!$MODULE PASSING2C
!!$  INTERFACE
!!$     SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS_C(ID,DRESSINGFIELDS_INDICES,ND,U_FD,E_DRESSED,INFO)
!!$       ! ID        (in)    :: TYPE(ATOM) system ID
!!$       ! DRESSINGFIELDS_INDICES (in) :: integer array indicating the indices of the dressing modes
!!$       ! MODES_NUM (in)    :: integer array indicating the number of harmonics of all driving modes 
!!$       ! FIELDS    (in)    :: Array of TYPE(MODE) of dimension 
!!$       ! U_FD      (out)   :: complex*16 matrix fourier decomposition of the micromotion operator of the dressed basis
!!$       ! E_DRESSED (out)   :: dressed energies
!!$       ! INFO      (inout) :: error flag
!!$       USE TYPES_C
!!$       USE SUBINTERFACE
!!$       USE MODES_4F
!!$       USE TYPES
!!$       IMPLICIT NONE
!!$       
!!$       TYPE(ATOM_C),                       INTENT(IN)  :: ID
!!$       INTEGER,                            INTENT(IN)  :: ND!DF,NM,ND,NF
!!$       INTEGER,          DIMENSION(:),     INTENT(IN)  :: DRESSINGFIELDS_INDICES
!!$       ! INTEGER,          DIMENSION(:),     INTENT(IN)  :: MODES_NUM
!!$       ! TYPE(MODE_C),     DIMENSION(:),     INTENT(IN)  :: FIELDS
!!$       COMPLEX*16,       DIMENSION(ND,ND), INTENT(OUT) :: U_FD
!!$       DOUBLE PRECISION, DIMENSION(ND),    INTENT(OUT) :: E_DRESSED
!!$       INTEGER,                            INTENT(INOUT) :: INFO
!!$     END SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS_C 
!!$  END INTERFACE
!!$END MODULE PASSING2C

SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS_C(ID,DF,DRESSINGFIELDS_INDICES,NM,MODES_NUM,NF,FIELDS,ND,U_FD,E_DRESSED,INFO)
!SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS_C(ID,DRESSINGFIELDS_INDICES,MODES_NUM,FIELDS,ND,U_FD,E_DRESSED,INFO)
!SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS_C(ID,DRESSINGFIELDS_INDICES,ND,U_FD,E_DRESSED,INFO)
! ID        (in)    :: TYPE(ATOM) system ID
! DRESSINGFIELDS_INDICES (in) :: integer array indicating the indices of the dressing modes
! MODES_NUM (in)    :: integer array indicating the number of harmonics of all driving modes 
! FIELDS    (in)    :: Array of TYPE(MODE) of dimension 
! U_FD      (out)   :: complex*16 matrix fourier decomposition of the micromotion operator of the dressed basis
! E_DRESSED (out)   :: dressed energies
! INFO      (inout) :: error flag
  USE TYPES_C
  USE SUBINTERFACE
  USE MODES_4F
  USE TYPES
  IMPLICIT NONE

  TYPE(ATOM_C),                       INTENT(IN)  :: ID
  INTEGER,                            INTENT(IN)  :: ND,DF,NM,NF
  INTEGER,          DIMENSION(DF),     INTENT(IN)  :: DRESSINGFIELDS_INDICES
  INTEGER,          DIMENSION(NM),     INTENT(IN)  :: MODES_NUM
  TYPE(MODE_C),     DIMENSION(NF),     INTENT(IN)  :: FIELDS
  COMPLEX*16,       DIMENSION(ND,ND), INTENT(OUT) :: U_FD
  DOUBLE PRECISION, DIMENSION(ND),    INTENT(OUT) :: E_DRESSED
  INTEGER,                            INTENT(INOUT) :: INFO

  WRITE(*,*) dressingfields_indices
  CALL MICROMOTIONFOURIERDRESSEDBASIS(ATOM_,DRESSINGFIELDS_INDICES,MODES_NUM,COUPLING, U_FD__,E_DRESSED__,INFO)
!  U_FD      = U_FD__
!  E_DRESSED = E_DRESSED__

END SUBROUTINE MICROMOTIONFOURIERDRESSEDBASIS_C

SUBROUTINE MICROMOTIONDRESSEDBASIS_C(ID,MODES_NUM,DRESSINGFIELDS_INDICES,FIELDS,T1,U,INFO) 

! ID (in)        :: TYPE(ATOM) system ID
! MODES_NUM (in) :: integer array indicating the number of harmonics of each driving mode
! DRESSINFIELDS_INDICES :: integer array indicating the indices of the dressing modes
! FIELDS         :: Array of TYPE(MODES) with NM components (all driving fields)
! U_F_MODES      :: complex*16 matrix of dimension DxD. Fourier decomposition of the micromotion operator of the dressed basis
! E_MULTIFLOQUET :: dressed energies
! T1             :: double precision, time
! U              :: complex*16 matrix of dimension D_BARE x D_BARE. micromotion operator at time T1
! INFO           :: error flag


  USE TYPES_C
  USE SUBINTERFACE
  USE MODES_4F
  USE TYPES
  IMPLICIT NONE
  TYPE(ATOM),                       INTENT(IN)    :: ID
  INTEGER,          DIMENSION(:),   INTENT(IN)    :: MODES_NUM
  INTEGER,          DIMENSION(:),   INTENT(IN)    :: DRESSINGFIELDS_INDICES
!  COMPLEX*16,       DIMENSION(:,:), INTENT(IN)    :: U_F_MODES
!  DOUBLE PRECISION, DIMENSION(:),   INTENT(IN)    :: E_MULTIFLOQUET
  TYPE(MODE),       DIMENSION(:),   INTENT(IN)    :: FIELDS
  DOUBLE PRECISION ,                INTENT(IN)    :: T1
  COMPLEX*16,       DIMENSION(:,:), INTENT(OUT)   :: U
  INTEGER,                          INTENT(INOUT) :: INFO
  

  CALL MICROMOTIONDRESSEDBASIS(ID,MODES_NUM,DRESSINGFIELDS_INDICES,FIELDS,U_FD__,E_DRESSED__,T1,U,INFO) 

END SUBROUTINE MICROMOTIONDRESSEDBASIS_C


