SUBROUTINE MULTIMODEFLOQUETMATRIX_SP_C(ATOM__C,NM,NF,MODES_NUM,FIELDS_C,INFO)!VALUES_,ROW_INDEX_,COLUMN_,SP,INFO)
! THIS SUBROUTINE BUILDS THE MULTIMODE FLOQUET MATRIX

!ATOM_      (IN)    : type of quantum system
!NM         (IN)    : number of modes
!NF         (IN)    : number of driving fields
!MODES_NUM  (IN)    : vector indicating the number of harmonics of each driving field (mode)
!FIELDS_C     (IN)    : Fields

! THE FOLLOWING VARIABLES ARE DECLARED AS GLOBAL ALLOCATABLE ARRAYS. THIS SUBROUTINE SET THEIR VALUES AND SIZE.

!VALUES_    (OUT)   : Hamiltonian values
!ROW_INDEX_ (OUT)   : vector indicating the row position of values
!COLUMN_    (OUT)   : vector indicating the column position of the values
!INFO       (INOUT) : error flag. INFO=0 means there is no error

  USE TYPES_C             !(modes.f90)
  USE MERGINGARRAYS     !(utils.f90)
  USE SPARSE_INTERFACE  !(sparse_utils.f90)
  USE MODES_4F          ! DEFINED IN modes_C.f90, declares atom_,coupling, values__, row_index__, column__
  
  IMPLICIT NONE
  INTEGER                  ,            INTENT(IN)    :: NM,NF
  TYPE(MODE_C), DIMENSION(NF),          INTENT(INout) :: FIELDS_C
  TYPE(ATOM_C),                         INTENT(IN)    :: ATOM__C
  INTEGER,    DIMENSION(NM),            INTENT(IN)    :: MODES_NUM
  INTEGER,                              INTENT(INOUT) :: INFO
  
!  COMPLEX*16, DIMENSION(:), ALLOCATABLE  :: VALUES__
!  INTEGER,    DIMENSION(:), ALLOCATABLE  :: COLUMN__
!  INTEGER,    DIMENSION(:), ALLOCATABLE  :: ROW_INDEX__
  !INTEGER,                              INTENT(OUT)   :: SP
 
!    write(*,*) 'FORTRAN FLOQUETMATRIX_SP SAYS',NM,NF,MODES_NUM!, COULPLIG(3)%OMEGA
  IF (INFO.EQ.0) THEN      
    CALL MULTIMODEFLOQUETMATRIX_SP(ATOM_,NM,NF,MODES_NUM,COUPLING,VALUES__,ROW_INDEX__,COLUMN__,INFO)
!    WRITE(*,*) "FORTRAN MULTIMODEHAMILTONAISPC SAYS: SIZE(VALUES__,1) =)",SIZE(VALUES__,1),SIZE(ROW_INDEX__,1)  
    H_FLOQUET_SIZE = SIZE(ROW_INDEX__,1)-1
  END IF

END SUBROUTINE MULTIMODEFLOQUETMATRIX_SP_C ! _SP  sparse packing


SUBROUTINE MULTIMODEFLOQUETMATRIX_PYTHON_SP_C(ATOM__C,NM,NF,MODES_NUM,FIELDS_C,MMF_DIM,INFO)! RESULT(MMF_DIM)
! THIS SUBROUTINE BUILDS THE MULTIMODE FLOQUET MATRIX

!ATOM_      (IN)    : type of quantum system
!NM         (IN)    : number of modes
!NF         (IN)    : number of driving fields
!MODES_NUM  (IN)    : vector indicating the number of harmonics of each driving field (mode)
!FIELDS_C     (IN)    : Fields

! THE FOLLOWING VARIABLES ARE DECLARED AS GLOBAL ALLOCATABLE ARRAYS. THIS SUBROUTINE SET THEIR VALUES AND SIZE.

!VALUES_    (OUT)   : Hamiltonian values
!ROW_INDEX_ (OUT)   : vector indicating the row position of values
!COLUMN_    (OUT)   : vector indicating the column position of the values
!INFO       (INOUT) : error flag. INFO=0 means there is no error

  USE TYPES_C             !(modes.f90)
  USE MERGINGARRAYS     !(utils.f90)
  USE SPARSE_INTERFACE  !(sparse_utils.f90)
  USE MODES_4F          ! DEFINED IN modes_C.f90, declares atom_,coupling, values__, row_index__, column__
  
  IMPLICIT NONE
  INTEGER                  ,            INTENT(IN)    :: NM,NF
  TYPE(MODE_C), DIMENSION(NF),          INTENT(INout) :: FIELDS_C
  TYPE(ATOM_C),                         INTENT(IN)    :: ATOM__C
  INTEGER,    DIMENSION(NM),            INTENT(IN)    :: MODES_NUM
  INTEGER,    DIMENSION(3),             INTENT(OUT)   :: MMF_DIM
  INTEGER,                              INTENT(INOUT) :: INFO
  
  !INTEGER, DIMENSION(3) :: MMF_DIM
 ! INTEGER MMF_DIM
  IF (INFO.EQ.0) THEN      
    CALL MULTIMODEFLOQUETMATRIX_SP(ATOM_,NM,NF,MODES_NUM,COUPLING,VALUES__,ROW_INDEX__,COLUMN__,INFO)
    H_FLOQUET_SIZE    = SIZE(ROW_INDEX__,1)-1
!    MMF_DIM = H_FLOQUET_SIZE
    MMF_DIM(1)        = SIZE(VALUES__,1)
    MMF_DIM(2)        = SIZE(ROW_INDEX__,1)
    MMF_DIM(3)        = SIZE(COLUMN__,1)
  END IF
END SUBROUTINE MULTIMODEFLOQUETMATRIX_PYTHON_SP_C
!END FUNCTION  MULTIMODEFLOQUETMATRIX_PYTHON_SP_C ! _SP  sparse packing

SUBROUTINE GET_H_FLOQUET_SP_C(h_floquet_size_,VALUES,ROW_INDEX,COLUMN,INFO)
!SUBROUTINE GET_H_FLOQUET_SP_C(h_floquet_size_,VALUES)!, VALUES,ROW_INDEX,COLUMN,INFO)
  
  USE MODES_4F          ! DEFINED IN modes_C.f90, declares atom_,coupling, values__, row_index__, column__

  IMPLICIT NONE
  INTEGER, DIMENSION(3), INTENT(IN) :: h_floquet_size_
  COMPLEX*16,DIMENSION(h_floquet_size_(1)), intent(OUT)   :: VALUES
  INTEGER,   DIMENSION(h_floquet_size_(2)), intent(OUT)   :: ROW_INDEX
  INTEGER,   DIMENSION(h_floquet_size_(3)), intent(OUT)   :: COLUMN
  INTEGER,                                  INTENT(INOUT) :: INFO

  !write(*,*) h_floquet_size_
  !write(*,*) values__
  !write(*,*) size(values__,1),h_floquet_size_

  VALUES    = VALUES__
  ROW_INDEX = ROW_INDEX__
  COLUMN    = COLUMN__


END SUBROUTINE GET_H_FLOQUET_SP_C
