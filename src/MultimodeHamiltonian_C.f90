SUBROUTINE MULTIMODEFLOQUETMATRIX_C(ATOM__C,NM,NF,MODES_NUM,FIELD_C,INFO)
  !ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO
  !  USE FLOQUET
  !ATOM_C type atom, -> dimension of the bare Hilbert space
  !NM -> number of modes
  !NF -> Number of Fields
  !MODES_NUM -> number of harmonics of each mode
  !FIELD_C -> Field couplings
  !INFO

  USE ARRAYS
  USE TYPES_C
  USE MODES_4F

  IMPLICIT NONE
  INTEGER,                     INTENT(IN)    :: NM,NF
  INTEGER,                     INTENT(INOUT) :: INFO
  INTEGER,      DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE_C), DIMENSION(NF), INTENT(IN)    :: FIELD_C
  TYPE(ATOM_C),                INTENT(IN)    :: ATOM__C                     

!  write(*,*) "C WRAPPER:", info,NM,NF,MODES_NUM
  IF (INFO.EQ.0) THEN
 !    write(*,*) "C wrapper"
     if(info.eq.0) CALL MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,COUPLING,INFO)
     H_FLOQUET_SIZE = SIZE(H_FLOQUET,1)
!     write(*,*) "C wrapper"
     !write(*,*) h_floquet_size,size(H_floquet,1),coupling(1)%V(1,1)!, atom_%d_bare,h_floquet(1,1)
     !CALL WRITE_MATRIX(H_FLOQUET)
  ELSE
     WRITE(*,*) "THERE IS AN ERROR UPSTREAM :INFO.NE.0", INFO
  END IF
END SUBROUTINE MULTIMODEFLOQUETMATRIX_C


SUBROUTINE MULTIMODEFLOQUETMATRIX_C_PYTHON(ATOM__C,NM,NF,MODES_NUM,FIELD_C,MMF_DIM,INFO)
  !ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO
  !  USE FLOQUET
  !ATOM_C type atom, -> dimension of the bare Hilbert space
  !NM -> number of modes
  !NF -> Number of Fields
  !MODES_NUM -> number of harmonics of each mode
  !FIELD_C -> Field couplings
  !INFO

  USE ARRAYS
  USE TYPES_C
  USE MODES_4F

  IMPLICIT NONE
  INTEGER,                     INTENT(IN)    :: NM,NF
  INTEGER,                     INTENT(INOUT) :: INFO
  INTEGER,                     INTENT(OUT) :: MMF_DIM  
  INTEGER,      DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE_C), DIMENSION(NF), INTENT(IN)    :: FIELD_C
  TYPE(ATOM_C),                INTENT(IN)    :: ATOM__C                     

!  write(*,*) "C WRAPPER:", info,NM,NF,MODES_NUM
  IF (INFO.EQ.0) THEN
 !    write(*,*) "C wrapper"
     if(info.eq.0) CALL MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,COUPLING,INFO)
     H_FLOQUET_SIZE = SIZE(H_FLOQUET,1)
     MMF_DIM = H_FLOQUET_SIZE
!     write(*,*) "C wrapper"
     !write(*,*) h_floquet_size,size(H_floquet,1),coupling(1)%V(1,1)!, atom_%d_bare,h_floquet(1,1)
     !CALL WRITE_MATRIX(H_FLOQUET)
  ELSE
     WRITE(*,*) "THERE IS AN ERROR UPSTREAM :INFO.NE.0", INFO
  END IF
END SUBROUTINE MULTIMODEFLOQUETMATRIX_C_PYTHON
