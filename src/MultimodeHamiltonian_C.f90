SUBROUTINE MULTIMODEFLOQUETMATRIX_C(ATOM__C,NM,NF,MODES_NUM,FIELD_C,H_FLOQUET_,INFO)
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
  !COMPLEX*16, DIMENSION(:,:),allocatable,intent(out) :: H_FLOQUET_
  COMPLEX*16, DIMENSION(2,2),intent(out) :: H_FLOQUET_

!  write(*,*) "C WRAPPER:", info,NM,NF,MODES_NUM
  IF (INFO.EQ.0) THEN
 !    write(*,*) "C wrapper"
     if(info.eq.0) CALL MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,COUPLING,INFO)
     H_FLOQUET_SIZE = SIZE(H_FLOQUET,1)
    ! ALLOCATE(H_FLOQUET_(H_FLOQUET_SIZE,H_FLOQUET_SIZE))
     H_FLOQUET_ = H_FLOQUET(1:2,1:2)
    ! H_FLOQUET_C    => H_FLOQUET
    ! write(*,*) "Floquet pointer: ", H_FLOQUET_C(1,1)
!     write(*,*) "C wrapper"
     !write(*,*) h_floquet_size,size(H_floquet,1),coupling(1)%V(1,1)!, atom_%d_bare,h_floquet(1,1)
     !CALL WRITE_MATRIX(H_FLOQUET)
  ELSE
     WRITE(*,*) "THERE IS AN ERROR UPSTREAM :INFO.NE.0", INFO
  END IF
END SUBROUTINE MULTIMODEFLOQUETMATRIX_C

FUNCTION MULTIMODEFLOQUETMATRIX_C_PT(ATOM__C,NM,NF,MODES_NUM,FIELD_C,INFO) RESULT(H_FLOQUET_)
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

  COMPLEX*16, DIMENSION(2,2) :: H_FLOQUET_
!  write(*,*) "C WRAPPER:", info,NM,NF,MODES_NUM
  IF (INFO.EQ.0) THEN
 !    write(*,*) "C wrapper"
     if(info.eq.0) CALL MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,COUPLING,INFO)
     H_FLOQUET_SIZE = SIZE(H_FLOQUET,1)
     H_FLOQUET_ = H_FLOQUET(1:2,1:2)
     !H_FLOQUET_C    => H_FLOQUET
     write(*,*) "Floquet pointer: ", H_FLOQUET_(1,1)
!     write(*,*) "C wrapper"
     !write(*,*) h_floquet_size,size(H_floquet,1),coupling(1)%V(1,1)!, atom_%d_bare,h_floquet(1,1)
     !CALL WRITE_MATRIX(H_FLOQUET)
  ELSE
     WRITE(*,*) "THERE IS AN ERROR UPSTREAM :INFO.NE.0", INFO
  END IF
END FUNCTION MULTIMODEFLOQUETMATRIX_C_PT


FUNCTION MULTIMODEFLOQUETMATRIX_C_PYTHON(ATOM__C,NM,NF,MODES_NUM,FIELD_C,INFO) RESULT(MMF_DIM)
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
  !INTEGER,                     INTENT(OUT) :: MMF_DIM  
  INTEGER,      DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE_C), DIMENSION(NF), INTENT(IN)    :: FIELD_C
  TYPE(ATOM_C),                INTENT(IN)    :: ATOM__C                     

  integer mmf_dim

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
     MMF_DIM = 0
  END IF
END FUNCTION MULTIMODEFLOQUETMATRIX_C_PYTHON
