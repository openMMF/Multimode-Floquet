!!$PROGRAM CRCPACKINGTEST
!!$
!!$  INTEGER                       :: N,DIM
!!$  INTEGER                       :: INFO
!!$  CHARACTER                     :: UPLO
!!$  COMPLEX*16,DIMENSION(5,5)     :: A
!!$
!!$  COMPLEX*16, DIMENSION(13) :: VALUES
!!$  INTEGER,    DIMENSION(13) :: COLUMNS
!!$  INTEGER,    DIMENSION(6) :: ROWINDEX
!!$
!!$  N    = SIZE(A,1)
!!$  DIM  =  SIZE(VALUES)
!!$  UPLO = "L"
!!$  A = 0
!!$  A(1,:) = (/ 1,-1, 0,-3, 0/)
!!$  A(2,:) = (/-2, 5, 0, 0, 0/)
!!$  A(3,:) = (/ 0, 0, 4, 6, 4/)
!!$  A(4,:) = (/-4, 0, 2, 0, 0/)
!!$  A(5,:) = (/ 0, 8, 0, 0,-5/)
!!$
!!$  VALUES   = 0 
!!$  COLUMNS  = 0
!!$  ROWINDEX = 0
!!$  INFO     = 0
!!$
!!$  CALL VARCRCPACKING(N,DIM,UPLO,A,VALUES,COLUMNS,ROWINDEX,INFO)
!!$  CALL WRITE_MATRIX(N,N, REAL(A))
!!$  !WRITE(*,*) VALUES
!!$  CALL WRITE_MATRIX(1,DIM,REAL(VALUES))
!!$  CALL WRITE_MATRIX_INT(1,DIM,COLUMNS)
!!$  CALL WRITE_MATRIX_INT(1,N+1,ROWINDEX)
!!$
!!$END PROGRAM CRCPACKINGTEST


!!$SUBROUTINE WRITE_MATRIX(N1,N2,A)
!!$  DOUBLE PRECISION, DIMENSION(N1,N2) :: A
!!$  CHARACTER(LEN=55) STRING
!!$  CHARACTER(LEN=55) aux_char
!!$  integer :: aux
!!$
!!$  aux = int(UBOUND(A,2))
!!$  write(aux_char,"(I4)") aux
!!$  aux_char = trim(aux_char)
!!$  write(string,"(A1,I4,A6)") "(",aux,"E15.6)"
!!$  !!string = trim(string)
!!$  !WRITE(*,*) string, aux,LBOUND(A,1),UBOUND(A,1),LBOUND(A,2),UBOUND(A,2)
!!$
!!$  DO I = LBOUND(A,1), UBOUND(A,1)
!!$     WRITE(*,string) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
!!$  END DO
!!$  WRITE(*,*)
!!$!  WRITE(*,*)
!!$  !101 format(string)
!!$END SUBROUTINE WRITE_MATRIX
!!$
!!$
!!$SUBROUTINE WRITE_MATRIX_INT(N1,N2,A)
!!$  INTEGER, DIMENSION(N1,N2) :: A
!!$  WRITE(*,*)
!!$  DO I = LBOUND(A,1), UBOUND(A,1)
!!$     WRITE(*,*) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
!!$  END DO
!!$END SUBROUTINE WRITE_MATRIX_INT
!!$

SUBROUTINE VARCRCPACKING(N,DIM,UPLO,zero,A,VALUES,COLUMNS,ROWINDEX,INFO)

  IMPLICIT NONE
  INTEGER,                   INTENT(IN)    :: N
  INTEGER,                   INTENT(INOUT) :: INFO,DIM
  CHARACTER,                 INTENT(IN)    :: UPLO
  DOUBLE PRECISION,          INTENT(IN)    :: ZERO
  COMPLEX*16,DIMENSION(N,N), INTENT(IN)    :: A

  COMPLEX*16, DIMENSION(DIM), INTENT(OUT) :: VALUES
  INTEGER,    DIMENSION(DIM), INTENT(OUT) :: COLUMNS
  INTEGER,    DIMENSION(N+1), INTENT(OUT) :: ROWINDEX


  INTEGER I,J, counter,aux


  VALUES   = 0
  counter  =  1
  columns  = -1
  ROWINDEX = -1

  SELECT CASE (UPLO)
     
  CASE('F')
     DO I=1,N
        DO J=1,N
           IF(ABS(A(I,J)).GT.ZERO .OR. I.EQ.J ) THEN
              VALUES(counter)  = A(I,J)
              COLUMNS(counter) = J
              IF(ROWINDEX(I).EQ.-1) ROWINDEX(I) = counter
              counter = counter + 1
           END IF
        END DO
        IF(I.EQ.N-1) aux = COUNTER   
     END DO
     AUX = COUNTER - 1 - AUX
     ROWINDEX(N+1) = aux + ROWINDEX(N) + 1
     

  CASE ('U')
     DO I=1,N
        DO J=I,N
           IF(ABS(A(I,J)).GT.ZERO .OR. I.EQ.J) THEN
              VALUES(counter)  = A(I,J)
              COLUMNS(counter) = J
              IF(ROWINDEX(I).EQ.-1) ROWINDEX(I) = counter
              counter = counter + 1
           END IF
        END DO
        IF(I.EQ.N-1) aux =COUNTER   
     END DO
     AUX = COUNTER - 1 - AUX

     ROWINDEX(N+1) = aux + ROWINDEX(N) + 1
     
  CASE('L')
     DO I=1,N
        DO J=1,I
           IF(ABS(A(I,J)).GT.ZERO  .OR. I.EQ.J) THEN
              VALUES(counter)  = A(I,J)
              COLUMNS(counter) = J
              IF(ROWINDEX(I).EQ.-1) ROWINDEX(I) = counter
              counter = counter + 1
           END IF
        END DO
        IF(I.EQ.N-1) aux = COUNTER   
     END DO
     AUX = COUNTER - 1 - AUX

     ROWINDEX(N+1) = aux + ROWINDEX(N) + 1
     
  END SELECT

  DIM = counter - 1

 
END SUBROUTINE VARCRCPACKING


