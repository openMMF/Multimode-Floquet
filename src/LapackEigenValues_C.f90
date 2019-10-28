SUBROUTINE PACKINGBANDMATRIX_C(N,A,KD,AB,INFO)

! brute force packing of a banded matrix
  
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: INFO
  INTEGER, INTENT(IN)    :: N,KD
  COMPLEX*16, DIMENSION(N,N)    :: A
  COMPLEX*16, DIMENSION(KD+1,N) :: AB

  INTEGER i,j,N_,I_

  
  CALL PACKINGBANDMATRIX(N,A,KD,AB,INFO)
  
END SUBROUTINE PACKINGBANDMATRIX_C


SUBROUTINE LAPACK_FULLEIGENVALUESBAND_C(AB,Z,KD,N,W,INFO)
  !eigenvalues/vectors of banded matrix ab
!AB, inout, packed banded matrix
!Z, out,eigenvectors
!KD out, calcuated eigenvectors
!N, in,matrix dimension
!W, out, eigenvalues
!INFO,inout, error flag

  !H is COMPLEX*16 array, dimension (N, N)
  !  69 *>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
  !  70 *>          leading N-by-N upper triangular part of A contains the
  !  71 *>          upper triangular part of the matrix A.  If UPLO = 'L',
  !  72 *>          the leading N-by-N lower triangular part of A contains
  !  73 *>          the lower triangular part of the matrix A.
  !  74 *>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
  !  75 *>          orthonormal eigenvectors of the matrix A.
  !  76 *>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
  !  77 *>          or the upper triangle (if UPLO='U') of A, including the
  !  78 *>          diagonal, is destroyed.
  !
  ! The eigenvector H(:,r) corresponds to the eigenvalue W_SPACE(r)
  !
  IMPLICIT NONE
  INTEGER,                                INTENT(IN)    :: N,KD
  COMPLEX*16,       DIMENSION(KD+1,N), INTENT(INOUT)    :: AB
  COMPLEX*16,       DIMENSION(N,N),       INTENT(INOUT) :: Z
  DOUBLE PRECISION, DIMENSION(N),         INTENT(INOUT) :: W
  INTEGER,                                INTENT(OUT)   :: INFO

  CALL LAPACK_FULLEIGENVALUESBAND(AB,Z,KD,N,W,INFO)

END SUBROUTINE LAPACK_FULLEIGENVALUESBAND_C




SUBROUTINE LAPACK_FULLEIGENVALUES_C(U_F,N,W_SPACE_,INFO)
!eigenvalues/vectors of matrix ab
!H, inout, packed banded matrix
! , out,eigenvectors
!N, in,matrix dimension
!W_space, out, eigenvalues
!INFO,inout, error flag

  !H is COMPLEX*16 array, dimension (N, N)
  !  69 *>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
  !  70 *>          leading N-by-N upper triangular part of A contains the
  !  71 *>          upper triangular part of the matrix A.  If UPLO = 'L',
  !  72 *>          the leading N-by-N lower triangular part of A contains
  !  73 *>          the lower triangular part of the matrix A.
  !  74 *>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
  !  75 *>          orthonormal eigenvectors of the matrix A.
  !  76 *>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
  !  77 *>          or the upper triangle (if UPLO='U') of A, including the
  !  78 *>          diagonal, is destroyed.
  !
  ! The eigenvector H(:,r) corresponds to the eigenvalue W_SPACE(r)
  !
  USE ARRAYS
  IMPLICIT NONE
  INTEGER,                          INTENT(IN)    :: N
  COMPLEX*16,       DIMENSION(N,N), INTENT(INOUT) :: U_F
  DOUBLE PRECISION, DIMENSION(N),   INTENT(INOUT) :: W_SPACE_
  INTEGER,                          INTENT(OUT)   :: INFO

  !ALLOCATE(H_FLOQUET(N,N))
  !H_FLOQUET = 1.0
  !write(*,*) H_FLOQUET(1,1)!,SIZE(H_FLOQUET,1),SIZE(W_SPACE_,1)
  CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,N,W_SPACE_,INFO)
  
  U_F = H_FLOQUET
  DEALLOCATE(H_FLOQUET)
  
END SUBROUTINE LAPACK_FULLEIGENVALUES_C


SUBROUTINE LAPACK_SELECTEIGENVALUES_C(H,N,W_SPACE,L1,L2,Z,INFO)
!selected eigenvalues/vectors of hermitian matrix
!H, inout, packed banded matrix
! , out,eigenvectors
!N, in,matrix dimension
!W_space, out, eigenvalues
!L1 ordinal lowest eigenvalue
!L2 ordinal highest eigenvlaue
!Z : eigenvectors
!INFO,inout, error flag

  !USE FLOQUET
  IMPLICIT NONE
  INTEGER,                        INTENT(IN)    :: N,L1,L2
  COMPLEX*16, DIMENSION(:,:),     INTENT(INOUT) :: H
  COMPLEX*16, DIMENSION(:,:),     INTENT(OUT)   :: Z
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)   :: W_SPACE
  INTEGER,                        INTENT(OUT)   :: INFO

  CALL LAPACK_SELECTEIGENVALUES(H,N,W_SPACE,L1,L2,Z,INFO)

END SUBROUTINE LAPACK_SELECTEIGENVALUES_C

