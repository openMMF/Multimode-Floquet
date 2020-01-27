SUBROUTINE MULTIMODEFLOQUETMATRIX_SP(ATOM__,NM,NF,MODES_NUM,FIELDS,VALUES_,ROW_INDEX_,COLUMN_,INFO)

  !ATOM_      (IN)    : type of quantum system
  !NM         (IN)    : number of modes
  !NF         (IN)    : number of driving fields
  !MODES_NUM  (IN)    : vector indicating the number of harmonics of each driving field (mode)
  !FIELDS     (IN)    : Fields
  !VALUES_    (OUT)   : Hamiltonian values
  !ROW_INDEX_ (OUT)   : if INFO = 6 vector indicating the row position of values
  !ROW_INDEX_ (OUT)   : IF INFO = 0 vector indicating the row_index position of values, as in the three array variation of the csr format
  !COLUMN_    (OUT)   : vector indicating the column position of the values
  !INFO       (INOUT) : error flag. INFO=0 means there is no error

  USE TYPES         !(modes.f90)
  USE MERGINGARRAYS !(utils.f90)
  USE QUICKSORTINTERFACE !(Modules.f90)

  IMPLICIT NONE
  INTEGER                  ,            INTENT(IN)    :: NM,NF
  TYPE(MODE), DIMENSION(NF),            INTENT(INOUT) :: FIELDS
  TYPE(ATOM),                           INTENT(IN)    :: ATOM__
  INTEGER,    DIMENSION(NM),            INTENT(IN)    :: MODES_NUM
  INTEGER,                              INTENT(INOUT) :: INFO
  COMPLEX*16, DIMENSION(:), ALLOCATABLE,INTENT(OUT)   :: VALUES_
  INTEGER,    DIMENSION(:), ALLOCATABLE,INTENT(OUT)   :: COLUMN_
  INTEGER,    DIMENSION(:), ALLOCATABLE,INTENT(OUT)   :: ROW_INDEX_


  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: V_AUX
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: VALUES,ARRAY_AUX,VALUES_OLD,VALUES_OLD_
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ROW,COLUMN,ARRAYI_AUX,ROW_OLD,COLUMN_OLD,ROW_INDEX,INDEX_ORDERROW,INDEX_ORDERROW_
  INTEGER,          DIMENSION(NM)               :: N_FLOQUET
  DOUBLE PRECISION, DIMENSION(NM)               :: OMEGA

  INTEGER     :: D_MULTIFLOQUET
  INTEGER     :: m,r,c,index,D,counter, N_MODES,values_dim,D_bare,t,FIELD_INDEX,MODE_INDEX
  CHARACTER*1 :: UPLO

  ! ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  ! ----- PARAMETERS REQUIRED TO TEST THE DIAGONALIZATION ROUTINE
  ! ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  !COMPLEX*16,      DIMENSION(:,:),ALLOCATABLE :: U_F
  !DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: E_FLOQUET
  !DOUBLE PRECISION :: E_L,E_R
  ! ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  ! ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

  N_MODES = NM
  D_bare  = ATOM__%D_BARE

  INFO       = 0
  N_FLOQUET  = 0
  counter    = 1
  !write(*,*) N_MODES,D_BARE,INFO,N_FLOQUET
  ALLOCATE(VALUES(D_bare*D_bare))
  ALLOCATE(ROW(D_bare*D_bare))
  ALLOCATE(COLUMN(D_bare*D_bare))

  !DEFINITION OF THE MATRICES: BY HAND  
  D_MULTIFLOQUET = SIZE(FIELDS(1)%V,1)
  DO r=2,NM
     FIELD_INDEX  = 2+SUM(MODES_NUM(2:r-1))
     N_FLOQUET(r) = FIELDS(FIELD_INDEX)%N_Floquet
     OMEGA(r)     = FIELDS(FIELD_INDEX)%OMEGA
     !WRITE(*,*) r,N_FLOQUET(r),field_index,modes_num(r),NM,NF
     IF(MODES_NUM(r).GT.N_FLOQUET(r)+1) THEN
        WRITE(*,*) "TO BUILD THE EXTENDED HAMILTONIAN THE NUMBER OF FLOQUET MODES + 1 MUST BE DEFINED"
        WRITE(*,*) "LARGER THAN THE NUMBER OF FIELD MODES"
        WRITE(*,*) "IN MODE ",r, "WE HAVE"
        WRITE(*,*) "FLOQUET MANIFOLDS =", N_FLOQUET(r)
        WRITE(*,*) "NUMBER OF MODES   =",MODES_NUM(r)
        INFO = -10
     END IF
     D_MULTIFLOQUET = D_MULTIFLOQUET*(2*N_FLOQUET(r)+1)
  END DO
  !write(*,*) D_MULTIFLOQUET
  ! IF(COUPLINGALLOCATED .EQV. .TRUE. ) THEN
  !    DO r=1,NF
  !       DEALLOCATE(FIELDS(r)%VALUES)
  !       DEALLOCATE(FIELDS(r)%COLUMN)
  !       DEALLOCATE(FIELDS(r)%ROW)
  !!       DEALLOCATE(FIELDS(r)%V)
  !    END DO
  !    !DEALLOCATE(ATOM_%E_BARE)
  ! END IF

  IF(INFO.GE.0) THEN
     ! COORDINATE PACKING OF EACH FIELD

     DO m=1,NF
        CALL COORDINATEPACKING(D_bare,FIELDS(m)%V,VALUES,ROW,COLUMN,index,INFO)
        ALLOCATE(FIELDS(m)%VALUES(index))
        IF(m.eq.1) FIELDS(m)%VALUES = VALUES(1:index)
        IF(m.ne.1) FIELDS(m)%VALUES = VALUES(1:index)/2.0
        !write(*,*) ABS(fields(m)%values)
!        write(*,*) m,fields(m)%omega, index
        ALLOCATE(FIELDS(m)%ROW(index))
        FIELDS(m)%ROW = ROW(1:index)
!        WRITE(*,*)m, FIELDS(m)%ROW
        ALLOCATE(FIELDS(m)%COLUMN(index))
        FIELDS(m)%COLUMN = COLUMN(1:index)
!        WRITE(*,*)m, FIELDS(m)%COLUMN
     END DO

     ! WE HAVE TO REDEFINE VALUES,ROW,COLUMN FOR m > 2 BECAUSE IN THE MULTIMODE FLOQUET
     ! HAMILTONIAN, THEY ARE REPEATED IN DIAGONAL MATRICES OF INCREASING DIMENSION:

     D = 1
     D_BARE = ATOM__%D_BARE
     index = 1
     FIELD_INDEX = 2+MODES_NUM(2)
     DO m=3,NM

        D = D*(2*N_FLOQUET(m-1)+1)                  
        FIELD_INDEX = 2+SUM(MODES_NUM(2:m-1))
        !write(*,*) field_index,m

        DO r=1,MODES_NUM(m)

           DEALLOCATE(VALUES)
           DEALLOCATE(ROW)
           DEALLOCATE(COLUMN)
           ALLOCATE(    VALUES(SIZE(FIELDS(FIELD_INDEX)%VALUES,1)))
           ALLOCATE(       ROW(SIZE(FIELDS(FIELD_INDEX)%VALUES,1)))
           ALLOCATE(    COLUMN(SIZE(FIELDS(FIELD_INDEX)%VALUES,1)))
           ALLOCATE(ARRAYI_AUX(SIZE(FIELDS(FIELD_INDEX)%VALUES,1)))
           VALUES = FIELDS(FIELD_INDEX)%VALUES
           ROW    = FIELDS(FIELD_INDEX)%ROW
           COLUMN = FIELDS(FIELD_INDEX)%COLUMN
!           write(*,*) size(column,1)
           !WRITE(*,*) m,r,row,real(values),D!SIZE(FIELDS(FIELD_INDEX)%VALUES,1)
           DO c=2,D ! it starts in 2, because 1 is explicitly done in the previous three lines
              CALL APPENDARRAYS(VALUES,FIELDS(FIELD_INDEX)%VALUES,INFO)
              ARRAYI_AUX = FIELDS(FIELD_INDEX)%ROW + (c-1)*D_BARE
              CALL APPENDARRAYSI(ROW,ARRAYI_AUX,INFO)
              ARRAYI_AUX =  FIELDS(FIELD_INDEX)%COLUMN + (c-1)*D_BARE
              CALL APPENDARRAYSI(COLUMN,ARRAYI_AUX,INFO)
              !WRITE(*,*) m,r,c,D,2*N_FLOQUET(m-1)+1,SIZE(VALUES,1)
              !           write(*,*) c,real(values)
              !           write(*,*) row
              !write(*,*) "C:",c,column
              !write(*,*)
           END DO
           !        WRITE(*,*) m,r, SIZE(VALUES,1),field_index
           !        write(*,*) REAL(VALUES)
           !        WRITE(*,*) ROW
           !        WRITE(*,*) COLUMN
           DEALLOCATE(FIELDS(FIELD_INDEX)%VALUES)
           DEALLOCATE(FIELDS(FIELD_INDEX)%ROW)
           DEALLOCATE(FIELDS(FIELD_INDEX)%COLUMN)
           DEALLOCATE(ARRAYI_AUX)
           ALLOCATE(FIELDS(FIELD_INDEX)%VALUES(SIZE(VALUES,1)))
           ALLOCATE(   FIELDS(FIELD_INDEX)%ROW(SIZE(VALUES,1)))
           ALLOCATE(FIELDS(FIELD_INDEX)%COLUMN(SIZE(VALUES,1)))
           FIELDS(FIELD_INDEX)%VALUES = VALUES
           FIELDS(FIELD_INDEX)%ROW    = ROW
           FIELDS(FIELD_INDEX)%COLUMN = COLUMN
           FIELD_INDEX = FIELD_INDEX + 1
        END DO
        !     DO r=1,SIZE(VALUES,1)
        !        WRITE(*,*) r,ROW(r),COLUMN(r),REAL(VALUES(r))
        !     END DO
        D_BARE = D_BARE*D
     END DO
     DEALLOCATE(VALUES)
     DEALLOCATE(ROW)
     DEALLOCATE(COLUMN)


     ! BUILDING THE COORDINATE PACKING OF THE MULTIMODE HAMILTONIAN MATRIX
     D_BARE = ATOM__%D_BARE
     D = D_bare
     values_dim = 0  
     ALLOCATE(VALUES_OLD(SIZE(FIELDS(1)%VALUES,1)))
     ALLOCATE(ROW_OLD(SIZE(FIELDS(1)%VALUES,1)))
     ALLOCATE(COLUMN_OLD(SIZE(FIELDS(1)%VALUES,1)))
     ALLOCATE(VALUES(SIZE(FIELDS(1)%VALUES)))
     ALLOCATE(ROW(SIZE(FIELDS(1)%VALUES)))
     ALLOCATE(COLUMN(SIZE(FIELDS(1)%VALUES)))
     
     VALUES     = FIELDS(1)%VALUES
     ROW        = FIELDS(1)%ROW
     COLUMN     = FIELDS(1)%COLUMN
     VALUES_OLD = FIELDS(1)%VALUES
     ROW_OLD    = FIELDS(1)%ROW
     COLUMN_OLD = FIELDS(1)%COLUMN  

     FIELD_INDEX = 2
     DO m=2,NM
        ALLOCATE(ARRAYI_AUX(SIZE(FIELDS(m)%ROW,1)))
        ALLOCATE(VALUES_OLD_(SIZE(VALUES_OLD,1)))
        VALUES_OLD_ = VALUES_OLD

        DO c = 1,2*N_FLOQUET(m) + 1 !THIS LOOP DEFINES THE DIAGONAL MATRIX
           DO t=1,SIZE(VALUES_OLD,1)
              IF(COLUMN_OLD(t).EQ.ROW_OLD(t)) VALUES_OLD_(t) = VALUES_OLD(t) + 1.0*(-N_FLOQUET(m)+c-1)*OMEGA(m)
           END DO
           IF(c.NE.1) THEN
              CALL APPENDARRAYS(VALUES,VALUES_OLD_,INFO)
              CALL APPENDARRAYSI(COLUMN,COLUMN_OLD+(c-1)*D,INFO)
              CALL APPENDARRAYSI(ROW,ROW_OLD+(c-1)*D,INFO)
           ELSE
              VALUES = VALUES_OLD_
           END IF
        END DO
        !     WRITE(*,*) REAL(VALUES)
        !     WRITE(*,*) ROW
        DO c=1,MODES_NUM(m) !NOW WE POPULATE AS MANY UPPER AND LOWER OFF-DIAGONAL MATRICES AS HARMONICS. 
           DO r = 1,2*N_FLOQUET(m) + 1 - c           
              CALL APPENDARRAYS(VALUES,FIELDS(FIELD_INDEX)%VALUES,INFO)           
              CALL APPENDARRAYS(VALUES,CONJG(FIELDS(FIELD_INDEX)%VALUES),INFO)           

              ARRAYI_AUX =  FIELDS(FIELD_INDEX)%COLUMN + (r-1+c)*D
              CALL APPENDARRAYSI(COLUMN,ARRAYI_AUX,INFO)           
              ARRAYI_AUX =  FIELDS(FIELD_INDEX)%ROW + (r-1)*D
              CALL APPENDARRAYSI(COLUMN,ARRAYI_AUX,INFO)           

              ARRAYI_AUX =  FIELDS(FIELD_INDEX)%ROW + (r-1)*D
              CALL APPENDARRAYSI(ROW,ARRAYI_AUX,INFO)
              ARRAYI_AUX =  FIELDS(FIELD_INDEX)%COLUMN + (r-1+c)*D
              CALL APPENDARRAYSI(ROW,ARRAYI_AUX,INFO)

           END DO
           FIELD_INDEX = FIELD_INDEX + 1
        END DO
        !WRITE(*,*) REAL(VALUES)
        !     WRITE(*,*) ROW

        D = D*(2*N_FLOQUET(m)+1)
        DEALLOCATE(VALUES_OLD_)
        DEALLOCATE(VALUES_OLD)
        DEALLOCATE(ROW_OLD)
        DEALLOCATE(COLUMN_OLD)
        DEALLOCATE(ARRAYI_AUX)

        ALLOCATE(VALUES_OLD(SIZE(VALUES,1)))
        ALLOCATE(ROW_OLD(SIZE(VALUES,1)))
        ALLOCATE(COLUMN_OLD(SIZE(VALUES,1)))
        VALUES_OLD = VALUES
        ROW_OLD    = ROW
        COLUMN_OLD = COLUMN
     END DO

     DEALLOCATE(VALUES_OLD)
     DEALLOCATE(ROW_OLD)
     DEALLOCATE(COLUMN_OLD)

     IF(info.eq.6) THEN ! values, row, column representation of the H matrix
        VALUES_    = VALUES
        ROW_INDEX_ = ROW
        COLUMN_    = COLUMN
     ELSE
     ! BUILDING THE VARIATION CRS PACKING OF THE MULTIMODE HAMILTONIAN MATRIX, USING THE COORDINATE PACKING

     !  DO r=1,size(VALUES,1)
     !     WRITE(*,*) r,ROW(r),COLUMN(r),REAL(VALUES(r))
     !  END DO
     !    WRITE(*,*) SIZE(VALUES,1)
!        WRITE(*,*) REAL(VALUES)
!     WRITE(*,*) ROW
     !WRITE(*,*) COLUMN
!     WRITE(*,*) 

        ALLOCATE(INDEX_ORDERROW(SIZE(ROW,1)))
        CALL QUICK_SORT_INTEGERS(ROW,INDEX_ORDERROW,SIZE(ROW,1))
        !write(*,*) ROW
        !write(*,*) INDEX_ORDERROW
        !write(*,*) ROW(INDEX_ORDERROW)

        ROW    = ROW(INDEX_ORDERROW)
        COLUMN = COLUMN(INDEX_ORDERROW)
        
        !     write(*,*) index_orderrow
        !write(*,*) real(values)
        !     write(*,*) row
        !WRITE(*,*)
        !WRITE(*,*)
        !     WRITE(*,*) size(values,1), size(column,1), size(row,1)
        !     WRITE(*,*) D_MULTIFLOQUET
        
        
        ALLOCATE(ROW_INDEX(D_MULTIFLOQUET+1))
        !ALLOCATE(ROW_INDEX(SIZE(VALUES,1)+1))
        ROW_INDEX = -1
        ROW_INDEX(D_MULTIFLOQUET+1) = SIZE(VALUES,1)+1
        !ROW_INDEX(SIZE(VALUES,1)+1) = SIZE(VALUES,1)+1
        
        counter = 1
        
        D = 1
        ROW_INDEX(1)=1
        DO r = 2,SIZE(ROW,1)
           IF(ROW(r).EQ.ROW(r-1)) THEN
              counter = counter +1    
              !           write(*,*) r,row(r), counter,D
           ELSE
              D = D + counter
              ROW_INDEX(ROW(r)) =D
              !if(row(r).eq.2) write(*,*)"m3",r,ROW_INDEX(ROW(r)),row(r),D
              counter = 1
           END IF
        END DO
        m = 1
        do r=2,size(row_index,1)
           counter = row_index(r)-row_index(r-1)-1

!           write(*,*) row_index(r)-row_index(r-1),m,m+counter
           ALLOCATE(INDEX_ORDERROW_(counter+1))           
!           write(*,*) INDEX_ORDERROW(m:m+counter)
           CALL QUICK_SORT_INTEGERS(INDEX_ORDERROW(m:m+counter),INDEX_ORDERROW_,counter+1)
           INDEX_ORDERROW(m:m+counter) = INDEX_ORDERROW((m-1)+INDEX_ORDERROW_)
           COLUMN(m:m+counter)         = COLUMN((m-1)+INDEX_ORDERROW_)
           m = m + counter + 1
           DEALLOCATE(INDEX_ORDERROW_)
        end do
        VALUES = VALUES(INDEX_ORDERROW)
 
        !write(*,*) INDEX_ORDERROW
        !write(*,*) COLUMN
        
        !  WRITE(*,*) D_MULTIFLOQUET,nf
        !E_L = -60.0
        !E_R =  60.0
        !ALLOCATE(E_FLOQUET(D_MULTIFLOQUET))
        !ALLOCATE(U_F(D_MULTIFLOQUET,D_MULTIFLOQUET))
        !CALL MKLSPARSE_FULLEIGENVALUES(D_MULTIFLOQUET,SIZE(VALUES,1),VALUES,ROW_INDEX,COLUMN,E_L,E_R,E_FLOQUET,U_F,INFO)
        !write(*,*) E_FLOQUET
        
        ALLOCATE(VALUES_(SIZE(VALUES,1)))
        ALLOCATE(ROW_INDEX_(SIZE(ROW_INDEX,1)))
        ALLOCATE(COLUMN_(SIZE(COLUMN,1)))
        
        VALUES_    = VALUES
        ROW_INDEX_ = ROW_INDEX
        COLUMN_    = COLUMN
        
        DO m=1,NF
           DEALLOCATE(FIELDS(m)%VALUES)
           DEALLOCATE(FIELDS(m)%ROW)
           DEALLOCATE(FIELDS(m)%COLUMN)
        END DO
        
        
        DEALLOCATE(VALUES)
        DEALLOCATE(ROW)
        DEALLOCATE(COLUMN)
        DEALLOCATE(ROW_INDEX)
     END IF
  END IF
     !  write(*,*) column_(3080:3200)
END SUBROUTINE MULTIMODEFLOQUETMATRIX_SP ! _SP  sparse packing

