!export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"; 
SUBROUTINE CRS2Square(VD,RD,DC,D,VALUES,ROW_INDEX,COLUMN,H_F,INFO)
  IMPLICIT NONE
  INTEGER,                    INTENT(IN)    :: RD,DC,VD,D
  COMPLEX*16, DIMENSION(VD),  INTENT(IN)    :: VALUES
  COMPLEX*16, DIMENSION(D,D), INTENT(OUT)   :: H_F
  INTEGER,    DIMENSION(RD),  INTENT(IN)    :: ROW_INDEX
  INTEGER,    DIMENSION(DC),  INTENT(IN)    :: COLUMN
  INTEGER,                    INTENT(INOUT) :: INFO
 
 INTEGER i, C,R,j,l

! write(*,*) "reconstruction", VD,RD,DC,D,row_index  
 R = 0
 h_f = 0.0
 l = 1
 H_F = 0.0
 DO j=1,RD-1
    DO i= ROW_INDEX(j),row_index(j+1)-1
       C = COLUMN(l)
       R = j
       if(R.GT.D .OR. C.GT.D) write(*,*) "Here a problem:", RD,j,i,l,R,C,column(l),ROW_INDEX(j+1),row_index(j),&
            & ROW_INDEX(j+1)-row_index(j)
       H_F(R,C) = values(l)
       if(c.gt.d) write(*,*) "warning: col"
       if(r.gt.d) write(*,*) "warning: row"
       if(l.gt.size(values,1)) write(*,*) "warning: val"
       l = l+1
    end do
 END DO
END SUBROUTINE CRS2Square

PROGRAM MULTIMODEFLOQUET

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SPARSE_INTERFACE
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE FLOQUETINITINTERFACE
  USE ARRAYS 

  IMPLICIT NONE
  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES
  INTEGER                                          INFO,m,INDEX0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,U_F1,U_F2,U_F1_red,U_F2_red,h_f
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  DOUBLE PRECISION                              :: T1,T2

  complex*16, dimension(:), allocatable:: val
  integer, dimension(:), allocatable :: row,col



  !PARAMETERS NEEDED TO DEFINE THE SPARSE MATRIX
  INTEGER,    DIMENSION(:), ALLOCATABLE :: ROW_INDEX,COLUMN
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VALUES
  !PARAMETERS TO DIAGONALIZE THE SPARSE MATRIX WITH MKL
  DOUBLE PRECISION :: E_L,E_R
  INTEGER          :: D_MULTIFLOQUET,r


!  OPEN(UNIT=3,FILE="qubit_oscillation_SP.dat",ACTION="WRITE")
!  OPEN(UNIT=4,FILE="qubit_avgtransition_SP.dat",ACTION="WRITE")


  INFO = 0
  CALL FLOQUETINIT(ID,'qubit',INFO)
  ALLOCATE(ENERGY(SIZE(J_Z,1)))
  ALLOCATE(H__(SIZE(J_z,1),SIZE(J_z,1)))
  ALLOCATE(P_AVG(SIZE(J_z,1),SIZE(J_z,1)))
  ALLOCATE(U_AUX(SIZE(J_z,1),SIZE(J_z,1)))
  H__ = J_z
  CALL LAPACK_FULLEIGENVALUES(H__,TOTAL_STATES_LSI,Energy,INFO)
  DEALLOCATE(H__)
  DEALLOCATE(ENERGY)
  
  
  ALLOCATE(MODES_NUM(3))

  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 1 !(DRIVING BY TWO HARMONICS)  
  MODES_NUM(3) = 1 !(DRIVING BY TWO HARMONICS)

  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(ID%D_BARE,ID%D_BARE))
  END DO
  
  FIELDS(1)%X    = 0.0
  FIELDS(1)%Y    = 0.0
  FIELDS(1)%Z    = 1.0
  FIELDS(1)%phi_x = 0.0
  FIELDS(1)%phi_y = 0.0
  FIELDS(1)%phi_z = 0.0
  FIELDS(1)%omega = 0.0
  FIELDS(1)%N_Floquet = 0

  FIELDS(2)%X     = 2.0
  FIELDS(2)%Y     = 0.0
  FIELDS(2)%Z     = 0.0
  FIELDS(2)%phi_x = 0.0
  FIELDS(2)%phi_y = 0.0
  FIELDS(2)%phi_z = 0.0
  FIELDS(2)%omega = 1.0
  FIELDS(2)%N_Floquet = 9

  FIELDS(3)%X     = 2.0
  FIELDS(3)%Y     = 0.0
  FIELDS(3)%Z     = 0.0
  FIELDS(3)%phi_x = 0.0
  FIELDS(3)%phi_y = 0.0
  FIELDS(3)%phi_z = 0.0
  FIELDS(3)%omega = 1.0
  FIELDS(3)%N_Floquet = 8

  D_MULTIFLOQUET = ID%D_BARE
  DO r=1,TOTAL_FREQUENCIES
     D_MULTIFLOQUET = D_MULTIFLOQUET*(2*FIELDS(r)%N_Floquet+1)
  END DO
!!$  allocate(h_f(5,5))
!!$  allocate(val(13))
!!$  allocate(row(6))
!!$  allocate(col(13))
!!$  val = (/1,-1,-3,-2,5,4,6,4,-4,2,7,8,-5/)
!!$  row = (/1,4,6,9,12,14/)
!!$  col = (/1,2,4,1,2,3,4,5,1,3,4,2,5/)
!!$  CALL CRS2SQUARE(13,6,13,5,val,row,col,h_f,info)
!!$  call write_matrix(real(h_f))
!!$  deallocate(h_f)
  DO m=1,1!128,4

     ! --- SET DRIVING PARAMETERS 
     FIELDS(2)%omega = 0.2 + (m-1)*2.0/128
     CALL SETHAMILTONIANCOMPONENTS(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)

     !--- FIND THE MULTIMODE FLOQUET SPECTRUM 
     CALL MULTIMODEFLOQUETMATRIX_SP(ID,SIZE(MODES_NUM,1),total_frequencies,MODES_NUM,FIELDS,VALUES,ROW_INDEX,COLUMN,INFO)
     ALLOCATE(H_F(D_MULTIFLOQUET,D_MULTIFLOQUET))
     !write(*,*) row_index(1:10)
     !CALL CRS2Square(size(VALUES,1),size(ROW_INDEX,1),size(COLUMN,1),D_MULTIFLOQUET,VALUES,ROW_INDEX,COLUMN,H_F,INFO)
     E_L = -50.0
     E_R =  50.0
     ALLOCATE(E_FLOQUET(D_MULTIFLOQUET))
     ALLOCATE(U_F(D_MULTIFLOQUET,D_MULTIFLOQUET))
     E_FLOQUET = 0.0
     U_F = 0.0
     CALL MKLSPARSE_FULLEIGENVALUES(D_MULTIFLOQUET,SIZE(VALUES,1),VALUES,ROW_INDEX,COLUMN,E_L,E_R,E_FLOQUET,U_F,INFO)

     !--- EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
     P_AVG = 0.0
     CALL MULTIMODETRANSITIONAVG(SIZE(U_F,1),size(MODES_NUM,1),FIELDS,MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,P_AVG,INFO) 
 !    WRITE(4,*) FIELDS(2)%omega,P_AVG
         
     !--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
     T1 = 0.0
     DO r=1,1!28
        T2 = r*32.0*4.0*atan(1.0)/128
        CALL MULTIMODETIMEEVOLUTINOPERATOR(SIZE(U_F,1),SIZE(MODES_NUM,1),MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,FIELDS,T1,T2,U_AUX,INFO) 
        P_AVG = ABS(U_AUX)**2
!        WRITE(*,*) t2,FIELDS(2)%OMEGA,ABS(U_AUX)**2
     END DO
!     WRITE(*,*)
     DEALLOCATE(E_FLOQUET)
     DEALLOCATE(U_F)

  END DO     
END PROGRAM MULTIMODEFLOQUET

