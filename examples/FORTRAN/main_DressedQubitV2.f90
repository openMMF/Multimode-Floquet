PROGRAM MULTIMODEFLOQUET

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE ARRAYS 
  USE FLOQUETINITINTERFACE



  IMPLICIT NONE
  TYPE(ATOM)                                       ID
  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_MULTIFLOQUET
  INTEGER                                          INFO,m,INDEX0,r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,U_F1,U_F2

  DOUBLE PRECISION                              :: T1,T2
 

  ! ===================================================
  !PARAMETERS REQUIRED TO DEFINE THE DRESSED BASIS
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: U_FD           ! TRANSFORMATION FROM THE BARE TO THE DRESSED BASIS
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: E_DRESSED      ! DRESSED SPECTRUM
  INTEGER                                       :: DRESSINGFIELDS ! NUMBER OF DRESSING FIELDS
  INTEGER,          DIMENSION(:), ALLOCATABLE   :: DRESSINGFIELDS_INDICES ! IDENTITY OF THE DRESSING FIELDS




  ! ===================================================

 
  OPEN(UNIT=3,file="qubit_bareoscillation_V2.dat", action="write")
  OPEN(UNIT=4,file="qubit_dressedoscillation_V2.dat", action="write")


  INFO = 0
  CALL FLOQUETINIT(ID,'qubit',INFO)
  ALLOCATE(ENERGY(SIZE(J_Z,1)))
  ALLOCATE(H__(SIZE(J_z,1),SIZE(J_z,1)))
  ALLOCATE(U_AUX(SIZE(J_z,1),SIZE(J_z,1)))
  DEALLOCATE(ENERGY)
  




  
  ALLOCATE(MODES_NUM(3))

  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 1 !(DRIVING BY TWO HARMONICS)
  MODES_NUM(3) = 1 !(DRIVING BY A SECOND FREQUENCY)
  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(ID%D_BARE,ID%D_BARE))
  END DO
  
  FIELDS(1)%X         = 0.0
  FIELDS(1)%Y         = 0.0
  FIELDS(1)%Z         = 1.0
  FIELDS(1)%phi_x     = 0.0
  FIELDS(1)%phi_y     = 0.0
  FIELDS(1)%phi_z     = 0.0
  FIELDS(1)%omega     = 0.0
  FIELDS(1)%N_Floquet = 0

  FIELDS(2)%X         = 0.125/2.0
  FIELDS(2)%Y         = 0.0
  FIELDS(2)%Z         = 0.0
  FIELDS(2)%phi_x     = 0.0
  FIELDS(2)%phi_y     = 0.0
  FIELDS(2)%phi_z     = 0.0
  FIELDS(2)%omega     = 1.0
  FIELDS(2)%N_Floquet = 5
  
  FIELDS(3)%X         = 0.2*FIELDS(2)%X
  FIELDS(3)%Y         = 0.0
  FIELDS(3)%Z         = FIELDS(2)%X
  FIELDS(3)%phi_x     = 0.0
  FIELDS(3)%phi_y     = 0.0
  FIELDS(3)%phi_z     = 0.0
  FIELDS(3)%omega     = FIELDS(2)%X/2.0
  FIELDS(3)%N_Floquet = 5

  D_MULTIFLOQUET = ID%D_BARE
  DO r=1,TOTAL_FREQUENCIES
     D_MULTIFLOQUET = D_MULTIFLOQUET*(2*FIELDS(r)%N_Floquet+1)
  END DO

  !=================================================================================
  !==== DEFINITION OF THE DRESSING FIELDS AND DRESSED BASIS AND VARIABLES NEEDED TO DEFINE THE MICROMOTION OPERATOR
  !=================================================================================

  DRESSINGFIELDS = 2   ! NUMBER OF DRESSING FIELDS
  ALLOCATE(DRESSINGFIELDS_INDICES(DRESSINGFIELDS)) ! ARRAY THAT TELL US WHICH OF THE FIELD DEFINED ABOVE ARE THE DRESSING ONES
  DRESSINGFIELDS_INDICES(1) = 1 
  DRESSINGFIELDS_INDICES(2) = 2


























  



  CALL MICROMOTIONFOURIERDRESSEDBASIS(ID,DRESSINGFIELDS_INDICES,MODES_NUM,FIELDS,U_FD,E_DRESSED,INFO)
  write(*,*) e_dressed, size(E_dressed,1)
  
  !=================================================================================
  !== MULTIMODE FLOQUET DRESSED BASIS AND TIME-EVOLUTION OPERATOR IN THE BARE BASIS
  !=================================================================================
  CALL SETHAMILTONIANCOMPONENTS(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
  
  !==== ALLOCATE SPACE FOR THE MICROMOTION OPERATOR 
  ALLOCATE(U_F1(ID%D_BARE,ID%D_BARE))
  ALLOCATE(U_F2(ID%D_BARE,ID%D_BARE))
  

  DO r=1,64,4

!!$!========= FIND THE MULTIMODE FLOQUET SPECTRUM 
      

     FIELDS(3)%omega     = FIELDS(1)%Z - FIELDS(2)%X + 2.0*(r-1)*FIELDS(2)%X/64
     CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)          
     ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
     E_FLOQUET = 0.0  
     U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
     CALL LAPACK_FULLEIGENVALUES(U_F,SIZE(U_F,1),E_FLOQUET,INFO)
     DEALLOCATE(H_FLOQUET)

     
     ! ===== EVALUATE TIME-EVOLUTION OPERATOR 

     T1 = 0.0
     DO m=1,512,4
        T2 = (m-1)*16.0*100.0/128.0


        
        ! ===== EVALUATE TIME-EVOLUTION OPERATOR  IN THE BARE BASIS
        U_aux = 0.0
        CALL MULTIMODETIMEEVOLUTINOPERATOR(SIZE(U_F,1),SIZE(MODES_NUM,1),MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,FIELDS,T1,T2,U_AUX,INFO) 
        WRITE(3,*) FIELDS(3)%OMEGA,t2,ABS(U_AUX)**2
        
!!$     !=================================================================================
!!$     !== TRANSFORM THE TIME-EVOLUTION OPERATOR TO THE DRESSED BASIS
!!$     !=================================================================================
!!$        
!!$     !== BUILD THE TIME-DEPENDENT TRANSFORMATIONO BETWEEN THE BARE AND THE RF DRESSED BASIS       
        info =0         
        CALL MICROMOTIONDRESSEDBASIS(ID,MODES_NUM,DRESSINGFIELDS_INDICES,FIELDS,U_FD,&
             & E_DRESSED,T1,U_F1,INFO) 
        CALL MICROMOTIONDRESSEDBASIS(ID,MODES_NUM,DRESSINGFIELDS_INDICES,FIELDS,U_FD,&
             & E_DRESSED,T2,U_F2,INFO) 
        
        ! ---- CALCULATE THE TIME-EVOLUTION OPERATOR IN THE DRESSED BASIS USING THE PREVIOUSLY CALCULATED IN THE BARE BASIS
        U_AUX = MATMUL(TRANSPOSE(CONJG(U_F2)),MATMUL(U_AUX,U_F1)) 
        WRITE(4,*) FIELDS(3)%OMEGA,t2,ABS(U_AUX)**2

     END DO
!     write(*,*)
     DEALLOCATE(E_FLOQUET)
     WRITE(3,*)
     WRITE(4,*)
  END DO
  
END PROGRAM MULTIMODEFLOQUET

