SUBROUTINE MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,FIELD,INFO)
  !ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO
  !  USE FLOQUET
  !ATOM_ type atom, -> dimension of the bare Hilbert space
  !NM -> number of modes
  !NF -> Number of Fields
  !MODES_NUM -> number of harmonics of each mode
  !FIELD -> Field couplings
  !INFO


  USE ARRAYS
  USE ATOMIC_PROPERTIES
  USE TYPES
  !USE SUBINTERFACE_LAPACK

  IMPLICIT NONE
  INTEGER,                  INTENT(IN)    :: NM,NF
  INTEGER,                  INTENT(INOUT) :: INFO
  INTEGER,   DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE),DIMENSION(NF), INTENT(IN)    :: FIELD
  TYPE(ATOM),               INTENT(IN)    :: ATOM_                       

  INTEGER m,n,D,r,o,D_OLD,index_aux
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H_TEMP,H_STATIC,COUPLING,Z_M_COPY,H_FLOQUET_ROW
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: E_DRESSED
  INTEGER index_l_lower,index_l_upper,index_r_lower,index_r_upper,FIELD_INDEX

  INTEGER, DIMENSION(NM) :: N_FLOQUET

!  write(*,*) "# Creating the multimode hamiltonian for a", ID_name
!  write(*,*) "# with ", NM," modes and ", NF, "frequencies"
  
  INFO      = 0
  N_FLOQUET = 0
  !write(*,*) "# FORTRAN MULTIMODEFLOQUETMATRIX SAYS:",nm 
  DO n=2,NM
     FIELD_INDEX = 2+SUM(MODES_NUM(2:n-1))
     N_FLOQUET(n)=FIELD(FIELD_INDEX)%N_Floquet
 !    write(*,*) n,N_FLOQUET(n),field_index,modes_num(n),NM,NF,info,n
     IF(modes_num(n).GT.N_FLOQUET(n)+1) THEN
        WRITE(*,*) "# TO BUILD THE EXTENDED HAMILTONIAN THE NUMBER OF FLOQUET MODES MUST BE DEFINED"
        WRITE(*,*) "# LARGER THAN THE NUMBER OF FIELD MODES"
        INFO = -10
     END IF
  END DO
!  write(*,*) "# Floquet modes:", N_FLOQUET

  IF(INFO.EQ.0) THEN
     D_OLD    = 1
     D        = ATOM_%D_BARE
     !write(*,*) D,info
     ALLOCATE(H_FLOQUET_COPY(D,D))
     H_FLOQUET_COPY = 0.0     
     H_FLOQUET_COPY = FIELD(1)%V  ! STATIC HAMILTONIAN

     FIELD_INDEX = 2 !?

     DO n=2,NM  ! RUN OVER EACH FREQUENCY


        ! D : UPDATED AT THE ENDO OF THE LOOP. DIMENSION OF THE MULTIMODE FLOQUET MATRIX
        ALLOCATE(H_STATIC(SIZE(H_FLOQUET_COPY,1),SIZE(H_FLOQUET_COPY,1)))
        H_STATIC  = H_FLOQUET_COPY  
        DEALLOCATE(H_FLOQUET_COPY)

        ALLOCATE(IDENTITY(D,D))
        IDENTITY  = 0.0
        DO m= 1,D 
           IDENTITY(m,m) = 1.0
        END DO

        ALLOCATE(COUPLING(D,D))
        ALLOCATE(H_FLOQUET_ROW(D,D*(2*N_FLOQUET(n)))) ! a row is made of blocks of COUPLING, defined for each harmonic
        H_FLOQUET_ROW = 0.0
        COUPLING  = 0.0

        !        IF(NF.NE.NM) THEN
        FIELD_INDEX =2+SUM(MODES_NUM(2:n-1))
        DO o=1,MODES_NUM(n) ! loop to define the blocks of  H_FLOQUET_ROW

           DO r=1,(2*N_Floquet(n-1)+1)!,(2*N_Floquet(n-1)+1)

              index_l_lower = ATOM_%D_BARE*(r - 1) + 1
              index_l_upper = ATOM_%D_BARE*(r - 1) + ATOM_%D_BARE
              index_r_lower = index_l_lower
              index_r_upper = index_l_upper
              COUPLING(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
                   &     FIELD(FIELD_INDEX)%V  ! COUPLING MATRIX OF MODE n
!              write(*,*) r,field_index
           END DO
           H_FLOQUET_ROW(1:D,(o-1)*D+1:o*D) = COUPLING/2.0
           FIELD_INDEX = FIELD_INDEX + 1
           !CALL WRITE_MATRIX(ABS(H_FLOQUET_ROW))
        END DO
        FIELD_INDEX =2+SUM(MODES_NUM(2:n-1))

!        write(*,*) n,N_FLOQUET(n),field_index,field(field_index)%omega
        !        write(*,*) n,FIELD_INDEX,N_Floquet(n),FIELD(FIELD_INDEX)%OMEGA
        D_OLD = D
        D     = D*(2*N_FLOQUET(n)+1)
        ALLOCATE(H_FLOQUET(D,D))
        H_FLOQUET = 0.0
        DO m=-N_floquet(n),N_Floquet(n)

           index_r_lower = (m+N_Floquet(n))*D_OLD + 1
           index_r_upper = index_r_lower + D_OLD - 1

           index_aux = 2*N_floquet(n)*D_OLD - D_OLD*(m+N_Floquet(n))
           index_l_lower = index_r_lower + D_OLD 
           index_l_upper = index_l_lower + index_aux - 1
           !     WRITE(*,*) m,index_r_lower,index_r_upper,index_L_lower,index_l_upper,D_OLD,index_aux,size(H_FLOQUET_ROW,2)
           IF(m.LT.N_Floquet(n)) THEN
              H_FLOQUET(index_r_lower:index_r_upper,index_l_lower:index_l_upper) = H_FLOQUET_ROW(1:D_OLD,1:index_aux)
              H_FLOQUET(index_l_lower:index_l_upper,index_r_lower:index_r_upper) = &
                   & TRANSPOSE(CONJG(H_FLOQUET_ROW(1:D_OLD,1:index_aux)))
           END IF
           H_FLOQUET(index_r_lower:index_r_upper,index_r_lower:index_r_upper) = H_STATIC + m*FIELD(FIELD_INDEX)%OMEGA*identity
 !          write(*,*) m,n,m*FIELD(FIELD_INDEX)%OMEGA,field_index,n_floquet(n)
        END DO
        !END IF
!!$        ELSE IF(NF.EQ.NM) THEN
!!$
!!$           DO r=1,(2*N_Floquet(n-1)+1)
!!$        
!!$              index_l_lower = ATOM_%D_BARE*(r - 1) + 1
!!$              index_l_upper = ATOM_%D_BARE*(r - 1) + ATOM_%D_BARE
!!$              index_r_lower = index_l_lower
!!$              index_r_upper = index_l_upper
!!$!              COUPLING(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
!!$!                   &     FIELD(n)%V  ! COUPLING MATRIX OF MODE n
!!$              !       write(*,*) n,r,index_l_lower,index_l_upper,BASIS
!!$           END DO
!!$  
!!$           DO m=-N_floquet(n),N_Floquet(n)
!!$              
!!$              index_l_lower = (m + N_Floquet(n)    )*SIZE(COUPLING,1) + 1
!!$              index_l_upper = index_l_lower + SIZE(COUPLING,1) - 1
!!$              index_r_lower =  index_l_lower
!!$              index_r_upper =  index_l_upper
!!$              H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
!!$                   &  1.0*H_STATIC + 1.0*m*FIELD(n)%OMEGA*IDENTITY
!!$              !        call write_matrix(ABS(H_STATIC))
!!$              !call write_matrix(ABS(H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper)))   
!!$              !        write(*,*) n,m,hbar,FIELD(n)%OMEGA,A
!!$              IF(m.LT.N_floquet(n)) THEN
!!$                 
!!$                 index_l_lower =  (m + N_Floquet(n) + 1)*SIZE(COUPLING,1) + 1
!!$                 index_l_upper =  index_l_lower + SIZE(COUPLING,1) - 1
!!$                 index_r_lower =  (m + N_Floquet(n)    )*SIZE(COUPLING,1) + 1
!!$                 index_r_upper =  index_r_lower + SIZE(COUPLING,1) - 1
!!$                 H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
!!$                      &     0.5*COUPLING
!!$                 !           write(*,*) index_l_lower,index_l_upper,index_r_lower,index_r_upper
!!$                 
!!$                 index_l_lower =  (m + N_Floquet(n)    )*SIZE(COUPLING,1) + 1
!!$                 index_l_upper =  index_r_lower + SIZE(COUPLING,1)  - 1        
!!$                 index_r_lower =  (m + N_Floquet(n) + 1)*SIZE(COUPLING,1) + 1
!!$                 index_r_upper =  index_l_lower + SIZE(COUPLING,1)  - 1         
!!$                 H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
!!$                      &     0.5*TRANSPOSE(CONJG(COUPLING))
!!$              END IF
!!$              
!!$           END DO
        !        END IF

        DEALLOCATE(IDENTITY)
        DEALLOCATE(COUPLING)
        DEALLOCATE(H_FLOQUET_ROW)

        IF(n.LT.NM) THEN
           ALLOCATE(H_FLOQUET_COPY(D,D))
           H_FLOQUET_COPY = H_FLOQUET
           DEALLOCATE(H_FLOQUET)
           DEALLOCATE(H_STATIC)
        END IF

     END DO
!     write(*,*) "# Dimension of the matrix:", size(H_FLOQUET,1)

!     CALL WRITE_MATRIX(real(H_FLOQUET))
  ELSE

  END IF
END SUBROUTINE MULTIMODEFLOQUETMATRIX

