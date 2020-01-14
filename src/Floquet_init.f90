SUBROUTINE SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,JTOTAL_,ID,INFO)
  ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
  ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
  ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
  !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES

  USE ATOMIC_PROPERTIES
  USE TYPES

  IMPLICIT NONE
  CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: ATOMICSPECIE
  CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: MANIFOLD  !
  !  INTEGER,          OPTIONAL, INTENT(IN) :: JTOTAL
  DOUBLE PRECISION,  OPTIONAL, INTENT(IN) :: JTOTAL_
  TYPE(ATOM),INTENT(OUT) :: ID
  INTEGER, INTENT(INOUT) :: INFO

  DOUBLE PRECISION  JTOTAL
!  write(*,*) atomicspecie,jtotal_
  JTOTAL = JTOTAL_/2.0
  !  IF(JTOTAL.NE.JTOTAL_) JTOTAL = JTOTAL_

  INFO = 4 !UNDEFINED ID
  IF(PRESENT(ATOMICSPECIE)) THEN
     SELECT CASE (ATOMICSPECIE)
     CASE("87Rb")
        mass_at = 87*amu
        I       = I_87Rb
        g_I     = g_I_87Rb
        J       = J_87Rb
        g_J     = g_J_87Rb
        A       = A_87Rb
        a_s     = a_s_87Rb
        alpha_E = alpha_E_87Rb
        Fup     = Fup_87Rb
        Fdown   = Fdown_87Rb
        INFO    = 1
        ID_name = ID_name_87Rb
!        write(*,*) "# 87Rb"
     CASE ("6Li")
        mass_at = 6*amu
        I       = I_6Li
        g_I     = g_I_6Li
        J       = J_6Li
        g_J     = g_J_6Li
        A       = A_6Li
        a_s     = a_s_6Li
        alpha_E = alpha_E_6Li
        Fup     = Fup_6Li
        Fdown   = Fdown_6Li     
        INFO    = 1  ! ITS AN ATOM
        ID_name = ID_name_6Li
!        write(*,*) "# 6Li"
     CASE("qubit")
        mass_at = amu
        I       = I_qubit
        g_I     = I_qubit
        J       = J_qubit
        g_J     = J_qubit
        A       = A_qubit
        a_s     = a_s_qubit
        alpha_E = alpha_E_qubit
        Fup     = Fup_qubit
        Fdown   = Fdown_qubit
        INFO    = 2 ! ITS A QUBIT
        ID_name = ID_name_qubit

        !write(*,*) "# qubit"
     CASE("spin")
        mass_at = amu
        I       = I_spin
        g_I     = I_spin
        J       = J_spin
        g_J     = J_spin
        A       = A_spin
        a_s     = a_s_spin
        alpha_E = alpha_E_spin
        Fup     = Fup_spin
        Fdown   = Fdown_spin
        INFO    = 3 ! ITS A SPIN
        ID_name = ID_name_spin

        write(*,*) "# spin"
     CASE("lattice")
        INFO    = 5 ! ITS A lattice
        ID_name = ID_name_lattice

     END SELECT
  END IF
  IF(INFO.EQ.1) THEN
     ! HERE WE DEFINE THE ATOMIC PARAMETES ENTERING THE HAMILTONIAN.
     ! H = A I \cdot J + sum_fields mu_B g_J B \cdot J + mu_B g_I B \cdot I , when using the two mainifolds of the ground state
     ! H = sum_fields mu_B g_F B \cdot F

     J = L+S
     F    = Fup
     gF_2 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
     F    = Fdown
     gF_1 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
     G_F  = (g_J-g_I)/16.0
     IF(PRESENT(MANIFOLD)) THEN
        IF(MANIFOLD.EQ."U") THEN

           Ftotal           = Fup
           total_states_lsi = 2*Fup + 1
           gF               = gF_2
           ID%id_system = 1

        ELSE  IF(MANIFOLD.EQ."L") THEN

           Ftotal           = Fdown
           total_states_lsi = 2*Fdown + 1
           gF               = gF_1        
           ID%id_system = 2

        ELSE IF(MANIFOLD.EQ."B") THEN
           Ftotal = Fup

           !----- Counting the number of states
           J = L+S
           Total_states_LSI = 0
           DO WHILE(J.GE.ABS(L-S))
              F = I+J
              DO WHILE(F.GE.ABS(J-I))
                 !        write(*,*) Total_states_LSI,F,I,J
                 Total_states_LSI = Total_states_LSI + 2*int(F) + 1
                 F = F - 1
              END DO
              J= J - 1 ! only two values allowed for J: J = L+1/2 and J=L-1/2
           END DO
           J = L+S  ! reseting value of J modifeid before
           IF(s.eq.0.5 .and. i.LT.0.5) total_states_lsi = 2
           ID%id_system = 3
        END IF
     END IF
  ELSE IF(INFO .EQ. 2) THEN
     ! HERE WE DEFINE THE PARAMETERS OF THE QUBIT HAMILTONIAN
     ! H = HBAR OMEGA_0 S_x + HBAR OMEGA_T S_x + HBAR OMEGA_RABI cos(omega_1 t + phi)
     total_states_lsi = 2
     J  = 0.0
     S  = 0.5
     F  = 0.5
     L  = 0.0
     gF = 1.0
     Ftotal = 0.5
     ID%id_system = 4
  ELSE IF(INFO .EQ. 3) THEN
     ! HERE WE DEFINE THE PARRAMETES OF THE SPIN HAMILTONIAN
     ! H = HBAR OMEGA_0 S_x + HBAR OMEGA_T S_x + HBAR OMEGA_RABI cos(omega_1 t + phi)
     total_states_lsi = 2*JTOTAL + 1
     J      = 0.0
     S      = JTOTAL
     F      = JTOTAL
     L      = 0.0
     gF     = 1.0
     Ftotal = JTOTAL
     ID%id_system = 5
  ELSE IF(INFO .EQ. 4) THEN
     ! HERE WE DEFINE THE PARAMETES OF THE LATTICE HAMILTONIAN
     ID%id_system = 6     
     ID_name      = "unnamed"
  ELSE IF(INFO.EQ.5) THEN
     ID%id_system    = 7
     ID_name         = "lattice"
     IF(PRESENT(MANIFOLD)) PERIODIC        =  manifold ! boundary conditionss
     total_states_lsi = JTOTAL_   ! number of sites
  END IF
  !WRITE(*,*) 'set_atomic_properties',JTOTAL,JTOTAL_,TOTAL_STATES_LSI
  ID%D_bare = total_states_lsi

END SUBROUTINE SET_ATOMIC_PARAMETERS

SUBROUTINE DEALLOCATEALL(ID)

  USE ARRAYS
  IMPLICIT NONE
  INTEGER, INTENT(IN)::ID
  !  write(*,*) ID
  IF(ALLOCATED(H_IJ)) DEALLOCATE(H_IJ)
  IF(ALLOCATED(U_ZEEMAN)) DEALLOCATE(U_ZEEMAN)
  !  DEALLOCATE(U_RF)
  IF(ALLOCATED(HAMILTONIAN)) DEALLOCATE(HAMILTONIAN)
  !  DEALLOCATE(H_AUX) 
  !  DEALLOCATE(H_RF)
  !  DEALLOCATE(H_RF_DAGGER)
  !  DEALLOCATE(H_MW)
  !  DEALLOCATE(H_ALPHA)
  !  DEALLOCATE(H_ALPHA_DAGGER)

  IF(ALLOCATED(j_x)) DEALLOCATE(j_x)
  IF(ALLOCATED(j_y)) DEALLOCATE(j_y)
  IF(ALLOCATED(j_z)) DEALLOCATE(j_z)
  IF(ALLOCATED(Identity)) DEALLOCATE(Identity)
  IF(ALLOCATED(H_FLOQUET)) DEALLOCATE(H_FLOQUET)
  IF(ALLOCATED(H_FLOQUET_COPY)) DEALLOCATE(H_FLOQUET_COPY)
  SELECT CASE(ID)
  CASE(1)

  CASE(2)

  CASE(3)
     !     DEALLOCATE(jz_dash)
     DEALLOCATE(I_x)
     DEALLOCATE(I_y)
     DEALLOCATE(I_z)
     !     DEALLOCATE(CLEBSH_GORDAN_JtoF)

     !    DEALLOCATE(g_F_matrix)
     !    DEALLOCATE(Fx)
     !    DEALLOCATE(Fy)
     !    DEALLOCATE(Fz)
     !    DEALLOCATE(Fz_DASH)
     !    DEALLOCATE(Hamiltonian_F)
     !    DEALLOCATE(Identity_F)
     !    DEALLOCATE(F_t)
     !    DEALLOCATE(H_w)
     !    DEALLOCATE(H_M)
     !    DEALLOCATE(H_J)
  CASE(4)

  CASE(5)
     IF(ALLOCATED(j_x)) DEALLOCATE(j_x)
     IF(ALLOCATED(j_y)) DEALLOCATE(j_y)
     IF(ALLOCATED(j_z)) DEALLOCATE(j_z)
  END SELECT


END SUBROUTINE DEALLOCATEALL

MODULE FLOQUETINITINTERFACE
  INTERFACE FLOQUETINIT
     MODULE PROCEDURE FLOQUETINIT_QUBIT, FLOQUETINIT_SPIN,FLOQUETINIT_ALKALI
  END INTERFACE FLOQUETINIT
contains
  SUBROUTINE FLOQUETINIT_QUBIT(ID,atomicspecie,INFO)
    ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
    ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
    ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
    !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES
    ! calculate the dimenson of the Hilbert space
    ! initialize all the matrices required for a full Floquet calcuations
    ! Calculate the nuclear, electron and total angular momentum operators
    
    USE physical_constants ! Standard Module with constants
    USE ATOMIC_PROPERTIES  ! gF, F , etc. factors for several species
  USE subinterface       ! To ubroutines for representation of I and J operators
  USE ARRAYS
  USE SUBINTERFACE_LAPACK
  USE TYPES
  IMPLICIT NONE
  
  TYPE(ATOM),        INTENT(INOUT) :: ID
  CHARACTER (LEN=*), INTENT(IN)    :: ATOMICSPECIE
  INTEGER,           INTENT(INOUT) :: INFO
  
  INTEGER  r,D_F2,P,r_,p_
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: Energy

!  write(*,*) atomicspecie,info,id%d_bare
    
  INFO = 4
  CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,'U',0.1D1,ID,INFO)
  !------ ALLOCATE NEEDED ARRAYS: Hamiltonian and Lapack
  ALLOCATE(Energy(TOTAL_STATES_LSI))
  ALLOCATE(H_IJ(Total_states_LSI,Total_states_LSI))
  ALLOCATE(U_ZEEMAN(Total_states_LSI,Total_states_LSI))
  ALLOCATE(HAMILTONIAN(Total_states_LSI,Total_states_LSI))
  U_zeeman = 0
  
  IF(INFO.EQ.2) THEN
     ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
     CALL F_REPRESENTATION(j_x,j_y,j_z,0.05D1)
  END IF
  
END SUBROUTINE FLOQUETINIT_QUBIT

SUBROUTINE FLOQUETINIT_SPIN(ID,atomicspecie,jtotal,info)
  ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
  ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
  ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
  !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES


  ! calculate the dimenson of the Hilbert space
  ! initialize all the matrices required for a full Floquet calcuations
  ! Calculate the nuclear, electron and total angular momentum operators

  USE physical_constants ! Standard Module with constants
  USE ATOMIC_PROPERTIES  ! gF, F , etc. factors for several species
  USE subinterface       ! To ubroutines for representation of I and J operators
  USE ARRAYS
  USE SUBINTERFACE_LAPACK
  USE TYPES
  IMPLICIT NONE

  TYPE(ATOM),        INTENT(INOUT) :: ID
  CHARACTER (LEN=*), INTENT(IN)    :: ATOMICSPECIE
  DOUBLE PRECISION,  intent(in)    :: jtotal
  INTEGER,           INTENT(INOUT) :: INFO
  
  INTEGER  r,D_F2,P,r_,p_
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: Energy

!  write(*,*) atomicspecie
  INFO = 4
  CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,'B',JTOTAL,ID,INFO)
  !------ ALLOCATE NEEDED ARRAYS: Hamiltonian and Lapack
  ALLOCATE(Energy(TOTAL_STATES_LSI))
  ALLOCATE(H_IJ(Total_states_LSI,Total_states_LSI))
  ALLOCATE(U_ZEEMAN(Total_states_LSI,Total_states_LSI))
  ALLOCATE(HAMILTONIAN(Total_states_LSI,Total_states_LSI))
  U_zeeman = 0
  IF(INFO.EQ.3) THEN
     ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
     CALL F_REPRESENTATION(j_x,j_y,j_z,1.0D0*Ftotal)     
  ELSE 
     ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
     j_x = 0.0
     j_y = 0.0
     j_z = 0.0
  END IF

END SUBROUTINE FLOQUETINIT_SPIN

SUBROUTINE FLOQUETINIT_ALKALI(ID,atomicspecie,manifold,info)
  ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
  ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
  ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
  !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES


  ! calculate the dimenson of the Hilbert space
  ! initialize all the matrices required for a full Floquet calcuations
  ! Calculate the nuclear, electron and total angular momentum operators

  USE physical_constants ! Standard Module with constants
  USE ATOMIC_PROPERTIES  ! gF, F , etc. factors for several species
  USE subinterface       ! To ubroutines for representation of I and J operators
  USE ARRAYS
  USE SUBINTERFACE_LAPACK
  USE TYPES
  IMPLICIT NONE

  TYPE(ATOM),        INTENT(INOUT) :: ID
  CHARACTER (LEN=*), INTENT(IN)    :: ATOMICSPECIE
  CHARACTER (LEN=1), INTENT(IN)    :: MANIFOLD  
  INTEGER,           INTENT(INOUT) :: INFO
  
  INTEGER  r,D_F2,P,r_,p_
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: Energy

!  write(*,*) atomicspecie, ' alkali ', manifold
  INFO = 4
  CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,0.2D1,ID,INFO)
  !------ ALLOCATE NEEDED ARRAYS: Hamiltonian and Lapack
!  write(*,*) total_states_lsi,info
  ALLOCATE(Energy(TOTAL_STATES_LSI))
  ALLOCATE(H_IJ(Total_states_LSI,Total_states_LSI))
  ALLOCATE(U_ZEEMAN(Total_states_LSI,Total_states_LSI))
  ALLOCATE(HAMILTONIAN(Total_states_LSI,Total_states_LSI))
  U_zeeman = 0
!  IF(INFO.EQ.1 .AND. PRESENT(MANIFOLD)) THEN
     IF(INFO.EQ.1 .AND. MANIFOLD.EQ.'B') THEN
        ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
        ALLOCATE(I_x(Total_states_LSI,Total_states_LSI))
        ALLOCATE(I_y(Total_states_LSI,Total_states_LSI))
        ALLOCATE(I_z(Total_states_LSI,Total_states_LSI))
        j_x           = 0.0
        j_y           = 0.0
        j_z           = 0.0
        I_x           = 0.0
        I_y           = 0.0
        I_z           = 0.0

        !---------------- Build the angular momentum operators in the JImImJ basis,
        !---------------- in decresing order of m=mI+mJ, with mJ=1/2,-1/2
        CALL I_and_J_representations(j_x,j_y,j_z,I_x,I_y,I_z,L,S,I)

    ELSE IF(INFO.EQ.1 .AND. MANIFOLD.NE.'B') THEN
        ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
        CALL F_REPRESENTATION(j_x,j_y,j_z,1.0D0*Ftotal)
     END IF
!  END IF

END SUBROUTINE FLOQUETINIT_ALKALI
END MODULE FLOQUETINITINTERFACE


SUBROUTINE SETHAMILTONIANCOMPONENTS(ID,NM,NF,MODES_NUM,FIELD,INFO)
  ! ID  tYPE OF ATOM
  ! MODES_NUM, VECTOR. THE SIZE OF THE VECTOR TELL US THE NUMBER OF FREQUENCIES, AND THE VALUE OF EACH COMPONENT INDICATES THE NUMBER OF HARMONICS OF EACH FREQUENCI
  ! FIELDS : IN AND OUTPUT THE MATRICES
  ! INFO

  USE ARRAYS
  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE_LAPACK ! write_matrix interface

  IMPLICIT NONE
  INTEGER,                   INTENT(IN)    :: NM,NF
  TYPE(ATOM),                INTENT(IN)    :: ID
  INTEGER,    DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE), DIMENSION(NF), INTENT(INOUT) :: FIELD
  INTEGER,                   INTENT(INOUT) :: INFO

  INTEGER m,TOTAL_FREQUENCIES,N_FLOQUET_,KD
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: mFs
  DOUBLE PRECISION, DIMENSION(id%d_bare) :: E_ZEEMAN
  DOUBLE PRECISION :: RESONANTrfFREQUENCY

  !  WRITE(*,*) "# Setting the Hamiltonian components for a ",ID_name
  !  WRITE(*,*) "# with ",nm, " modes and ",nf," fields"
  TOTAL_FREQUENCIES = NF
  !  write(*,*) "# total frequencies:", total_frequencies,ID%id_system
  SELECT CASE(ID%id_system)
  CASE(3) ! ATOM, BOTH HYPERFINE MANIFOLDS
     U_ZEEMAN = 0.0
     DO m=1,TOTAL_FREQUENCIES
        IF(m.EQ.1) THEN
           FIELD(m)%V = A*(MATMUL(I_x,J_x) - MATMUL(I_y,J_y) + MATMUL(I_z,J_z)) + &
                &       mu_B*(g_J*(FIELD(m)%X*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)*J_x  + &
                & DCMPLX(0.0,-1.0)*FIELD(m)%Y*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)*J_y  + &
                &                  FIELD(m)%Z*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)*J_z) + &
                &            (g_I*(FIELD(m)%X*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)*I_x  + &
                & DCMPLX(0.0,-1.0)*FIELD(m)%Y*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)*I_y  + &
                &                  FIELD(m)%Z*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)*I_z)))
           FIELD(m)%V = FIELD(m)%V/A
           U_ZEEMAN = FIELD(m)%V 
           CALL LAPACK_FULLEIGENVALUES(U_ZEEMAN,SIZE(FIELD(m)%V,1),E_ZEEMAN,INFO) ! FIND THE BASIS OF ZEEMAN STATES
           FIELD(m)%V = MATMUL(TRANSPOSE(CONJG(U_ZEEMAN)),MATMUL(FIELD(m)%V,U_ZEEMAN))
        END IF

        IF(m.GT.1) THEN
           FIELD(m)%V = mu_B*(g_J*(FIELD(m)%X*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)*J_x  + &
                & DCMPLX(0.0,-1.0)*FIELD(m)%Y*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)*J_y  + &
                &                  FIELD(m)%Z*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)*J_z) + &
                &            (g_I*(FIELD(m)%X*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)*I_x  + &
                & DCMPLX(0.0,-1.0)*FIELD(m)%Y*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)*I_y  + &
                &                  FIELD(m)%Z*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)*I_z)))
           FIELD(m)%V = FIELD(m)%V/A
           FIELD(m)%V = MATMUL(TRANSPOSE(CONJG(U_ZEEMAN)),MATMUL(FIELD(m)%V,U_ZEEMAN))!WRITE THE COUPLINGS IN THE ZEEMAN BASIS
        END IF
     END DO

     KD = SIZE(U_ZEEMAN,1)
     DO m=2,TOTAL_FREQUENCIES-1
        N_FLOQUET_  = MODES_NUM(m)
        KD = KD*(2*N_FLOQUET_+1)
     END DO
     KD = KD + SIZE(U_ZEEMAN,1) - 1

  CASE(7) ! ONE DIMENSIONAL TIGHT BINDING LATTICE

  CASE(8) ! TWO-DIMENSIONAL TIGHT BINDING LATTICE

  CASE DEFAULT ! alkali ATOM with ONE MANIFOLD, QUBIT, SPIN

     U_ZEEMAN = 0.0
     DO m=1,TOTAL_FREQUENCIES
        IF(m.EQ.1) THEN
           FIELD(m)%V =            (FIELD(m)%X*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)*J_x  + &
                &  DCMPLX(0.0,-1.0)*FIELD(m)%Y*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)*J_y  + &
                FIELD(m)%Z*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)*J_z) 
           U_ZEEMAN   = FIELD(m)%V 
           CALL LAPACK_FULLEIGENVALUES(U_ZEEMAN,SIZE(FIELD(m)%V,1),E_ZEEMAN,INFO)
        END IF

        IF(m.GT.1) THEN
           FIELD(m)%V =            (FIELD(m)%X*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)*J_x  + &
                &  DCMPLX(0.0,-1.0)*FIELD(m)%Y*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)*J_y  + &
                FIELD(m)%Z*EXP(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)*J_z) 
           FIELD(m)%V = MATMUL(TRANSPOSE(CONJG(U_ZEEMAN)),MATMUL(FIELD(m)%V,U_ZEEMAN))
           !           write(*,*) m, field(m)%X

        END IF
     END DO
     KD = SIZE(U_ZEEMAN,1)
     DO m=2,TOTAL_FREQUENCIES-1
        N_FLOQUET_  = MODES_NUM(m)
        KD = KD*(2*N_FLOQUET_+1)
     END DO
     KD = KD + SIZE(U_ZEEMAN,1) - 1

  END SELECT

END SUBROUTINE SETHAMILTONIANCOMPONENTS


SUBROUTINE FLOQUETINIT_OLD(ID,jtotal,manifold,atomicspecie)
  ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
  ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
  ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
  !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES


  ! calculate the dimenson of the Hilbert space
  ! initialize all the matrices required for a full Floquet calcuations
  ! Calculate the nuclear, electron and total angular momentum operators

  USE physical_constants ! Standard Module with constants
  USE ATOMIC_PROPERTIES  ! gF, F , etc. factors for several species
  USE subinterface       ! To ubroutines for representation of I and J operators
  USE ARRAYS
  USE SUBINTERFACE_LAPACK
  USE TYPES
  IMPLICIT NONE

  TYPE(ATOM),        INTENT(INOUT)       :: ID
  CHARACTER (LEN=*), INTENT(IN),optional :: ATOMICSPECIE
  !integer, intent(in),optional :: jtotal
    DOUBLE PRECISION,  INTENT(IN),optional :: JTOTAL
  CHARACTER (LEN=1), INTENT(IN),optional :: MANIFOLD  
  ! INTEGER,           INTENT(INOUT)       :: INFO
  
  INTEGER  r,D_F2,P,r_,p_,info
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: Energy

  INFO = 4
  IF( PRESENT(ATOMICSPECIE) .AND. PRESENT(MANIFOLD) .AND. PRESENT(JTOTAL)) THEN
     CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,JTOTAL,ID,INFO)
  ELSE IF( PRESENT(ATOMICSPECIE) .AND. PRESENT(JTOTAL)) THEN
     CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,'B',JTOTAL,ID,INFO)
  ELSE IF(PRESENT(ATOMICSPECIE) .AND. PRESENT(MANIFOLD)) THEN
     CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,0.2D1,ID,INFO)
  ELSE IF (PRESENT(ATOMICSPECIE)) THEN 
     CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,'U',0.1D1,ID,INFO)
  END IF
  !------ ALLOCATE NEEDED ARRAYS: Hamiltonian and Lapack
  ALLOCATE(Energy(TOTAL_STATES_LSI))
  ALLOCATE(H_IJ(Total_states_LSI,Total_states_LSI))
  ALLOCATE(U_ZEEMAN(Total_states_LSI,Total_states_LSI))
  !  ALLOCATE(U_RF(Total_states_LSI,Total_states_LSI))
  ALLOCATE(HAMILTONIAN(Total_states_LSI,Total_states_LSI))
  !  ALLOCATE(H_AUX(Total_states_LSI,Total_states_LSI)) 
  !  ALLOCATE(H_RF(Total_states_LSI,Total_states_LSI))
  !  ALLOCATE(H_RF_DAGGER(Total_states_LSI,Total_states_LSI))
  !  ALLOCATE(H_MW(Total_states_LSI,Total_states_LSI))
  !  ALLOCATE(H_ALPHA(Total_states_LSI,Total_states_LSI))
  !  ALLOCATE(H_ALPHA_DAGGER(Total_states_LSI,Total_states_LSI))
  !  write(*,*) "Floquet inint says:", total_states_lsi
  U_zeeman = 0
  !  write(*,*) real(u_zeeman(1,8))
  !  call write_matrix(abs(u_zeeman)) 
  IF(INFO.EQ.1 .AND. PRESENT(MANIFOLD)) THEN
     IF(INFO.EQ.1 .AND. MANIFOLD.EQ.'B') THEN
        ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
        !ALLOCATE(jz_dash(Total_states_LSI,Total_states_LSI))
        ALLOCATE(I_x(Total_states_LSI,Total_states_LSI))
        ALLOCATE(I_y(Total_states_LSI,Total_states_LSI))
        ALLOCATE(I_z(Total_states_LSI,Total_states_LSI))
        !ALLOCATE(CLEBSH_GORDAN_JtoF(Total_states_LSI,Total_states_LSI))
        j_x           = 0.0
        j_y           = 0.0
        j_z           = 0.0
        I_x           = 0.0
        I_y           = 0.0
        I_z           = 0.0
        !CLEBSH_GORDAN_JtoF = 0.0

        !---------------- Build the angular momentum operators in the JImImJ basis,
        !---------------- in decresing order of m=mI+mJ, with mJ=1/2,-1/2
        CALL I_and_J_representations(j_x,j_y,j_z,I_x,I_y,I_z,L,S,I)

!!$        !(BECAUSE SOME ROUTINES ARE WRITEN IN THE F BASIS, BUT THIS WILL BE ELLIMINATED)
!!$        ALLOCATE(g_F_matrix(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(Fx(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(Fy(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(Fz(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(Fz_DASH(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(Hamiltonian_F(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(Identity_F(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(F_t(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(H_w(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(H_M(Total_states_LSI,Total_states_LSI))
!!$        ALLOCATE(H_J(Total_states_LSI,Total_states_LSI))
!!$        Hamiltonian_F = 0.0
!!$        Fx            = 0.0
!!$        Fy            = 0.0
!!$        Fz            = 0.0
!!$        F_t           = 0
!!$        Identity_F = 0.0
!!$        DO r=1,SIZE(Identity_F,1)
!!$           Identity_F(r,r) = 1.0
!!$        END DO
!!$
!!$
!!$        !CALL WRITE_MATRIX(j_x)
!!$        !  CALL WRITE_MATRIX(j_y)
!!$        !  CALL WRITE_MATRIX(j_z)
!!$        !
!!$        !  CALL WRITE_MATRIX(I_x)
!!$        !  CALL WRITE_MATRIX(I_y)
!!$        !  CALL WRITE_MATRIX(I_z)
!!$        !  CALL IJtoF_Matrix(L,S,I,INT(Total_states_LSI),CLEBSH_GORDAN_JtoF)
!!$        !CALL F_representation(Fx,Fy,Fz,Ftotal)
!!$
!!$        Fx = J_x+I_x
!!$        Fy = J_y+I_y
!!$        Fz = J_z+I_z
!!$
!!$        !  CALL WRITE_MATRIX(J_Z+I_Z)
!!$
!!$        !IF(TOTAL_STATES_LSI.EQ.8) THEN
!!$        F_t = 0
!!$        DO r=1,2*Fdown+1
!!$           F_t(r,r) = 1
!!$        END DO
!!$        DO r=2*Fdown+1,total_states_lsi
!!$           F_t(r,r) = -1
!!$        END DO
!!$        !END IF
!!$
!!$
!!$        !IF(TOTAL_STATES_LSI.EQ.8) THEN
!!$        g_F_matrix = 0
!!$        DO r=1,2*Fdown+1
!!$           g_F_matrix(r,r) = gF_1
!!$        END DO
!!$        DO r=2*Fdown+1,total_states_lsi
!!$           g_F_matrix(r,r) = gF_2
!!$        END DO
!!$        !END IF
!!$
!!$        HAMILTONIAN_F = 0.00001*(g_J*J_z+g_I*I_z) + 2*(MATMUL(I_Z,J_Z) +MATMUL(I_x,J_x) - MATMUL(I_y,J_y) )
!!$
!!$        CALL LAPACK_FULLEIGENVALUES(HAMILTONIAN_F,TOTAL_STATES_LSI,Energy,INFO)
!!$
!!$        Fx=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(Fx,HAMILTONIAN_F))
!!$        Fy=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(Fy,HAMILTONIAN_F))
!!$        Fz=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(Fz,HAMILTONIAN_F))
!!$        H_IJ = MATMUL(I_X,J_X)+MATMUL(I_Z,J_Z)-MATMUL(I_Y,J_Y)
!!$        H_IJ=MATMUL(TRANSPOSE(CONJG(HAMILTONIAN_F)),MATMUL(H_IJ,HAMILTONIAN_F))
!!$
!!$        DO r=1,TOTAL_STATES_lsi
!!$           DO P=1,TOTAL_STATES_lsi
!!$              IF(ABS(HAMILTONIAN_F(r,p)).LT.1E-2) HAMILTONIAN_F(r,p) = 0.0
!!$              IF(ABS(H_IJ(r,p)).LT.1E-2) H_IJ(r,p) = 0.0
!!$              IF(ABS(Fx(r,p)).LT.1E-2) Fx(r,p) = 0.0
!!$              IF(ABS(Fy(r,p)).LT.1E-2) Fy(r,p) = 0.0
!!$              IF(ABS(Fz(r,p)).LT.1E-2) Fz(r,p) = 0.0
!!$           END DO
!!$        END DO
!!$        CLEBSH_GORDAN_JtoF = real(HAMILTONIAN_F)
!!$
!!$        ! CALL WRITE_MATRIX(FZ)
!!$
!!$        DO r=1,TOTAL_STATES_LSI
!!$           IF(r.le.2*Fdown+1) r_ = -(r-2)
!!$           IF(r.gt.2*Fdown+1) r_ =   r-6
!!$           DO p=1,TOTAL_STATES_LSI
!!$              IF(p.le.2*Fdown+1) p_ = -(p-2)
!!$              IF(p.gt.2*Fdown+1) p_ =   p-6
!!$              H_W(r,p) = r_ - p_        
!!$           END DO
!!$        END DO
!!$
!!$        Fz_dash = 0.0
!!$        DO r=1,TOTAL_STATES_LSI
!!$           ! DC DRESSED HAMILTONIAN IN A FRAME OF REFERENCE WHERE THE TOP FIELD IS STATIC
!!$           FZ_DASH(r,r) = int(g_F_matrix(r,r)/abs(g_F_matrix(r,r)))*Fz(r,r)
!!$        END DO
!!$        !  CALL WRITE_MATRIX(1.0D0*(Fz))
!!$        !  CALL WRITE_MATRIX(1.0D0*(Fz_dash))
!!$        !  WRITE(*,*) Fz(3,3),Fz_dash(3,3)
!!$
!!$        DO r=1,TOTAL_STATES_LSI
!!$           IF(r.le.2*Fdown+1) r_ = (r-Fdown)
!!$           IF(r.gt.2*Fdown+1) r_ =  r-(2*Fdown+1 + Fup + 1)
!!$           DO p=1,TOTAL_STATES_LSI
!!$              IF(p.le.2*Fdown+1) p_ = (p-Fdown)
!!$              IF(p.gt.2*Fdown+1) p_ =  p-(2*Fdown+1 + Fup + 1)
!!$              H_M(p,r) = p_ - r_        
!!$              !        write(*,*) p,r,Fz_dash(p,p) , Fz_dash(r,r),Fz_dash(p,p) - Fz_dash(r,r)
!!$              !        H_M(p,r) = Fz_dash(p,p) - Fz_dash(r,r)
!!$           END DO
!!$        END DO
!!$
!!$        !  CALL WRITE_MATRIX(1.0D0*H_M)
!!$
!!$        Jz_dash = 0.0
!!$
!!$        DO r=1,TOTAL_STATES_LSI  
!!$           Jz_dash(r,r) = INT(g_F_matrix(r,r)/abs(g_F_matrix(r,r)))
!!$        END DO
!!$
!!$        !  CALL WRITE_MATRIX(1.0D0*(Jz_DASH))
!!$        DO r=1,TOTAL_STATES_LSI
!!$           !     IF(r.le.3) r_ = (r-2)
!!$           !     IF(r.gt.3) r_ =  r-6
!!$           DO p=1,TOTAL_STATES_LSI
!!$              !        IF(p.le.3) p_ = (p-2)
!!$              !        IF(p.gt.3) p_ =  p-6
!!$              !        H_M(p,r) = p_ - r_        
!!$              H_J(p,r) = Jz_dash(p,p) - Jz_dash(r,r)  
!!$           END DO
!!$        END DO
!!$
!!$        !  CALL WRITE_MATRIX(1.0D0*H_J)
!!$        !  CALL write_MATRIX(MATMUL(TRANSPOSE(CLEBSH_GORDAN_JtoF),CLEBSH_GORDAN_JtoF))
!!$        !  CALL write_MATRIX(MATMUL(CLEBSH_GORDAN_JtoF,TRANSPOSE(CLEBSH_GORDAN_JtoF)))
!!$  CALL WRITE_MATRIX(CLEBSH_GORDAN_JtoF)
!!$  CALL WRITE_MATRIX(Fx)
!!$  CALL WRITE_MATRIX(Fy)
!!$  call WRITE_MATRIX(Fz)
!!$  CALL WRITE_MATRIX(real(H_IJ))
    ELSE IF(INFO.EQ.1 .AND. MANIFOLD.NE.'B') THEN
        ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
        ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
        CALL F_REPRESENTATION(j_x,j_y,j_z,1.0D0*Ftotal)
     END IF
  END IF

  IF(INFO.EQ.2) THEN
     ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
     CALL F_REPRESENTATION(j_x,j_y,j_z,0.05D1)
  ELSE IF(INFO.EQ.3) THEN
     ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
     CALL F_REPRESENTATION(j_x,j_y,j_z,1.0D0*Ftotal)     
  ELSE IF(INFO.EQ.4) THEN
  ELSE IF(INFO.EQ.5) THEN
     ALLOCATE(j_x(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_y(Total_states_LSI,Total_states_LSI))
     ALLOCATE(j_z(Total_states_LSI,Total_states_LSI))
  END IF

END SUBROUTINE FLOQUETINIT_OLD


