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
        !write(*,*) "# 6Li"
     CASE ("85Rb")
        mass_at = mass_at_85Rb
        I       = I_85Rb
        g_I     = g_I_85Rb
        J       = J_85Rb
        g_J     = g_J_85Rb
        A       = A_85Rb
        a_s     = a_s_85Rb
        alpha_E = alpha_E_85Rb
        Fup     = Fup_85Rb
        Fdown   = Fdown_85Rb     
        INFO    = 1  ! ITS AN ATOM
        ID_name = ID_name_85Rb
        !write(*,*) "# 6Li"
     CASE ("123Cs")
        mass_at = mass_at_123Cs
        I       = I_123Cs
        g_I     = g_I_123Cs
        J       = J_123Cs
        g_J     = g_J_123Cs
        A       = A_123Cs
        a_s     = a_s_123Cs
        alpha_E = alpha_E_123Cs
        Fup     = Fup_123Cs
        Fdown   = Fdown_123Cs     
        INFO    = 1  ! ITS AN ATOM
        ID_name = ID_name_123Cs
        !write(*,*) "# 6Li"
     CASE ("41K")
        mass_at = mass_at_41K
        I       = I_41K
        g_I     = g_I_41K
        J       = J_41K
        g_J     = g_J_41K
        A       = A_41K
        a_s     = a_s_41K
        alpha_E = alpha_E_41K
        Fup     = Fup_41K
        Fdown   = Fdown_41K     
        INFO    = 1  ! ITS AN ATOM
        ID_name = ID_name_41K
        !write(*,*) "# 6Li"
     CASE ("23Na")
        mass_at = mass_at_23Na
        I       = I_23Na
        g_I     = g_I_23Na
        J       = J_23Na
        g_J     = g_J_23Na
        A       = A_23Na
        a_s     = a_s_23Na
        alpha_E = alpha_E_23Na
        Fup     = Fup_23Na
        Fdown   = Fdown_23Na     
        INFO    = 1  ! ITS AN ATOM
        ID_name = ID_name_23Na
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

        !write(*,*) "# spin",Fup,Fdown
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
           !write(*,*) Ftotal,total_states_lsi,gF
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
     !write(*,*) JTOTAL,total_states_lsi
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
     MODULE PROCEDURE FLOQUETINIT_QUBIT, FLOQUETINIT_SPIN,FLOQUETINIT_ALKALI,FLOQUETINIT_MBH
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
    write(*,*)

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
    !write(*,*)"Floquet_qubit in floquet_init.f90 ", atomicspecie,info,id%d_bare
    !call write_matrix(J_x)
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

    !write(*,*) atomicspecie,JTOTAL
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

  SUBROUTINE FLOQUETINIT_MBH(ID,atomicspecie,NP,L_,stats,info)
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
    INTEGER,           intent(in)    :: NP,L_ ! 
    CHARACTER (LEN=*), INTENT(IN)    :: stats
    INTEGER,           INTENT(INOUT) :: INFO

    INTEGER  r,D_F2,P,r_,p_
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: Energy

    !  write(*,*) atomicspecie
    INFO = 6

    TOTAL_STATES_LSI =  1.0!D_H(L_,NP,stats)
    CALL SET_ATOMIC_PARAMETERS('lattice','B',1.0D0*total_states_lsi,ID,INFO)

    !------ ALLOCATE NEEDED ARRAYS: Hamiltonian and Lapack
    ALLOCATE(Energy(TOTAL_STATES_LSI))
    ALLOCATE(HAMILTONIAN(Total_states_LSI,Total_states_LSI))

  END SUBROUTINE FLOQUETINIT_MBH


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
    !write(*,*) total_states_lsi,info
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

  TOTAL_FREQUENCIES = NF

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
  DOUBLE PRECISION,  INTENT(IN),optional :: JTOTAL
  CHARACTER (LEN=1), INTENT(IN),optional :: MANIFOLD  

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
  ALLOCATE(HAMILTONIAN(Total_states_LSI,Total_states_LSI))

  U_zeeman = 0

  IF(INFO.EQ.1 .AND. PRESENT(MANIFOLD)) THEN
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


