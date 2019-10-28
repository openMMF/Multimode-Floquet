SUBROUTINE SET_ATOMIC_PARAMETERS_C(ATOMICSPECIE,MANIFOLD,JTOTAL,ID_C,INFO)
  ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
  ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
  ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
  !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES
  USE TYPES_C
  USE MODES_4F
  IMPLICIT NONE
  CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: ATOMICSPECIE
  CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: MANIFOLD  !
  INTEGER,          OPTIONAL, INTENT(IN) :: JTOTAL
  TYPE(ATOM_C),INTENT(OUT) :: ID_C
  INTEGER, INTENT(INOUT) :: INFO
  INFO = 0

  CALL SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,1.0D0*JTOTAL,ATOM_,INFO)

END SUBROUTINE SET_ATOMIC_PARAMETERS_C



SUBROUTINE DEALLOCATEALL_C(ID)
  USE ARRAYS
  IMPLICIT NONE
  INTEGER, INTENT(IN)::ID
  CALL DEALLOCATEALL(ID)
END SUBROUTINE DEALLOCATEALL_C




!MODULE FLOQUETINITINTERFACE_C
!  INTERFACE FLOQUETINIT_C
!     MODULE PROCEDURE FLOQUETINIT_QUBIT_C, FLOQUETINIT_SPIN_C,FLOQUETINIT_ALKALI_C
!  END INTERFACE FLOQUETINIT_C
!CONTAINS
  SUBROUTINE FLOQUETINIT_QUBIT_C(ID_C,length_name,ATOMICSPECIE,INFO)
    ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
    ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
    ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
    !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES
    ! calculate the dimenson of the Hilbert space
    ! initialize all the matrices required for a full Floquet calcuations
    ! Calculate the nuclear, electron and total angular momentum operators
    
    USE TYPES_C
    USE MODES_4F
    USE FLOQUETINITINTERFACE
    IMPLICIT NONE


    CHARACTER(length_name), INTENT(IN)    :: ATOMICSPECIE
    INTEGER,                    INTENT(IN)    :: length_name!,JTOTAL
    !DOUBLE PRECISION,           INTENT(IN)    :: JTOTAL
    TYPE(ATOM_C),               INTENT(OUT)   :: ID_C
    INTEGER,                    INTENT(INOUT) :: INFO

    CHARACTER(length_name) atomicspecie_F

    atomicspecie_F = atomicspecie(1:length_name)
    
    !CALL FLOQUETINIT_OLD(ATOM_,atomicspecie_F,manifold,1.0D0*JTOTAL,info)
    CALL FLOQUETINIT_QUBIT(ATOM_,atomicspecie_F,INFO)
    
    ID_C%id_system = ATOM_%id_system
    ID_C%D_BARE    = ATOM_%D_BARE

  END SUBROUTINE FLOQUETINIT_QUBIT_C

  SUBROUTINE FLOQUETINIT_SPIN_C(ID_C,length_name,atomicspecie,jtotal,info)
    ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
    ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
    ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
    !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES


    ! calculate the dimenson of the Hilbert space
    ! initialize all the matrices required for a full Floquet calcuations
    ! Calculate the nuclear, electron and total angular momentum operators

    USE TYPES_C
    USE MODES_4F
    USE FLOQUETINITINTERFACE
    IMPLICIT NONE


    CHARACTER(length_name), INTENT(IN)    :: ATOMICSPECIE
    INTEGER,                    INTENT(IN)    :: length_name!,JTOTAL
    DOUBLE PRECISION,           INTENT(IN)    :: JTOTAL
    TYPE(ATOM_C),               INTENT(OUT)   :: ID_C
    INTEGER,                    INTENT(INOUT) :: INFO

    CHARACTER(length_name) atomicspecie_F
    
    atomicspecie_F = atomicspecie(1:length_name)
    CALL FLOQUETINIT_SPIN(ATOM_,atomicspecie_F,jtotal,info)
    
    ID_C%id_system = ATOM_%id_system
    ID_C%D_BARE    = ATOM_%D_BARE



  END SUBROUTINE FLOQUETINIT_SPIN_C

  SUBROUTINE FLOQUETINIT_ALKALI_C(ID_C,length_name,atomicspecie,length_name2,manifold,info)
    ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
    ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
    ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
    !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES


    ! calculate the dimenson of the Hilbert space
    ! initialize all the matrices required for a full Floquet calcuations
    ! Calculate the nuclear, electron and total angular momentum operators

    USE TYPES_C
    USE MODES_4F
    USE FLOQUETINITINTERFACE
    IMPLICIT NONE


    CHARACTER(length_name), INTENT(IN)    :: ATOMICSPECIE
    CHARACTER(length_name2), INTENT(IN)    :: manifold
    INTEGER,                    INTENT(IN)    :: length_name,length_name2
    TYPE(ATOM_C),               INTENT(OUT)   :: ID_C
    INTEGER,                    INTENT(INOUT) :: INFO

    CHARACTER(length_name) atomicspecie_F
    CHARACTER(length_name2) manifold_F
     
    atomicspecie_F = atomicspecie(1:length_name)
    manifold_F     = manifold(1:length_name2)
    CALL FLOQUETINIT_ALKALI(ATOM_,atomicspecie_F,manifold,info)
    
    ID_C%id_system = ATOM_%id_system
    ID_C%D_BARE    = ATOM_%D_BARE

  END SUBROUTINE FLOQUETINIT_ALKALI_C

!END MODULE FLOQUETINITINTERFACE_C

SUBROUTINE FLOQUETINIT_OLD_C(length_name,atomicspecie,manifold,JTOTAL,ID_C,info)
  ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
  ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
  ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
  !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES


  ! calculate the dimenson of the Hilbert space
  ! initialize all the matrices required for a full Floquet calcuations
  ! Calculate the nuclear, electron and total angular momentum operators

  USE TYPES_C
  USE MODES_4F
  USE FLOQUETINITINTERFACE
  IMPLICIT NONE


  CHARACTER(length_name), INTENT(IN)    :: ATOMICSPECIE
  CHARACTER(1),           INTENT(IN)    :: MANIFOLD  !
  INTEGER,                    INTENT(IN)    :: length_name!,JTOTAL
  DOUBLE PRECISION,           INTENT(IN)    :: JTOTAL
  TYPE(ATOM_C),               INTENT(OUT)   :: ID_C
  INTEGER,                    INTENT(INOUT) :: INFO



  CHARACTER(length_name) atomicspecie_F

  atomicspecie_F = atomicspecie(1:length_name)

  CALL FLOQUETINIT_OLD(ATOM_,atomicspecie_F,manifold,1.0D0*JTOTAL,info)

  ID_C%id_system = ATOM_%id_system
  ID_C%D_BARE    = ATOM_%D_BARE


END SUBROUTINE FLOQUETINIT_OLD_C


SUBROUTINE SETHAMILTONIANCOMPONENTS_C(ATOM__C,NM,NF,MODES_NUM,COUPLING_C,INFO)
  ! ID  tYPE OF ATOM
  ! MODES_NUM, VECTOR. THE SIZE OF THE VECTOR TELL US THE NUMBER OF FREQUENCIES, AND THE VALUE OF EACH COMPONENT INDICATES THE NUMBER OF HARMONICS OF EACH FREQUENCI
  ! FIELDS : IN AND OUTPUT THE MATRICES
  ! INFO

  !USE ARRAYS

  USE TYPES_C
  USE MODES_4F

  IMPLICIT NONE
  INTEGER,                     INTENT(IN)    :: NM,NF
  TYPE(ATOM_C),                INTENT(IN)    :: ATOM__C
  INTEGER,      DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE_C), DIMENSION(NF), INTENT(INOUT) :: COUPLING_C
  INTEGER,                     INTENT(INOUT) :: INFO


  CALL COUPLINGINIT_C(ATOM__C%D_BARE,NF,ATOM__C,COUPLING_C,INFO) ! IN THIS SUBROUTINE WE USE THE C STRUCTURES ATOM__C AND COUPLING_C
  ! TO DEFINE AND INITIALISE THE GLOBAL FORTRAN TYPES ATOM_ AND COUPLING
  ! THESE TWO TYPES ARE THEN PASSED ON TO THE ALL OTHER ROUTINES
  !write(*,*) info,nm,nf
  CALL SETHAMILTONIANCOMPONENTS(ATOM_,NM,NF,MODES_NUM,COUPLING,INFO)

END SUBROUTINE SETHAMILTONIANCOMPONENTS_C

