MODULE TYPES

  TYPE :: MODE
     DOUBLE PRECISION :: OMEGA
     COMPLEX*16       :: X,Y,Z
     DOUBLE PRECISION :: phi_x,phi_y,phi_z
     INTEGER          :: N_Floquet
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: V
     COMPLEX*16, DIMENSION(:),   ALLOCATABLE :: VALUES
     INTEGER,    DIMENSION(:),   ALLOCATABLE :: ROW,COLUMN
  END TYPE MODE
  
  TYPE :: ATOM
     INTEGER          :: id_system
     INTEGER          :: D_BARE
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_BARE
  END TYPE ATOM

  TYPE :: HARMONIC_FACTORS
     COMPLEX*16,DIMENSION(:,:), ALLOCATABLE :: U,U_r,U_AVG
     INTEGER,   DIMENSION(:),   ALLOCATABLE :: n
  END type HARMONIC_FACTORS

END MODULE TYPES



MODULE subinterface
!!$  IMPLICIT NONE
!!$  INTERFACE
!!$     SUBROUTINE SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,JTOTAL,ID,INFO)
!!$     !SUBROUTINE SET_ATOMIC_PARAMETERS(ATOMICSPECIE,MANIFOLD,JTOTAL,JTOTAL_HALF,ID,INFO)
!!$       ! ATOMICSPECIE: 87Rb,6Li,Cs,41K,qubit,lattice, SPIN
!!$       ! MANIFOLD : "U" UPPER HYPERFINE MANIFOLD, "L" LOWER HYPERFIND MANIFOLD, "B" BOTH
!!$       ! JTOTAL   :  IF ATOMICSPECIE .EQ. SPIN THEN JTOTAL IS THE TOTAL ANGULAR MOMENTUM OF THE SPIN
!!$       !             IF ATOMICSPECIE .EQ. LATTICE, THEN JTOTAL IS THE NUMBER OF SITES
!!$       USE TYPES
!!$       IMPLICIT NONE
!!$       CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: ATOMICSPECIE
!!$       CHARACTER (LEN=*),OPTIONAL, INTENT(IN) :: MANIFOLD  !
!!$       !INTEGER,          OPTIONAL, INTENT(IN) :: JTOTAL
!!$       DOUBLE PRECISION,          OPTIONAL, INTENT(IN) :: JTOTAL
!!$       TYPE(ATOM),OPTIONAL,INTENT(OUT) :: ID       
!!$       INTEGER, INTENT(INOUT) :: INFO
!!$     END SUBROUTINE SET_ATOMIC_PARAMETERS
!!$  END INTERFACE
  implicit none
  interface
     subroutine mesh_exp(r_min, r_max, a, N,mesh) 
       use types
       double precision, intent(in), value :: r_min
       double precision, intent(in), value :: r_max
       double precision, intent(in), value :: a
       integer, intent(in), value :: N
       !type(atom), intent(in) :: id
       double precision, dimension(N),intent(out) :: mesh 
     END subroutine mesh_exp
  end interface
END MODULE subinterface

subroutine mesh_exp(r_min, r_max, a, N,mesh) 
  USE TYPES
  double precision, intent(in), value :: r_min
  double precision, intent(in), value :: r_max
  double precision, intent(in), value :: a
  integer, intent(in), value :: N
  !TYPE(ATOM) :: ID
  double precision, dimension(N),intent(out) :: mesh 
  
  integer i;
  
  do i=1,N
     mesh(i) = a + (i-1)*(r_max-r_min)/N
  end do
  
end subroutine mesh_exp


subroutine mesh_exp_c(id,r_min, r_max, a,N,mesh) !bind(c,name='mesh_exp_c')
  USE SUBINTERFACE
  USE ISO_C_BINDING
  USE TYPES
  
  IMPLICIT NONE
  real(c_double), intent(in), value :: r_min
  real(c_double), intent(in), value :: r_max
  real(c_double), intent(in), value :: a
  integer(c_int), intent(in), value :: N
  TYPE(ATOM),intent(inout) :: ID
  real(c_double), intent(out) :: mesh(N)
  
  call mesh_exp(r_min,r_max,a,N,mesh)
  allocate(ID%E_bare(7))
  !ID%E_bare=0
  !ID%E_bare(3) = 4
  !write(*,*) ID%E_bare
end subroutine mesh_exp_c

subroutine derivedType_c(INFO) bind(c, name='derivedType_c')
  USE ISO_C_BINDING
  USE TYPES
  
  IMPLICIT NONE
  !TYPE(ATOM),intent(out) :: ID
  INTEGER(c_int), intent(in) :: INFO
  !ID%id_system = 1
  !ID%D_BARE = 2
  !ALLOCATE(ID%E_BARE(2))
  !info = 0
  write(*,*) info

end subroutine derivedType_c


!!$PROGRAM TEST
!!$  
!!$  
!!$!  USE ATOMIC_PROPERTIES
!!$  USE TYPES
!!$  USE SUBINTERFACE
!!$  
!!$  IMPLICIT NONE
!!$!  TYPE(ATOM) ID
!!$  INTEGER INFO
!!$  
!!$  !CALL SET_ATOMIC_PARAMETERS('87Rb','L',0.2D1,ID,INFO)
!!$  double precision :: r_min,r_max,a
!!$  integer :: N
!!$  double precision, allocatable, dimension(:) :: mesh
!!$  TYPE(ATOM) :: ID
!!$  !double precision, dimension(8) :: mesh
!!$
!!$  r_min = 1.0
!!$  r_max = 2.0
!!$  a  = 3.0
!!$  N = 8
!!$
!!$  allocate(mesh(N))
!!$  call mesh_exp(r_min,r_max,a,N,mesh)
!!$  write(*,*) mesh
!!$END PROGRAM TEST
