
MODULE ATOMIC_PROPERTIES
  USE physical_constants
  IMPLICIT NONE
  ! Parameters values by default
  DOUBLE PRECISION :: L = 0.0
  DOUBLE PRECISION :: S = 0.5
  DOUBLE PRECISION :: J,F,gf,mf  
  DOUBLE PRECISION :: gF_2,gF_1,G_F
  INTEGER          :: Total_states_LSI
  INTEGER          :: Fup,Fdown
  DOUBLE PRECISION :: mass_at,Ftotal,I,g_I,g_J,A,a_s,alpha_E
  CHARACTER(LEN=7) :: ID_name
   
  !87Rb
  DOUBLE PRECISION, PARAMETER :: mass_at_87Rb = 87*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_87Rb  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_87Rb       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_87Rb       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_87Rb     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_87Rb     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_87Rb       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_87Rb     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_87Rb = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_87Rb     =  2
  INTEGER,          PARAMETER :: Fdown_87Rb   =  1
  CHARACTER(LEN=7) :: ID_name_87Rb = "87Rb"
 
  !85Rb
  DOUBLE PRECISION, PARAMETER :: mass_at_85Rb = 85*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_85Rb  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_85Rb       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_85Rb       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_85Rb     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_85Rb     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_85Rb       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_85Rb     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_85Rb = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_85Rb     =  2
  INTEGER,          PARAMETER :: Fdown_85Rb   =  1
  CHARACTER(LEN=7) :: ID_name_85Rb = "85Rb"


  !6Li  
  DOUBLE PRECISION, PARAMETER :: mass_at_6Li = 6*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6Li  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_6Li       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6Li       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6Li     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6Li     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6Li       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6Li     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6Li = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6Li     =  2
  INTEGER,          PARAMETER :: Fdown_6Li   =  1
  CHARACTER(LEN=7) :: ID_name_6Li = "6Li"


  !Cs
  DOUBLE PRECISION, PARAMETER :: mass_at_6Cs = 87*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6Cs  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_Cs       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6Cs       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6Cs     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6Cs     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6Cs       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6Cs     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6Cs = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6Cs     =  2
  INTEGER,          PARAMETER :: Fdown_6Cs   =  1
  CHARACTER(LEN=7) :: ID_name_Cs = "Cs"


  !41K
  DOUBLE PRECISION, PARAMETER :: mass_at_6K = 41*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6K  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_K       =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6K       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6K     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6K     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6K       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6K     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6K = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6K     =  2
  INTEGER,          PARAMETER :: Fdown_6K   =  1
  CHARACTER(LEN=7) :: ID_name_41K = "41K"

  !Na
  DOUBLE PRECISION, PARAMETER :: mass_at_6Na = 87*amu
  DOUBLE PRECISION, PARAMETER :: Ftotal_6Na  = 2.0
  DOUBLE PRECISION, PARAMETER :: J_Na        =  0.5  
  DOUBLE PRECISION, PARAMETER :: I_6Na       =  1.5  
  DOUBLE PRECISION, PARAMETER :: g_J_6Na     =  2.0
  DOUBLE PRECISION, PARAMETER :: g_I_6Na     = -0.000995
  DOUBLE PRECISION, PARAMETER :: A_6Na       =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION, PARAMETER :: a_s_6Na     = 5.77E-9
  DOUBLE PRECISION, PARAMETER :: alpha_E_6Na = 2*pi*hbar*0.0794*1E-4
  INTEGER,          PARAMETER :: Fup_6Na     =  2
  INTEGER,          PARAMETER :: Fdown_6Na   =  1
  CHARACTER(LEN=7) :: ID_name_Na = "Na"


  !Qubit
  DOUBLE PRECISION, PARAMETER :: mass_at_qubit =  amu
  DOUBLE PRECISION :: Ftotal_qubit  =  0.0
  DOUBLE PRECISION :: J_qubit       =  0.0
  DOUBLE PRECISION :: I_qubit       =  0.0
  DOUBLE PRECISION, PARAMETER :: g_J_qubit     =  1.0
  DOUBLE PRECISION, PARAMETER :: g_I_qubit     =  0.0
  DOUBLE PRECISION, PARAMETER :: A_qubit       =  1.0
  DOUBLE PRECISION, PARAMETER :: a_s_qubit     =  0.0
  DOUBLE PRECISION, PARAMETER :: alpha_E_qubit =  0.0
  INTEGER,          PARAMETER :: Fup_qubit     =  0
  INTEGER,          PARAMETER :: Fdown_qubit   =  0.5
  CHARACTER(LEN=7) :: ID_name_qubit = "qubit"

  !spin
  DOUBLE PRECISION :: I_spin   =  0.0
  DOUBLE PRECISION :: J_spin   =  0.0  
  DOUBLE PRECISION :: gJ_spin  =  1.0
  DOUBLE PRECISION :: gI_spin  =  0.0
  DOUBLE PRECISION :: A_spin   =  1.0
  DOUBLE PRECISION :: a_s_spin =  0.0
  DOUBLE PRECISION :: alpha_E_spin = 0.0
  INTEGER          :: Fup_spin     =  1
  INTEGER          :: Fdown_spin   =  1
  CHARACTER(LEN=7) :: ID_name_spin = "spin"


  !lattice
  CHARACTER        :: PERIODIC      
  CHARACTER(LEN=7) :: ID_name_lattice = "lattice"

END MODULE ATOMIC_PROPERTIES
