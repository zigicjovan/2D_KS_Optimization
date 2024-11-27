!===================================================================================================
! MODULE CONTAINS GLOBAL VARIABLES TO BE DECLARED TO RUN THE FORWARD NAVIER-STOKES SOLVER
!
! Author: Pritpal Matharu                             
! Department of Mathematics and Statistics            
! McMaster University                                 
! Date: 2020/12/24                                    
!
!===================================================================================================
MODULE global_variables
  USE, INTRINSIC :: iso_c_binding ! Standard module to define equivalents of C types
  IMPLICIT NONE                   ! Prevent using implicit typing rules

  ! Define parameters to use throughout codes (Note: variables that are initialized when declared have an implicit save attribute)
  INTEGER, PARAMETER          :: pr           = KIND (0.0d0)              ! Integer label for compiler, used to represent double precision
  REAL, PARAMETER             :: MACH_EPSILON = 1.0e-16                   ! Define machine epsilon
  REAL(pr), PARAMETER         :: PalinIV      = 0.0_pr                 ! Value of initial cost functional at time 0 (for H1 semi norm)
  REAL(pr), PARAMETER         :: iniTime      = 0.0_pr, endTime = 10.0_pr  ! Initial and final time
  CHARACTER(len=*), PARAMETER :: normconstr   = "H1semi"                  ! Type of norm constraint to enforce on problem
  CHARACTER(len=*), PARAMETER :: Grad_type    = "H1"                      ! Type of gradient used in optimization scheme
  REAL, PARAMETER             :: visc    = 1.0     ! Kinematic viscosity
  INTEGER, PARAMETER          :: RESOL   = 512       ! Number of discretization points in one direction
  REAL(pr), PARAMETER         :: dt      = 1.0e-3  ! Time step size
  REAL(pr), PARAMETER         :: ell     = 1.0_pr     ! Sobolev parameter for H1 Gradient
  CHARACTER(len=*), PARAMETER :: IC_type = "sineKS"    ! Type of initial vorticity to use (sinusoidal)
  !CHARACTER(len=*), PARAMETER :: IC_type = "gaussianKS"       ! Type of initial vorticity to use (gaussian)

  INTEGER, PARAMETER          :: RESOLP   = RESOL     ! Number of discretization points from previous optimization for bootstrapping
  REAL, PARAMETER             :: viscP    = visc      ! Viscosity from previous optimization for bootstrapping
  REAL, PARAMETER             :: endTimeP = endTime   ! Final time from previous optimization for bootstrapping
  REAL(pr), PARAMETER         :: ellP     = ell       ! Sobolev parameter from previous optimization for bootstrapping
  LOGICAL, PARAMETER          :: BS_flag  = .FALSE.   ! Flag to indicate if this is more than first bootstrap (first bs, set to false; more, set to true)

  ! Initialize variables
  LOGICAL                        :: diagn_flag             ! Flag to indicate if code should save data
  LOGICAL                        :: bin_flag               ! Flag to indicate if code should save bin files
  LOGICAL                        :: vid_flag               ! Flag to indicate if code should data should be saved for video
  REAL(pr), SAVE                 :: PI                     ! Value of pi
  REAL(pr), SAVE                 :: T1                     ! Length of time window
  REAL(pr), SAVE                 :: Lx, Ly                 ! Length of domain in x and y
  REAL(pr), SAVE                 :: dV                     ! Volume step size
  REAL(pr), SAVE                 :: Kcut                   ! Frequency cutoff for dealiasing
  INTEGER                        :: Time_iter              ! Counter for number of time iterations
  INTEGER,  SAVE                 :: numel_T                ! Number of time steps in time window
  INTEGER,  SAVE                 :: dtnp                   ! Number of to divide up the time window for processors 
  INTEGER,  SAVE                 :: ny_dim                 ! Number of elements of array on LOCAL processor
  INTEGER,  SAVE                 :: nx_dim                 ! Number of elements of array on LOCAL processor
  INTEGER,  SAVE                 :: nh                     ! Number of discretization points for x stored in Fourier space
  INTEGER,  DIMENSION(2),   SAVE :: n_nse                  ! Number of discretization points for (x, y)
  REAL(pr), DIMENSION(2),   SAVE :: dx                     ! Step size in spatial directions in (x, y)
  REAL(pr), DIMENSION(1:4), SAVE :: Alpha, BetaI, BetaE, Gamma    ! Values for the IMEX Method

  REAL(pr), DIMENSION (:),     ALLOCATABLE, SAVE :: K1, K2 ! Wavenumbers in x and y
  REAL(pr), DIMENSION (:),     ALLOCATABLE, SAVE :: t_vec  ! Time vector
  REAL(pr), DIMENSION (:),     ALLOCATABLE, SAVE :: Enst   ! Enstrophy vector
  REAL(pr), DIMENSION (:),     ALLOCATABLE, SAVE :: KinEn  ! Kinetic Energy vector
  REAL(pr), DIMENSION (:),     ALLOCATABLE, SAVE :: Palin  ! Palinstrophy vector
  REAL(pr), DIMENSION (:,:),   ALLOCATABLE, SAVE :: vort0  ! Vorticity Field
  REAL(pr), DIMENSION (:,:),   ALLOCATABLE, SAVE :: ksq    ! Vorticity Field (laplace operator)
  REAL(pr), DIMENSION (:,:),   ALLOCATABLE, SAVE :: lin_hat    ! Vorticity Field (linear term)

  ! Directories for saving
  CHARACTER(len=*), PARAMETER :: work_pathname    = "/home/zigicj/2DKS_optimization/2D_KS_Optimization/2DKS_data/"
  CHARACTER(len=*), PARAMETER :: scratch_pathname = "/home/zigicj/2DKS_optimization/2D_KS_Optimization/bin_files/"

  ! Variables for saving
  CHARACTER(4)  :: Nchar        ! Resolution as character
  CHARACTER(8)  :: lchar        ! Sobolev parameter as character
  CHARACTER(10) :: tchar        ! Final time as character
  CHARACTER(10) :: viscchar     ! Viscosity as character

  !==========================================================================
  !                            MPI VARIABLES
  !==========================================================================
  INTEGER, SAVE :: rank, Statinfo, np ! Local processor rank, error flag, and number of processors

  !==========================================================================
  !                            FFTW VARIABLES
  !==========================================================================
  INTEGER(C_INTPTR_T), DIMENSION(2), SAVE :: C_n                             ! Number of points in x and y (defined for C interface, to use with FFTW)
  INTEGER(C_INTPTR_T), SAVE :: C_nh                                          ! Half number of points in x (defined for C interface, to use with FFTW)
  INTEGER(C_INTPTR_T), SAVE :: C_local_y_alloc, C_local_Ny, C_local_y_offset ! Local array data size on processor, potion of y data on processor, and reference index for correct positioning (for C interface)
  INTEGER(C_INTPTR_T), SAVE :: C_local_x_alloc, C_local_Nx, C_local_x_offset ! Local array data size on processor, potion of x data on processor, and reference index for correct positioning (for C interface)
  INTEGER,             SAVE :: local_Ny, local_y_offset                      ! Portion of y data on processor and reference index for correct positioning
  INTEGER,             SAVE :: local_Nx, local_x_offset                      ! Portion of x data on processor and reference index for correct positioning

  TYPE(C_PTR),         SAVE :: fwdplan, bwdplan, tmpdataR, tmpdataC          ! Forward plan, backward plan, and temporary memory allocation for C interface


END MODULE
