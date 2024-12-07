!===================================================================================================
! INITIALIZE AND SET VARIABLES FOR THE FORWARD KURAMOTO-SIVASHINSKY SOLVER
!
! Author: Jovan Zigic (inherited from Pritpal Matharu)                                          
! McMaster University                                 
! Date: 2024/12/06                                    
!
!===================================================================================================
SUBROUTINE KS_initialize
  USE, INTRINSIC :: iso_c_binding     ! Standard module to define equivalents of C types
  USE mpi                             ! Use MPI module (binding works well with fftw libraries)
  USE global_variables                ! Declares variables and defines parameters of the system
  USE fftwfunction                    ! Use FFTW Version 3
  IMPLICIT NONE                       ! Prevent using implicit typing rules

  ! Temporary integer for loops to create wavenumbers
  INTEGER :: i, i1, i2

  ! Time window
  T1      = endTime - iniTime
  ! Number of time steps in window
  numel_T = int( T1/dt )
  ! Number of time steps to divide window amoung processors
  dtnp    = int( numel_T/np )

  ! Values for the IMEX method
  ! Coefficients from https://doi.org/10.1007/s10898-019-00855-1
  ! Journal of Global Optimization (2020) - Alimo, Cavaglieri, Beyhaghi, Bewley
  AlphaI  = dt*[343038331393.0_pr/1130875731271.0_pr, 288176579239.0_pr/1140253497719.0_pr,  253330171251.0_pr/677500478386.0_pr,   189462239225.0_pr/1091147436423.0_pr]
  BetaI  = dt*[35965327958.0_pr/140127563663.0_pr,   19632212512.0_pr/2700543775099.0_pr,  -173747147147.0_pr/351772688865.0_pr,   91958533623.0_pr/727726057489.0_pr]
  AlphaE  = dt*[14.0_pr/25.0_pr,                      777974228744.0_pr/1346157007247.0_pr,  251277807242.0_pr/1103637129625.0_pr,  113091689455.0_pr/220187950967.0_pr]
  BetaE  = dt*[0.0_pr,                              -251352885992.0_pr/790610919619.0_pr,  -383714262797.0_pr/1103637129625.0_pr, -403360439203.0_pr/1888264787188.0_pr]

  ! Define pi
  PI = 4.0_pr*ATAN2(1.0_pr,1.0_pr)

  ! Set number of discretization points in each direction (defined for C interface, to use with FFTW)
  C_n(1) = RESOL           ! x direction
  C_n(2) = RESOL           ! y direction
  C_nh   = int(C_n(1)/2)+1 ! x direction, in Fourier space

  ! Number of points in each grid direction
  n_nse(1) = RESOL               ! Number of points in the x direction
  n_nse(2) = RESOL               ! Number of points in the y direction
  nh       = int(n_nse(1)/2) + 1 ! Half the number of points in x with padding, for R2C and C2R transforms (using conjugate symmetry)

  ! Size of periodic domain
  !!Lx = (2.0_pr*PI)*0.2_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*0.5_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*0.8_pr                ! Length in the x direction
  Lx = (2.0_pr*PI)*1.0_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*1.2_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*1.5_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*1.8_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*2.0_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*5.0_pr                ! Length in the x direction
  !!Lx = (2.0_pr*PI)*10.0_pr               ! Length in the x direction

  !!Ly = (2.0_pr*PI)*0.5_pr                ! Length in the y direction
  Ly = (2.0_pr*PI)*1.0_pr                ! Length in the y direction
  !!Ly = (2.0_pr*PI)*1.5_pr                ! Length in the y direction
  !!Ly = (2.0_pr*PI)*2.0_pr                ! Length in the y direction
  !!Ly = (2.0_pr*PI)*5.0_pr                ! Length in the y direction
  !!Ly = (2.0_pr*PI)*10.0_pr               ! Length in the y direction

  ! Step sizes in each direction
  dx(1) = Lx/REAL(n_nse(1), pr) ! Step size in x direction
  dx(2) = Ly/REAL(n_nse(2), pr) ! Step size in y direction
  dV    = dx(1)*dx(2)           ! Volume step size

  ! Get local data size and allocate for each calling processor. Note data is C_n(1) x C_n(2), i.e. [X, Y],
  ! but dimension reversal due to the different memory storage convention between C and Fortran
  C_local_y_alloc = fftw_mpi_local_size_2d(C_n(2), C_nh, MPI_COMM_WORLD, C_local_Ny, C_local_y_offset)
  ! Determine number of elements in the y-direction that are on each calling processor
  local_Ny = int( C_local_Ny )
  ! Offset in the y-direction, to obtain correct element in stride
  local_y_offset = int( C_local_y_offset )
  ! Number of elements of from array on LOCAL processor
  ny_dim = n_nse(1)*local_Ny

  ! Get local data size and allocate for each calling processor. Note data is C_nh x C_n(2), since this
  ! this represents the data in Fourier space, where we are using transposed distributions to save a global
  ! transpose
  C_local_x_alloc = fftw_mpi_local_size_2d(C_nh, C_n(2), MPI_COMM_WORLD, C_local_Nx, C_local_x_offset)
  ! Determine number of elements in the x-direction that are on each calling processor
  local_Nx = int( C_local_Nx )
  ! Offset in the y-direction, to obtain correct element in stride
  local_x_offset = int( C_local_x_offset )
  ! Number of elements of from array on LOCAL processor
  nx_dim = n_nse(2)*local_Nx

  ! Allocate storage for arrays that will be dynamically created
  ALLOCATE (K1(1:nh))                       ! Wavenumbers in x (r2c transform)
  ALLOCATE (K2(1:n_nse(2)))                 ! Wavenumbers in y
  ALLOCATE (t_vec(1:numel_T+1))             ! Time vector
  ALLOCATE (KinEn(1:numel_T+1))             ! Kinetic Energy vector
  ALLOCATE (InnerProduct_L2(1:numel_T+1))   ! L^2 inner product vector
  ALLOCATE (InnerProduct_H1(1:numel_T+1))   ! H^1 inner product vector
  ALLOCATE (InnerProduct_H2(1:numel_T+1))   ! H^2 inner product vector
  ALLOCATE (InnerProduct_Hn1(1:numel_T+1))  ! H^(-1) inner product vector
  ALLOCATE (Enst(1:numel_T+1))              ! Enstrophy vector
  ALLOCATE (Palin(1:numel_T+1))             ! Palinstrophy vector
  ALLOCATE (vort0(1:n_nse(1), 1:local_Ny))  ! Vorticity field
  ALLOCATE (ksq(1:n_nse(2),1:local_Nx))     ! Vorticity field (laplace operator)
  ALLOCATE (lin_hat(1:n_nse(2),1:local_Nx))     ! Vorticity field (linear term)

  ! Set up vectors
  ! Wavenumbers in x
!  DO i = 0, n_nse(1)-1
  DO i = 0, nh-1
    IF (i<=n_nse(1)/2) THEN
      K1(i+1) = ((2.0_pr*PI)/Lx)*REAL(i,pr) 
    ELSE
      K1(i+1) = ((2.0_pr*PI)/Lx)*REAL(i-n_nse(1),pr)
    END IF
  END DO
  ! Wavenumbers in y
  DO i = 0,n_nse(2)-1
    IF (i <= n_nse(2)/2) THEN
      K2(i+1) = ((2.0_pr*PI)/Ly)*REAL(i,pr) 
    ELSE
      K2(i+1) = ((2.0_pr*PI)/Ly)*REAL(i-n_nse(2),pr)
    END IF
  END DO
  ! Local array, for the squared magnitude of the wavevector
  DO i2=1,local_Nx
    DO i1=1,n_nse(2)
      ksq(i1, i2) = ( K1(i2+local_x_offset)**2 + K2(i1)**2 ) ! 2DNS linear term + Laplace operator
      lin_hat(i1, i2) = ( ( K1(i2+local_x_offset)**2 + K2(i1)**2 ) - ( K1(i2+local_x_offset)**4 + K2(i1)**4 )) ! 2DKS linear term
    END DO
  END DO
  ! Time vector
  DO i = 0, numel_T
    t_vec(i+1) = i*dt
  END DO
  ! Frequency cutoff for dealiasing
  Kcut = SQRT( K1(n_nse(1)/2+1)**2 + K2(n_nse(2)/2+1)**2 ) * (2.0_pr/3.0_pr) ! 2DNS

  ! Resolution as a character
  WRITE(Nchar, '(i4.4)') RESOL
  ! Viscosity as a character
  WRITE(viscchar, '(ES10.2)') visc
  ! Final time as a character
  WRITE(tchar, '(ES10.2)') endTime
  ! Sobolev parameter as a character
  WRITE(lchar, '(ES8.1)') ell

END SUBROUTINE
