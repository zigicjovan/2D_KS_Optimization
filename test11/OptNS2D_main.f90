!-----------------------------------------------------!
! Fortran code using MPI to solve the optimal IC of   !
! the 2D Navier-Stokes system.                        !
!                                                     !
! -Parallel version using MPI, Real-Complex FFT       !
! -Uses tranposed distributions in Fourier space      !
! -Uses pseudo-spectral Galerkin approach with        !
!   dealiasing and a third-order, four step IMEX time !
!   stepping method                                   !
!                                                     !
! Author: Pritpal Matharu                             !
! Department of Mathematics and Statistics            !
! McMaster University                                 !
! Date: 2021/04/27                                    !
!-----------------------------------------------------!

PROGRAM OptNS2D_main
  USE mpi                 ! Use MPI module (binding works well with fftw libraries)
  USE global_variables    ! Declares variables and defines parameters of the system
  USE fftwfunction        ! Contains FFTW Version 3, the FFTW wrappers, and subroutines for performing the FFT
  USE data_ops            ! Contains functions/subroutines that modify data or read/write data
  USE function_ops        ! Contains functions/subroutines to operate on data
  USE solvers             ! Contains the time stepping methods, to solve the PDES
  USE optimization        ! Contains methods and schemes for solving the optimization problem
  IMPLICIT NONE           ! Prevent using implicit typing rules

  !=============================================
  ! MPI
  !=============================================
  CALL MPI_INIT(Statinfo)                           ! Initialize MPI environment
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,Statinfo)  ! Determine the rank of the calling processor
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,Statinfo)    ! Determine the total number of processors
  CALL FFTW_MPI_INIT()                              ! Initialize MPI-FFTW environment
  !=============================================
  CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)         ! Ensure all processors sync before continuing

  ! Set variables and arrays required to solve the 2D system
  CALL nse_initialize

  ! Flag for saving bin files
  bin_flag   = .TRUE.
  ! Flag for saving diagnostics
  diagn_flag = .TRUE.
  ! Flag for saving files for video (need diagn_flag to be true as well)
  vid_flag   = .TRUE.
  !vid_flag   = .FALSE.

  ! Print parameters to .out file
  IF (rank==0) THEN
    PRINT *, "============================================== "
    PRINT *, "                 Optimization                  "
    PRINT *, "============================================== "
    PRINT'(a,I5)',     " Resolution Nx  = ", n_nse(1)
    PRINT'(a,I5)',     " Resolution Ny  = ", n_nse(2)
    PRINT'(a,ES15.4)', " Time Step Size = ", dt
    PRINT'(a,ES15.4)', " Final Time T   = ", endTime
    PRINT'(a,ES15.4)', " Viscosity      = ", visc
    PRINT'(a,I5)',     " Processors     = ", np
    PRINT *,            "Initial guess  = ", IC_type
    PRINT *,            "Norm Constr    = ", normconstr
    PRINT *,            "Grad Type      = ", Grad_type
    PRINT'(a,ES15.4)', " Sobolev L      = ", ell
    PRINT'(a,ES15.4)', " Prev Sobolev L = ", ellP
    PRINT'(a,ES15.4)', " Prev Time T    = ", endTimeP
  END IF

  ! Ensure all processors sync before continuing
  CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

  ! Create plans for FFT
  CALL init_fft                             ! Found in fftwfunction.f90
  CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo) ! Ensure all processors sync before continuing

  ! Print statement, for testing purposes
  IF (rank==0) THEN
    PRINT *, "============================================== "
    PRINT *, " Successful FFT Initialization. "
  END IF

  ! Create initial condition
  CALL initial_guess(IC_type)                 ! Found in function_ops.f90
  CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)   ! Ensure all processors sync before continuing

  ! Print statement, for testing purposes
  IF (rank==0) THEN
    PRINT *, "============================================== "
    PRINT *, " Successful Initial Condition. "
    PRINT *, "============================================== "
  END IF

  ! Solve the optimal initial condition
  CALL OptScheme(vort0)    ! Found in optimization.f90

  ! Deallocate allocatable arrays
  DEALLOCATE(vort0)
  DEALLOCATE(K1)
  DEALLOCATE(K2)
  DEALLOCATE(t_vec)
  DEALLOCATE(Enst)
  DEALLOCATE(InnerProduct_L2)
  DEALLOCATE(InnerProduct_H1)
  DEALLOCATE(InnerProduct_H2)
  DEALLOCATE(InnerProduct_Hn1)
  DEALLOCATE(KinEn)
  DEALLOCATE(Palin)
  DEALLOCATE(ksq)

  ! Print statement, for testing purposes
  IF (rank==0) THEN
    PRINT *, "============================================== "
    PRINT *, " Successful Deallocation. "
    PRINT *, "============================================== "
  END IF

  ! Deallocate plan and temporary memory for C interface (before cleanup)
  CALL fftw_destroy_plan(fwdplan)
  CALL fftw_destroy_plan(bwdplan)
  CALL fftw_free(tmpdataR)
  CALL fftw_free(tmpdataC)

  ! Get rid of all memory and all resources allocated for FFTW
  CALL fftw_mpi_cleanup()

  ! Terminates the MPI environment
  CALL MPI_FINALIZE (Statinfo)

  ! Print statement, to indicate that everything executed correctly
  IF (rank == 0) then
    PRINT *, "Exit Normally!"
  END IF

END PROGRAM OptNS2D_main

