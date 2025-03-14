!===========================================================================================================
! MODULE CONTAINS ROUTINES TO PERFORM REAL-COMPLEX FFT WITH TRANSPOSED DISTRIBUTIONS IN 2D
!
! Author: Jovan Zigic (inherited from Pritpal Matharu)                                          
! McMaster University                                 
! Date: 2024/12/06                                 
!
! SEE FFTW VERSION 3 OFFICIAL DOCUMENTATION FOR MORE INFORMATION ON FOURIER TRANFORMS USED HERE
! http://www.fftw.org/fftw3_doc/
! Section 7     - Calling FFTW from Modern Fortran
! Section 7.2   - Reversing array dimensions
! Section 6     - Distributed-memory FFTW with MPI
! Section 6.4.3 - Transposed distributions
! Section 6.5   - Multi-dimensional MPI DFTs of Real Data
! Section 6.10  - FFTW MPI Performance Tips
! Section 6.13  - FFTW MPI Fortran Interface

! CONTAINS:
! (*) init_fft - Initializes, and creates forward and backward plans for FFT using transposed distributions
! (*) fftfwd   - Performs forward Fourier transform (real to complex transform)
! (*) fftbwd   - Performs backward Fourier transform (complex to real transform)
!===========================================================================================================
MODULE fftwfunction
  USE, INTRINSIC :: iso_c_binding   ! Standard module to define equivalents of C types!
  USE mpi                           ! Use MPI module (binding works well with fftw libraries)
  IMPLICIT NONE                     ! Prevent using implicit typing rules for entire module
  INCLUDE 'fftw3-mpi.f03'           ! Use FFTW Version 3

  ! Initialize variables for FFT
  REAL(C_DOUBLE),            pointer :: tmppointerR(:,:)             ! Temporary pointer of physical space data for C interface
  COMPLEX(C_DOUBLE_COMPLEX), pointer :: tmppointerC(:,:)             ! Temporary pointer of Fourier space data for C interface

    CONTAINS
      !=======================================
      ! *** Forward Fourier transform ***
      ! Input:  u  - field in physical space
      ! Output: fu - field in Fourier space
      !=======================================
      SUBROUTINE fftfwd(u,fu)
        ! Load variables
        USE global_variables, ONLY: pr, n_nse, local_Ny, local_Nx, fwdplan
        ! Initialize variables
!        REAL(pr),    DIMENSION (1:n_nse(1),1:local_Ny), INTENT(IN)  :: u   ! Field in physical space
!        COMPLEX(pr), DIMENSION (1:n_nse(2),1:local_Nx), INTENT(OUT) :: fu  ! Field in Fourier space
        REAL(pr),    DIMENSION (:,:), INTENT(IN)  :: u                    ! Field in physical space
        COMPLEX(pr), DIMENSION (:,:), INTENT(OUT) :: fu                   ! Field in Fourier space       
        tmppointerR(:,:) = u(:,:)                                         ! Assign physical field in Fortran to a C pointer for FFT      
        CALL fftw_mpi_execute_dft_r2c(fwdplan, tmppointerR, tmppointerC)  ! Compute forward Fourier transform
        fu(:,:) = tmppointerC(:,:)                                        ! Assign Fourier field back to Fortran type from C pointer from FFT
      END SUBROUTINE fftfwd

      !=======================================
      ! *** Inverse Fourier transform ***
      ! Input:  fu - field in Fourier space
      ! Output: u  - field in physical space
      !=======================================
      SUBROUTINE fftbwd(fu,u)
        ! Load variables
        USE global_variables, ONLY: pr, n_nse, local_Ny, local_Nx, bwdplan
        ! Initialize variables
!        COMPLEX(pr), DIMENSION (1:n_nse(2),1:local_Nx), INTENT(IN)  :: fu  ! Field in Fourier space
!        REAL(pr),    DIMENSION (1:n_nse(1),1:local_Ny), INTENT(OUT) :: u   ! Field in physical space
        COMPLEX(pr), DIMENSION (:,:), INTENT(IN)  :: fu  ! Field in Fourier space
        REAL(pr),    DIMENSION (:,:), INTENT(OUT) :: u   ! Field in physical space
        ! Assign Fourier field in Fortran to a C pointer for FFT
        tmppointerC(:,:) = fu(:,:)
        ! Compute backward Fourier transform
        CALL fftw_mpi_execute_dft_c2r(bwdplan, tmppointerC, tmppointerR)
        ! Assign physical field back to Fortran type from C pointer from FFT
        u(:,:) = tmppointerR(:,:)
        ! Unnormalized FFT is computed, so normalize output by number of elements
        u = u/(PRODUCT(REAL(n_nse,pr)))
      END SUBROUTINE fftbwd

      !=======================================
      ! *** Initialize FFT plans ***
      !=======================================
      SUBROUTINE init_fft
        ! Load variables
        USE global_variables, ONLY: pr, C_n, C_nh, C_local_Ny, C_local_Nx, C_local_y_alloc, C_local_x_alloc, fwdplan, bwdplan, tmpdataR, tmpdataC

        ! Allocate local data for FFT (physical space)
        tmpdataR = fftw_alloc_real(2*C_local_y_alloc)
        ! Assign target of C pointer to Fortran pointer
        CALL c_f_pointer(tmpdataR, tmppointerR, [2*C_nh, C_local_Ny])
        ! Allocate local data for FFT (Fourier space)
        tmpdataC = fftw_alloc_complex(C_local_x_alloc)
        ! Assign target of C pointer to Fortran pointer
        CALL c_f_pointer(tmpdataC, tmppointerC, [C_n(2), C_local_Nx])
        ! Create forward and backward plans
        fwdplan = fftw_mpi_plan_dft_r2c_2d(C_n(2), C_n(1), tmppointerR, tmppointerC, MPI_COMM_WORLD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT)) ! Forward plan (fast initialization)
        bwdplan = fftw_mpi_plan_dft_c2r_2d(C_n(2), C_n(1), tmppointerC, tmppointerR, MPI_COMM_WORLD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN) ) ! Backward plan (fast initialization)
!        fwdplan = fftw_mpi_plan_dft_r2c_2d(C_n(2), C_n(1), tmppointerR, tmppointerC, MPI_COMM_WORLD, ior(FFTW_PATIENT, FFTW_MPI_TRANSPOSED_OUT)) ! Forward plan (use greater rigour from Wisdom when computing several FFTs)
!        bwdplan = fftw_mpi_plan_dft_c2r_2d(C_n(2), C_n(1), tmppointerC, tmppointerR, MPI_COMM_WORLD, ior(FFTW_PATIENT, FFTW_MPI_TRANSPOSED_IN) ) ! Backward plan (use greater rigour from Wisdom when computing several FFTs)

      END SUBROUTINE init_fft

END MODULE
