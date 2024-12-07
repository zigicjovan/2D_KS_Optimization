!==================================================================================================
! MODULE CONTAINING INTERFACES FOR KURAMOTO-SIVASHINSKY DIRECT NUMERICAL SIMULATION SOLVERS
!
! Author: Jovan Zigic (inherited from Pritpal Matharu)                                      
! McMaster University                                 
! Date: 2024/12/06                                    
!
! CONTAINS:
! (*) fwd_NS2D               - Forward DNS of 2D NS
! (*) fwd_2DKS               - Forward DNS of 2D KS
! (*) IMEX_fwd               - Third-order, four step IMEX time stepping method
! (*) IMEX_bwd               - Third-order, four step IMEX time stepping method for adjoint system
!==================================================================================================
MODULE solvers
  IMPLICIT NONE ! Prevent using implicit typing rules for entire module

  CONTAINS
    !=================================================================================
    ! Numerical DNS solver of 2D Navier Stokes Forward System using a pseudo-spectral
    ! Galerkin approach with dealiasing and a third order IMEX time stepping method
    ! Input: w0 - Initial condition
    !=================================================================================
    SUBROUTINE fwd_NS2D(w0)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Ny, Palin, Enst, KinEn, t_vec, rank, Time_iter, Statinfo, InnerProduct_L2, InnerProduct_H1, InnerProduct_H2, InnerProduct_Hn1
      ! Load subroutines
      USE data_ops          ! Contains functions/subroutines that read/write data
      USE mpi               ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
      REAL(pr), DIMENSION(:,:), INTENT(IN) :: w0 ! Vorticity in physical space

      ! Save initial vorticity field
      Time_iter = 0 ! Time iteration counter
      CALL save_NS_vorticity(w0, Time_iter, "FWD")
      ! Solve forward 2D DNS
      IF (rank==0) PRINT *, " Solving Forward DNS..."
      Time_iter = 1 ! Time iteration counter
      CALL IMEX_fwd(w0, Palin)
      IF (rank==0) PRINT *, " DNS solved."
      ! Save vorticity and diagnostics for MATLAB analysis
      CALL save_NS_DNS(w0, Palin, Enst, KinEn, t_vec, InnerProduct_L2, InnerProduct_H1, InnerProduct_H2, InnerProduct_Hn1)

      ! Ensure all processes have completed
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
    END SUBROUTINE fwd_NS2D

    !=================================================================================
    ! Numerical DNS solver of 2D Navier Stokes Forward System using a pseudo-spectral
    ! Galerkin approach with dealiasing and a third order IMEX time stepping method
    ! Input: w0 - Initial condition
    !=================================================================================
    SUBROUTINE fwd_2DKS(w0)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Ny, Palin, Enst, KinEn, t_vec, rank, Time_iter, Statinfo, InnerProduct_L2, InnerProduct_H1, InnerProduct_H2, InnerProduct_Hn1
      ! Load subroutines
      USE data_ops          ! Contains functions/subroutines that read/write data
      USE mpi               ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
      REAL(pr), DIMENSION(:,:), INTENT(IN) :: w0 ! Vorticity in physical space

      ! Save initial vorticity field
      Time_iter = 0 ! Time iteration counter
      CALL save_NS_vorticity(w0, Time_iter, "FWD")
      ! Solve forward 2D DNS
      IF (rank==0) PRINT *, " Solving Forward DNS..."
      Time_iter = 1 ! Time iteration counter
      CALL IMEX_fwd(w0, Palin)
      IF (rank==0) PRINT *, " DNS solved."
      ! Save vorticity and diagnostics for MATLAB analysis
      CALL save_NS_DNS(w0, Palin, Enst, KinEn, t_vec, InnerProduct_L2, InnerProduct_H1, InnerProduct_H2, InnerProduct_Hn1)

      ! Ensure all processes have completed
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
    END SUBROUTINE fwd_2DKS

    !=================================================================================
    ! Forward solver with third order IMEX Time stepping
    ! Input: w0 - vorticity field (initial condition) in physical space
    !=================================================================================
    SUBROUTINE IMEX_fwd(w0, palins)
      ! Load variables
      USE global_variables
      ! Load subroutines
      USE fftwfunction        ! FFT routines
      USE function_ops        ! Contains functions/subroutines to operate on data
      USE data_ops            ! Contains functions/subroutines that read/write data
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny), INTENT(IN) :: w0             ! Initial vorticity field
      REAL(pr),    DIMENSION(:,:),       INTENT(IN)          :: w0             ! Initial vorticity field
      REAL(pr),    DIMENSION(:),        INTENT(OUT)          :: palins         ! Palinstrophy of the system
      REAL(pr),    DIMENSION(1:n_nse(1),1:local_Ny)          :: w1             ! Final vorticity
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)          :: w_hat, w2_hat  ! Fourier transform of complex vorticity and RK storage term
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)          :: psi_hat        ! Streamfunction in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)          :: conv_hat       ! Convection term in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)          :: conv0_hat      ! Storage for RK convection term in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx, 2)       :: u_hat          ! Fourier transform of velocity field ! 2DNS
      !COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)       :: u_hat          ! Fourier transform of velocity field ! 2DKS
      INTEGER                                                :: rk, i1, i2     ! Integers for looping through values
      INTEGER                                                :: Ni, Nsave             ! Integer for storing diagn_flag save position
      !REAL(pr)                                               :: Nsave          ! Scalar for storing diagn_flag save interval
      REAL(pr)                                               :: mean_val       ! Scalar for storing the mean of vorticity
      REAL(pr) :: local_kin, local_enst, local_palins ! Temporary for storing values of local kinetic, local enstrophy, and palinstrophy
      REAL(pr) :: local_L2, local_H1, local_H2, local_Hn1 ! Temporary for storing values of local inner products
      
      IF (vid_flag) THEN ! If saving video flag true, start counter
        Ni = 1
      END IF

      CALL fftfwd(w0,w_hat) ! Compute Fourier transform of vorticity

      IF (diagn_flag) THEN ! If flag is true, save diagnostic quantities
        ! Saving interval for number of times to save vorticity (number of times is in the denominator)
        Nsave = (endTime/dt)/1000.0_pr; 

        IF (bin_flag) THEN ! If flag true, save bin file for adjoint solver
          CALL save_bin(w_hat)   ! Save vorticity
        END IF
 
        CALL vort2velFR(w_hat, u_hat) ! Transform vorticity to velocity in Fourier space
        Enst(1) = 0.5_pr*inner_product(w0, w0, "L2") ! 2DNS ! Compute enstrophy in physical space (could equivalently do this in Fourier space as well)
        ! Compute norms in physical space (could equivalently do this in Fourier space as well)
        InnerProduct_L2(1) = inner_product(w0, w0, "L2")   ! 2DKS
        InnerProduct_H1(1) = inner_product(w0, w0, "H1")   ! 2DKS
        InnerProduct_H2(1) = inner_product(w0, w0, "H2")   ! 2DKS
        InnerProduct_Hn1(1) = inner_product(w0, w0, "Hn1") ! 2DKS
        
        local_kin = kinetic_energy(u_hat) ! Compute the kinetic energy in Fourier space
        
        CALL MPI_REDUCE(local_kin,    KinEn(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo) ! Store kinetic energy (only on rank 0 needs a copy)
        palins(1) = palinstrophyreal(w_hat) ! Compute Palinstrophy

        IF (rank==0) THEN  ! Print information
          PRINT'(a,ES15.4)', " Initial time value      = ", t_vec(1)
          PRINT'(a,I6)',     " Number of elements      = ", numel_T
          PRINT'(a,ES16.4)', " Initial Enstrophy       = ", Enst(1)
          PRINT'(a,ES16.4)', " Initial L^2 In.Prod.    = ", InnerProduct_L2(1)
          PRINT'(a,ES16.4)', " Initial H^1 In.Prod.    = ", InnerProduct_H1(1)
          PRINT'(a,ES16.4)', " Initial H^2 In.Prod.    = ", InnerProduct_H2(1)
          PRINT'(a,ES16.4)', " Initial H^(-1) In.Prod. = ", InnerProduct_Hn1(1)
          PRINT'(a,ES16.4)', " Initial Kinetic Energy  = ", KinEn(1)
          PRINT'(a,ES16.4)', " Initial Palinstrophy    = ", palins(1)
          PRINT *, "============================================== "
        END IF
      ELSE
        ! If flag true, save bin file for adjoint solver
        IF (bin_flag) THEN ! Save vorticity in physical space, in binary format, to use in the adjoint solver
          CALL save_bin(w_hat) ! Save vorticity
        END IF

        palins(1) = palinstrophyreal(w_hat) ! Compute Palinstrophy
        ! MPI Reduce is a blocking communication, so all processes are good to proceed with timestepping
      END IF

      DO Time_iter = 2, numel_T+1 ! Main time-stepping loop
        conv0_hat = 0.0_pr ! Storage of convection term for RK method, and set to zero for first substep
        DO rk = 1, 4 ! IMEX time stepping
          CALL cal_stream(w_hat, psi_hat) ! 2DNS ! Determine the streamfunction in Fourier space: -lap(p) = w
          CALL J_bilinear_form(w_hat, psi_hat, conv_hat) ! Compute the Jacobian ! 2DNS ! Use stream function, and get the velocity and gradient of vorticity
          
          DO i2=1,local_Nx ! Compute Solution at the next step using IMEX time-stepping
            ! Compute vorticity
            DO i1=1,n_nse(2)
              !w2_hat(i1, i2) = ( ( 1.0_pr + BetaI(rk) * (-visc * ksq(i1, i2)) ) * w_hat(i1, i2) + AlphaE(rk) * conv_hat(i1, i2) &
               !                 + BetaE(rk) * conv0_hat(i1, i2) ) / ( 1.0_pr - AlphaI(rk) * (-visc * ksq(i1, i2)) ) ! 2DNS
              w2_hat(i1, i2) = ( ( 1.0_pr + BetaI(rk) * (lin_hat(i1, i2)) ) * w_hat(i1, i2) + AlphaE(rk) * conv_hat(i1, i2) &
                                + BetaE(rk) * conv0_hat(i1, i2) ) / ( 1.0_pr - AlphaI(rk) * (lin_hat(i1, i2)) ) ! 2DKS
            END DO
          END DO
          
          w_hat = w2_hat ! Update vorticity for next step ! 2DNS
          conv0_hat = conv_hat ! Update explicit part, for next substep
        END DO

        ! If flag is true, save diagnostic quantities
        IF (diagn_flag) THEN
          ! If flag true, save bin files for adjoint solver
          IF (bin_flag) THEN
            ! Save vorticity
            CALL save_bin(w_hat)
          END IF

          ! Save vorticity solution based on Nsave or every time step, if video flag active
          IF (vid_flag .AND. (mod(Time_iter, Nsave) == 1) ) THEN
          !IF (vid_flag) THEN
            ! Compute backward Fourier transform of vorticity to save in physical space
            CALL fftbwd(w_hat, w1) ! 2DNS            
            CALL save_NS_vorticity(w1, Ni, "FWD") ! Save vorticity for MATLAB analysis ! 2DNS

            ! Update video counter
            Ni = Ni + 1
          END IF
          
          local_enst = enstrophy(w_hat) ! Compute enstrophy in Fourier space
          !CALL vort2velFR(w_hat, u_hat) ! Transform vorticity to velocity in Fourier space         
          local_kin = kinetic_energy(u_hat) ! Compute the kinetic energy in Fourier space
          ! Compute inner products in Fourier space
          CALL fftbwd(w_hat, w1) ! 2DNS 
          local_L2 = inner_product(w1, w1, "L2")   ! 2DKS
          local_H1 = inner_product(w1, w1, "H1")   ! 2DKS
          local_H2 = inner_product(w1, w1, "H2")   ! 2DKS
          local_Hn1 = inner_product(w1, w1, "Hn1") ! 2DKS

          ! Store enstrophy (only on rank 0 needs a copy)
          CALL MPI_REDUCE(local_enst,   Enst(Time_iter),  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo)
          ! Store inner products (only on rank 0 needs a copy)
          CALL MPI_REDUCE(local_L2,   InnerProduct_L2(Time_iter),  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo)
          CALL MPI_REDUCE(local_H1,   InnerProduct_H1(Time_iter),  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo)
          CALL MPI_REDUCE(local_H2,   InnerProduct_H2(Time_iter),  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo)
          CALL MPI_REDUCE(local_Hn1,   InnerProduct_Hn1(Time_iter),  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo)
          ! Store kinetic energy (only on rank 0 needs a copy)
          CALL MPI_REDUCE(local_kin,    KinEn(Time_iter), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, Statinfo)
          ! Compute Palinstrophy
          palins(Time_iter) = palinstrophyreal(w_hat)
          ! MPI Reduce is a blocking communication, so all processes are good to proceed with next timestep
        ELSE
          ! If flag true, save bin files for adjoint solver
          IF (bin_flag) THEN
            ! Save vorticity
            CALL save_bin(w_hat)
          END IF
          ! Compute Palinstrophy
          palins(Time_iter) = palinstrophyreal(w_hat)
          ! MPI Reduce is a blocking communication, so all processes are good to proceed with next timestep
        END IF
      END DO
      ! Do loop adds an additional value to counter to exit loop, so we correct this
      Time_iter = Time_iter-1

      ! If flag is true, save diagnostic quantities
      IF (diagn_flag) THEN
        ! Transform to physical space
        CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
        CALL fftbwd(w_hat, w1)

        ! Compute mean of the final solution in physical space
        mean_val = mean_vort(w1)

        IF (rank==0) THEN
          PRINT'(a,ES15.4)', " Final time value       = ", t_vec(numel_T+1)
          PRINT'(a,I6)',     " Time iterations        = ", Time_iter
          PRINT'(a,ES16.4)', " Final Enstrophy        = ", Enst(Time_iter)
          PRINT'(a,ES16.4)', " Final L^2 In.Prod.     = ", InnerProduct_L2(Time_iter)
          PRINT'(a,ES16.4)', " Final H^1 In.Prod.     = ", InnerProduct_H1(Time_iter)
          PRINT'(a,ES16.4)', " Final H^2 In.Prod.     = ", InnerProduct_H2(Time_iter)
          PRINT'(a,ES16.4)', " Final H^(-1) In.Prod.  = ", InnerProduct_Hn1(Time_iter)
          PRINT'(a,ES16.4)', " Final Kinetic Energy   = ", KinEn(Time_iter)
          PRINT'(a,ES16.4)', " Final Palinstrophy     = ", palins(Time_iter)
          PRINT'(a,ES16.4)', " Mean of final solution = ", mean_val
          PRINT *, "============================================== "
        END IF

      ELSE
        ! If flag true, print values
        IF (bin_flag) THEN
          ! Transform to physical space
          CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
          CALL fftbwd(w_hat, w1) ! 2DNS
          ! CALL fftbwd(u_hat, w1) ! 2DKS
          ! Compute mean of the final solution in physical space
          mean_val = mean_vort(w1)

          IF (rank==0) THEN
            PRINT'(a,ES15.4)', " Final time value       = ", t_vec(numel_T+1)
            PRINT'(a,I6)',     " Time iterations        = ", Time_iter
            PRINT'(a,ES16.4)', " Final Palinstrophy     = ", palins(Time_iter)
            PRINT'(a,ES16.4)', " Mean of final solution = ", mean_val
            PRINT *, "============================================== "
          END IF

        END IF
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo) ! Wait to exit until all processors finished
    END SUBROUTINE IMEX_fwd

    !=================================================================================
    ! Backward solver with third order IMEX Time stepping
    ! Output: z_hat - adjoint vorticity field (at time 0) in Fourier space
    !! Output: z - adjoint vorticity field (at time 0) in physical space
    !=================================================================================
    SUBROUTINE IMEX_bwd(z_hat)
      ! Load variables
      USE global_variables
      ! Load subroutines
      USE fftwfunction        ! FFT routines
      USE function_ops        ! Contains functions/subroutines to operate on data
      USE data_ops            ! Contains functions/subroutines that read/write data
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      ! Initialize variables
!      REAL(pr),    DIMENSION(:,:),      INTENT(OUT)  :: z                   ! Initial vorticity field
!      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: z_hat, z2_hat       ! Complex adjoint vorticity, its Fourier transform, and RK storage term
      COMPLEX(pr), DIMENSION(:,:),      INTENT(OUT)  :: z_hat               ! Initial vorticity field
      REAL(pr),    DIMENSION(1:n_nse(1),1:local_Ny)  :: z                   ! Initial vorticity field
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: z2_hat              ! Complex adjoint vorticity, its Fourier transform, and RK storage term
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: w_hat, w1_hat       ! Complex vorticity storage from the forward solution
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: psi_hat, wpsi_hat   ! Adjoint and vorticity streamfunction in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: nonlin_hat          ! Adjoint nonlinear terms in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: nonlin0_hat         ! Storage for RK adjoint nonlinear terms in Fourier space
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx)  :: conv_hat            ! Convection term in Fourier space
      REAL(pr),    DIMENSION(1:n_nse(1),1:local_Ny)  :: w, wx, wy, u, v     ! Vorticity, partial derivatives, and velocity components from the forward solution
      INTEGER                                        :: rk, i1, i2          ! Integers for looping through values
      REAL(pr)                                       :: mean_val            ! Scalar for storing the mean of vorticity
      ! Set the terminal condition for adjoint system as the zero solution
      z_hat   = 0.0_pr
      psi_hat = 0.0_pr
      ! Read the terminal vorticity solution
      CALL read_bin(w_hat)
      ! Dealias (ensure no round-off errors were introduced from read/write)
      CALL dealiasing(w_hat)

      ! Ensure all processes have completed before proceeding with timestepping method
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

      ! Main time-stepping loop
      DO Time_iter = numel_T, 1, -1
        ! Storage of nonlinear terms for RK method, and set to zero for first substep
        nonlin0_hat = 0.0_pr
        ! Obtain the previous vorticity solution, to use for nonlinear term (previously treated implicitly, so it is now treated explicitly)
        w1_hat = w_hat
        ! Read the vorticity solution, for RHS calculations (treated implicitly)
        CALL read_bin(w_hat)
        ! Dealias (ensure no round-off errors were introduced from read/write)
        CALL dealiasing(w_hat)

        ! Determine the forward streamfunction in Fourier space: -lap(p) = w
        CALL cal_stream(w1_hat, wpsi_hat)
        ! Determine the components of the forward velocity field
        CALL cal_deriv(wpsi_hat, v, u)
        ! Partial derivatives of vorticity solution (explicitly treated)
        CALL cal_deriv(w1_hat, wx, wy)

        ! IMEX time stepping
        DO rk = 1, 4
          ! Determine adjoint convection term
          CALL adj_conv(z_hat, u, v, conv_hat)
          ! Add adjoint streamfunction to nonlinear term, since it is treated explicity in calculations
          nonlin_hat = conv_hat - psi_hat

          ! Compute Solution at the next step using IMEX time-stepping
          DO i2=1,local_Nx
            DO i1=1,n_nse(2)
              ! Compute vorticity
              ! Note: RHS of the adjoint equation is simply a constant times the Laplacian of the vorticity solution
              ! hence the additional implicitly treated w_hat term independent of the RK scheme (dt/4 to account for the 4 step method).
              z2_hat(i1, i2) = (  ( 1.0_pr + BetaI(rk) * (-visc * ksq(i1, i2)) ) * z_hat(i1, i2) + &
                                  (dt/4.0_pr) * ( ksq(i1, i2) ) * w_hat(i1, i2) + &
                                  AlphaE(rk) * nonlin_hat(i1, i2) + &
                                  BetaE(rk) * nonlin0_hat(i1, i2) ) / &
                                  ( 1.0_pr - AlphaI(rk) * (-visc * ksq(i1, i2)) )
            END DO
          END DO
          ! Update vorticity for next step
          z_hat = z2_hat
          ! Transform back to physical space
          CALL fftbwd(z_hat, z)

          ! Determine adjoint streamfunction
          CALL adj_stream(z, wx, wy, psi_hat)

          ! Update explicit part, for next substep
          nonlin0_hat = nonlin_hat
        END DO
      END DO

      IF (rank==0) THEN
        PRINT'(a,I6)',    " Time iteration = ", Time_iter
        PRINT *, "============================================== "
      END IF
      ! Compute mean of the final solution in physical space
      mean_val = mean_vort(z)
      IF (rank==0) PRINT'(a,ES16.4)', " Mean of final adjoint solution = ", mean_val
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo) ! Wait to exit until all processors finished
    END SUBROUTINE IMEX_bwd

END MODULE

