!=====================================================================
! MODULE CONTAINS ROUTINES FOR THE OPTIMIZATION METHOD.
!
! Author: Pritpal Matharu
! Department of Mathematics and Statistics            
! McMaster University
! Date: 2020/12/24
!
! CONTAINS
! (*) OptScheme   - Determines optimal IC via optimization scheme
! (*) KappaTest   - Performs the Kappa Test
! (*) mnbrak      - Bracketing method, to maximize cost functional
! (*) brent       - Brent method, for determining optimal step length
! (*) costfuneval - Determines value of cost functional
!=====================================================================
MODULE optimization
  USE global_variables, ONLY: pr
  IMPLICIT NONE ! Prevent using implicit typing rules for entire module

  ! Set parameters for optimization routines
  REAL(pr), PARAMETER :: TOL      = 1e-6         ! Tolerance for optimization
  INTEGER, PARAMETER  :: max_iter = 1000         ! Maximum number of iterations for optimization scheme
  INTEGER, PARAMETER  :: fq       = 10           ! Frequency for clearing conjugate gradient method (clearing Polak-Ribière beta value)

  REAL(pr), PARAMETER :: GOLD  = 1.0_pr + ((1.0_pr + SQRT(5.0_pr))/2.0_pr) ! Golden ratio plus one, for expanding interval
  REAL(pr), PARAMETER :: CGOLD = 1.0_pr/GOLD                               ! Inverse golden ration, for shrinking interval
  REAL(pr), PARAMETER :: TOLBR = 1E-10                                      ! Relative error tolerance for bracketing method
  INTEGER,  PARAMETER :: ITMAX = 200                                       ! Max number of iterations for bracketing and brent method
  REAL(pr), PARAMETER :: TOLBT = 1E-4                                      ! Relative error tolerance for brent method
  REAL(pr), PARAMETER :: ZEPS  = 1E-4                                      ! Positive absolute error tolerance. Small number that protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero

  CONTAINS

    !===========================================================================
    ! Perform optimization scheme, to find the IC which maximize cost functional
    ! Input: w0 - Initial condition
    !===========================================================================
    SUBROUTINE OptScheme(w0)
      ! Load variables
      USE global_variables
      USE solvers             ! Contains the time stepping methods, to solve the PDES
      USE fftwfunction        ! FFT routines
      USE function_ops        ! Contains functions/subroutines to operate on data
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      USE data_ops            ! Contains functions/subroutines that modify data or read/write data

      ! Input initial guess of vorticity field
      REAL(pr), DIMENSION(:,:),       INTENT(IN)    :: w0                     ! Initial condition
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx) :: gradL2_hat             ! Fourier transform of L2 Gradient
      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny)    :: del_J, delk, PolRib, pPolRib   ! Gradient, previous iteration, conjugate gradient, and previous conjugate gradient (L2 or Sobolev, depending on Grad_type parameter)
      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny)    :: w1                     ! Optimal initial condition
      INTEGER                                       :: i1, i2, ii             ! Temporary integers for loops
      REAL(pr)                                      :: J1, J2                 ! Value of cost functional and previous cost functional
      REAL(pr)                                      :: mean_val, constr_val   ! Scalar for storing the mean of vorticity and constraint value
      REAL(pr)                                      :: tau0                   ! Starting step length for line search
      REAL(pr)                                      :: bPR                    ! Starting step length for line search
      INTEGER                                       :: iter                   ! Iteration counter

      REAL(pr), DIMENSION(1:max_iter)     :: Jiter ! Vector for storing values of cost functional dynamically
      REAL(pr), DIMENSION(:), ALLOCATABLE :: Jsave ! Vector for saving cost functional at final stage
      ! Set intial values for optimization scheme
      tau0     = 1.0_pr
      iter     = 1

      ! Verify constraint and mean of the initial condition in physical space
      constr_val = inner_product(w0, w0, normconstr)
      mean_val   = mean_vort(w0)
      IF (rank==0) PRINT'(a,ES16.4)', " Mean of IC       = ", mean_val
      IF (rank==0) PRINT'(a,ES16.4)', " Constraint value = ", constr_val
      ! MPI Reduce is a blocking communication, so all processes are good to proceed with timestepping
      IF (rank==0) PRINT *, " Solving Forward DNS..."
      ! Solve forward 2D DNS
      Time_iter = 1 ! Time iteration counter
      CALL IMEX_fwd(w0, Palin)
      bin_flag   = .FALSE.
      diagn_flag = .FALSE.
      IF (rank==0) PRINT *, " DNS solved."
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

      ! Compute value of cost functional
      J1 = trapzT(Palin)
      J2 = 0.0_pr ! Dummy value
      IF (rank==0) PRINT *, " Initial Cost Functional = ", J1
      ! Set initial condition to be iteration 1 IC
      w1 = w0

      ! Store values for initial guess
      Jiter(1) = J1

      !! Loop iteration
      ! Iterate until condition is met or maximum iteration is reached
      DO WHILE ( (ABS(J2 - J1)/ABS(J1) > TOL) .AND. (iter < max_iter) )
        IF (rank==0) PRINT *, " Solving Adjoint System..."
        ! Solve adjoint system
        CALL IMEX_bwd(gradL2_hat)
        IF (rank==0) PRINT *, " Adjoint solved."
        ! Compute the H1 Sobolev gradient
        CALL SobolevGrad(gradL2_hat, del_J)
        IF (rank==0) PRINT *, " Sobolev solved."

        ! *** Just for testing size of gradient
        constr_val = inner_product(del_J, del_J, "L2")
        IF (rank==0) PRINT'(a,ES16.4)', " L2 size of grad = ", constr_val
        constr_val = MAXVAL(ABS(del_J))
        IF (rank==0) PRINT'(a,ES16.4)', " sup-norm of grad = ", constr_val

        !! Conjugate Gradient method
        ! Use Conjugate gradient method, with the Polak-Ribiere method including frequency resetting
        IF ( (iter >= 2) .AND. (MOD(iter, fq)  /= 0) ) THEN
          IF (rank==0) PRINT *, " Solving bPR..."
          ! Using the Polak Ribere Method
          bPR = ( inner_product(del_J, del_J - delk, Grad_type) )/( inner_product(delk, delk, Grad_type) )

          ! Use value to create the conjugate gradient
          PolRib = del_J + bPR*pPolRib;
        ELSE
          ! Frequency clearing to reset the conjugate-gradient procedure
          bPR = 0.0_pr
          ! Use gradient (previous gradient cleared)
          PolRib = del_J

          IF (rank==0) PRINT *, " Saving checkpoint."
          ! Allocate array to save cost functional values
          ALLOCATE( Jsave(1:iter) )
          ! Allocate values of cost functional
          Jsave = Jiter(1:iter)
          CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
          ! Save values from optimization
          CALL save_NS_Opt(w0, w1, Palin, Enst, KinEn, t_vec, Jsave)
          ! Deallocate array
          DEALLOCATE(Jsave)
        END IF

        ! *** Just for testing size of gradient
        constr_val = inner_product(PolRib, PolRib, "L2")
        IF (rank==0) PRINT'(a,ES16.4)', " L2 size of grad = ", constr_val


        ! Set Cost functional equal to current iteration
        IF (iter /= 1) THEN ! Ensure that the cost does not equal zero (only for the first iteration)
          J1 = J2
        END IF

        ! Initialize bracketing procedure
        IF (rank==0) PRINT *, " Start bracketing."
        ! Bracketing the maximum from the right, for the Brent method
        CALL mnbrak(w1, PolRib, J1, tau0, J2)

        ! Determine optimal step length along the gradient, using line search method
        CALL brent(w1, PolRib, 0.0_pr, tau0, -J2, tau0, J2)
        J2 = -J2 ! Correct sign output from brent method
        IF (rank==0) PRINT'(a,ES16.4)', " Cost Functional = ", J2

        ! Ensure that the new value is larger than previous iteration. If not, force to exit optimization loop
        IF (J1 > J2) THEN
          IF (rank==0) PRINT *, " Can't bracket get a better value with current parameters"
          J2   = J1     ! This will force it to exit optimization loop

          ! Don't update initial condition, as this will force the optimization to exit optimization
          ELSE
          ! Update initial condition
          w1 = w1 + tau0*PolRib
        END IF

        ! Store values for the Polak Ribere method
        delk     = del_J;
        pPolRib  = PolRib

        ! Enforce the norm constraint
        CALL ProjConstr(w1)

        ! Verify constraint and mean of the initial condition in physical space
        constr_val = inner_product(w1, w1, normconstr)
        mean_val   = mean_vort(w1)
        IF (rank==0) PRINT'(a,ES16.4)', " Mean of IC       = ", mean_val
        IF (rank==0) PRINT'(a,ES16.4)', " Constraint value = ", constr_val
        ! See if we are going into another loop or going to be exiting
        IF ( (ABS(J2 - J1)/ABS(J1) <= TOL) .OR. (iter >= max_iter) ) THEN
          IF (rank==0) PRINT *, " Solving optimal DNS for saving..."
            Time_iter = 1       ! Reset time iteration counter
            diagn_flag = .TRUE.       ! Save diagnostics
          ELSE
            IF (rank==0) PRINT *, " Solving Forward DNS..."
            Time_iter = 1       ! Reset time iteration counter
            bin_flag  = .TRUE.  ! Save bin files for adjoint solver
        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

        ! Solve forward 2D DNS
        CALL IMEX_fwd(w1, Palin)
        bin_flag  = .FALSE. ! Turn off bin save
        IF (rank==0) PRINT *, " DNS solved."

        ! Compute value of cost functional
        J2 = trapzT(Palin)
        CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
        IF (rank==0) PRINT'(a,I5)',     " Iteration       = ", iter
        IF (rank==0) PRINT'(a,ES16.4)', " Cost Functional = ", J2

        ! Increment iteration counter
        iter = iter + 1
        ! Save cost functional value
        Jiter(iter) = J2
      END DO

      IF (rank==0) PRINT *, " Optimal IC solved, saving results."

      ! Allocate array to save cost functional values
      ALLOCATE( Jsave(1:iter) )
      ! Allocate values of cost functional
      Jsave = Jiter(1:iter)
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
      ! Save values from optimization
      CALL save_NS_Opt(w0, w1, Palin, Enst, KinEn, t_vec, Jsave)
      ! Deallocate array
      DEALLOCATE(Jsave)
    END SUBROUTINE OptScheme

    !==========================================================================
    ! *** Bracket the optimal step length, to maximize cost functional ***
    ! Input:     w    - current initial condition, in physical space
    !         grad    - current gradient
    !           J1    - previous value of cost functional
    !          tau    - current step length
    ! Output:  tau    - bracketed step length
    !           J2    - corresponding value of cost functional
    !==========================================================================
    SUBROUTINE mnbrak(w, grad, J1, tau, J2)
      ! Load variables
      USE global_variables, ONLY: pr, rank, MACH_EPSILON
      ! Initialize variables
      REAL(pr), DIMENSION(:,:), INTENT(IN)    :: w, grad  ! Current iteration of initial condition and gradient
      REAL(pr),                 INTENT(IN)    :: J1       ! Previous iteration cost functional
      REAL(pr),                 INTENT(INOUT) :: tau      ! Starting step length, and bracketed step length
      REAL(pr),                 INTENT(OUT)   :: J2       ! Value of cost functional at current step length/ bracketed cost functional
      REAL(pr)                                :: J3       ! Value of cost functional at previous step length
      INTEGER                                 :: FuncEval ! Counter to keep track of function evaluations

      ! Start counter for function evaluations
      FuncEval = 1
      ! Determine the cost functional for the given step length and gradient
      J2 = costfuneval(w, grad, tau, 0)
      J3 = 0.0_pr ! Dummy value

      ! Print value of cost functional
      if (rank==0) then
        print *, "      brak; do iteration FuncEval =", FuncEval
        print *, "      brak;                   tau =", tau
        print *, "      brak;                    J2 =", J2
      end if

      ! Modifying the interval for step lengths
      IF (J2 > J1) THEN
        ! If the cost functional at tau is greater than the previous iteration, we try to expand the bracket interval, to ensure we include the best possible stepsize
        DO WHILE ( (J2 > J1) .AND. (ABS(J3 - J2)/ABS(J2) > TOLBR) .AND. (FuncEval < ITMAX) )

          ! Store for check criterion
          J3 = J2

          ! Expand by the constant bracketing factor (note constant above, aggressive factor to reduce number of function evaluations here)
          tau = GOLD*tau
          ! Determine the cost functional for the given step length
          J2 = costfuneval(w, grad, tau, 0)
          FuncEval = FuncEval+1 ! Update the counter for number of function evaluations
          if (rank==0) then
            print *, "      brak; do iteration FuncEval =", FuncEval
            print *, "      brak;                   tau =", tau
            print *, "      brak;                    J2 =", J2
          end if
        END DO
      ELSE
        ! If the cost functional at tau is less than the previous iteration, we try to shrink the bracket interval, to reduce the number of evaluations in the line search
        DO WHILE ( (J2 < J1) .AND. (tau > MACH_EPSILON) .AND. (FuncEval < ITMAX) )
          ! Shrink by the constant bracketing factor (note constant above, aggressive factor to reduce number of function evaluations here)
          tau = CGOLD*tau
          ! Determine the cost functional for the given step length
          J2 = costfuneval(w, grad, tau, 0)
          FuncEval = FuncEval+1 ! Update the counter for number of function evaluations
          if (rank==0) then
            print *, "      brak; do iteration FuncEval =", FuncEval
            print *, "      brak;                   tau =", tau
            print *, "      brak;                    J2 =", J2
          end if
        END DO
        ! Ensure that we did not shrink the bracket too much!
        tau = GOLD*tau
      END IF

      ! Print error if max iterations reached
      IF (FuncEval==ITMAX) THEN
        IF (rank==0) PRINT *, " Maximum number of iterations reached in bracketing function!!!"
      ELSEIF ( ABS(J3 - J2)/ABS(J2) < TOLBR ) THEN
        IF (rank==0) PRINT *, " Bracketing method didn't really increase things...."
      END IF

    END SUBROUTINE mnbrak


    !==========================================================================
    ! *** Brent Algorithm for line minimization (note function calls) ***
    ! Based on Section 10.2-10.3 of Press et al. (2007), codes by Brent:
    ! https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
    ! Input:    phi   - current initial condition, in physical space
    !          grad   - current gradient
    !          tauA   - left end of bracketed step length interval
    !          tauB   - right end of bracketed step length interval
    !            FB   - cost functional at right end of step length bracket
    ! Output:     X   - optimal step length
    !            FX   - value of cost functional at optimal step length
    !==========================================================================
!    FUNCTION brent(phi, grad, tauA, tauB) RESULT(X)
    SUBROUTINE brent(phi, grad, tauA, tauB, FB, X, FX)
      ! Load variables
      USE global_variables, ONLY: pr, rank
      ! Initialize variables
      REAL(pr), DIMENSION(:,:), INTENT(IN)  :: phi, grad          ! Current value of optimal control variable and the current gradient
      REAL(pr),                 INTENT(IN)  :: tauA, tauB         ! Bracket of interval to search within for optimal step length
      REAL(pr),                 INTENT(IN)  :: FB                 ! Cost functional at right endpoint of bracket
      REAL(pr),                 INTENT(OUT) :: X                  ! Optimal step length and cost functional
      REAL(pr),                 INTENT(OUT) :: FX                 ! Optimal step length and cost functional
      REAL(pr) :: D, A, B, V, W, E, ETEMP, P, Q, R, U, XM         ! Temporary variable for determining optimal step length
      REAL(pr) :: FV, FW, FU                                      ! Temporary values of the functional
      REAL(pr) :: TOL1, TOL2                                      ! Tolerances
      INTEGER  :: FLAG, j                                         ! Flag for determining parabolic interpolation or golden section search, and iteration counter

      D = 0.0_pr          ! Variable for distance for new trial step
      A = MIN(tauA,tauB)  ! Left endpoint, for interval
      B = MAX(tauA,tauB)  ! Right endpoint, for interval
      X = tauB*CGOLD      ! Value of (previous) trial step (should be A + C(B-A), but A=0)
      W = X               ! Variable for parabolic fit
      V = X               ! Variable for parabolic fit
      E = 0.0_pr          ! Distance from old step to left or right endpoint

      ! Value of functional at trial step
      FX = costfuneval(phi, grad, X, 1)

      ! Store value for reference in loop
      FV = FX
      FW = FX
      ! Main loop (note stopping criterion in loop)
      DO j=1,ITMAX
        if (rank==0) then
          print *, "      brent; do iteration j =", j
          print *, "      brent;              X =", X
          print *, "      brent;             FX =", -FX
        end if

        XM = 0.5_pr*(A+B)         ! Midpoint of current interval
        TOL1 = TOLBT*ABS(X)+ZEPS    ! Tolerance for parabolic fit
        TOL2 = 2.0_pr*TOL1        ! Tolerance for stopping criterion

        ! Stopping criterion for main loop
        IF ( ABS(X-XM) <= (TOL2-0.5*(B-A)) ) EXIT
        ! Flag to default a golden step
        FLAG = 1
        ! Fit a parabola if satisfied
        IF ( ABS(E) > TOL1 ) THEN
          ! Construct a trial parabolic fit
          R = (X-W)*(FX-FV)
          Q = (X-V)*(FX-FW)
          P = (X-V)*Q - (X-W)*R
          Q = 2.0_pr*(Q-R)
          IF ( Q > 0.0_pr ) P = -P
          Q = ABS(Q)
          ETEMP = E
          E = D
          ! Determine the acceptability of the parabolic fit
          IF ( (ABS(P) >= ABS(0.5_pr*Q*ETEMP)) .OR. (P <= Q*(A-X)) .OR. (P >= Q*(B-X)) ) THEN
            ! Flag to do a golden step
            FLAG = 1
          ELSE
            ! Flag to do a parabolic interpolation fit
            FLAG = 2
          END IF
          if (rank==0) then
            print *, "      brent; do iteration; case 1; FLAG =", FLAG
          end if
        END IF

        ! Either perform golden section search or parabolic interpolation
        SELECT CASE (FLAG)
          CASE (1)
            ! A golden-section step
            IF (X >= XM) THEN
              if (rank==0) then
                print *, "      brent; do iteration; golden step, from A"
              end if
              E = A-X
            ELSE
              if (rank==0) then
                print *, "      brent; do iteration; golden step, from B"
              end if
              E=B-X
            END IF
            D = CGOLD*E
          CASE (2)
            ! Parabolic interpolation fit
            if (rank==0) then
              print *, "      brent; do iteration; parabolic interpolation step "
            end if
            ! Take parabolic interpolation step
            D = P/Q
            U = X+D
            ! Must ensure that functional is not evaluated too close to end points
            IF ( (U-A < TOL2) .OR. (B-U < TOL2) ) D = SIGN(TOL1, XM-X)
        END SELECT
        ! Ensure that functional is not evaluated too close to X
        IF ( ABS(D) >= TOL1 ) THEN
          U = X+D
        ELSE
          U = X + SIGN(TOL1,D)
        END IF

        ! Update value of functional
        FU = costfuneval(phi, grad, U, 1)

        ! Further shrink the search interval, and update values of A, B, V, W, and X, in addition to respective functionals
        IF ( FU <= FX ) THEN
          IF ( U >= X ) THEN
            if (rank==0) then
              print *, "      brent; do iteration; left interval is now old X step"
            end if
            A = X
          ELSE
            if (rank==0) then
              print *, "      brent; do iteration; right interval is now old X step"
            end if
            B = X
          END IF
          ! Update values
          V = W
          FV = FW
          W = X
          FW = FX
          X = U
          FX = FU
        ELSE
          IF ( U < X ) THEN
            if (rank==0) then
              print *, "      brent; do iteration; left interval is now new U step"
            end if
            A = U
          ELSE
            if (rank==0) then
              print *, "      brent; do iteration; right interval is now new U step "
            end if
            B = U
          END IF
          ! Perform housekeeping
          IF ( (FU <= FW) .OR. (W == X) ) THEN
            if (rank==0) then
              print *, "      brent; do iteration; housekeeping 1 "
            end if
            V = W
            FV = FW
            W = U
            FW = FU
          ELSEIF ( (FU <= FV) .OR. (V==X) .OR. (V==W)) THEN
            if (rank==0) then
              print *, "      brent; do iteration; housekeeping 2 "
            end if
            V = U
            FV = FU
          END IF
        END IF ! End of updating values

      END DO ! End of main loop

      ! Check if right bracket value if better
      IF ( FB < FX) THEN
        if (rank==0) then
          print *, "      brent; right bracket better! "
        end if
        ! Set right endpoint as optimal step length (algorithm cannot evaluate at this point)
        X = tauB
      END IF
    END SUBROUTINE brent
!    END FUNCTION brent

    !==========================================================================
    ! *** Calculate value of cost functional for a given step length ***
    ! Input:     f    - current initial condition, in physical space
    !         grad    - current gradient
    !          tau    - given step length
    !      minflag    - flag for indicating brent method
    ! Output:    J    - value of cost functional
    !==========================================================================
    FUNCTION costfuneval(f, grad, tau, minflag) RESULT (J)
      ! Load variables
      USE global_variables, ONLY: pr, n_nse, local_Ny, Palin
      ! Load subroutines
      USE solvers             ! Contains the time stepping methods, to solve the PDES
      USE function_ops        ! Contains functions/subroutines to operate on data
      ! Initialize variables
      REAL(pr), DIMENSION(:,:),       INTENT(IN) :: f,grad    ! Current initial condition and gradient
      REAL(pr),                       INTENT(IN) :: tau       ! Given step length
      INTEGER,                        INTENT(IN) :: minflag   ! Flag for brent method
      REAL(pr)                                   :: J         ! Cost functional
      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny) :: w         ! Storage for initial condition

      ! Apply given step length and gradient
      w = f+tau*grad
      ! Enforce the norm constraint
      CALL ProjConstr(w)
      ! Solve PDE to determine the cost functional for the given step length and gradient
      CALL IMEX_fwd(w, Palin)

      ! Compute value of cost functional
      SELECT CASE (minflag)
        CASE (1)
          ! If called from brent method, use negative value for "minimization"
          J = -trapzT(Palin)
        CASE DEFAULT
          ! Otherwise, compute regularly
          J = trapzT(Palin)
      END SELECT


    END FUNCTION costfuneval

    !===========================================================================
    ! *** Perform Kappa test to verify gradient expression ***
    ! Input: w0 - Initial condition
    !===========================================================================
    SUBROUTINE KappaTest(w0)
      ! Load variables
      USE global_variables
      USE solvers             ! Contains the time stepping methods, to solve the PDES
      USE fftwfunction        ! FFT routines
      USE function_ops        ! Contains functions/subroutines to operate on data
      USE mpi                 ! Use MPI module (binding works well with fftw libraries)
      USE data_ops            ! Contains functions/subroutines that modify data or read/write data

      ! Input initial guess of vorticity field
      REAL(pr), DIMENSION(:,:),       INTENT(IN)    :: w0                     ! Initial condition
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx) :: gradL2_hat             ! Fourier transform of L2 Gradient
      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny)    :: grad_H1                ! L2 Gradient
      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny)    :: wpert                  ! Perturbation
      REAL(pr), DIMENSION(1:n_nse(1),1:local_Ny)    :: wprime                 ! Perturbation system field
      COMPLEX(pr), DIMENSION(1:n_nse(2),1:local_Nx) :: w_hat                  ! Fourier transform of complex vorticity
      REAL(pr)                                      :: kap_denom, kap_numer   ! Numerator and denominator for Kappa Test
      REAL(pr)                                      :: J, J1                  ! Value of cost functional and perturbation cost functional
      REAL(pr)                                      :: ep                     ! Epsilon value for Kappa test
      REAL(pr)                                      :: X, Y                   ! Scalars for storing x and y coordinates
      REAL(pr)                                      :: mean_val               ! Scalar for storing the mean of vorticity
      INTEGER                                       :: i1, i2, ii             ! Temporary integers for loops


      ! Compute mean of the initial condition in physical space
      mean_val = mean_vort(w0)
      IF (rank==0) PRINT'(a,ES16.4)', " Mean of IC = ", mean_val
      ! MPI Reduce is a blocking communication, so all processes are good to proceed with timestepping
      IF (rank==0) PRINT *, " Solving Forward DNS..."
      ! Solve forward 2D DNS
      Time_iter = 1 ! Time iteration counter
      CALL IMEX_fwd(w0, Palin)
      bin_flag = .TRUE.
      IF (rank==0) PRINT *, " DNS solved."
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

      ! Compute value of cost functional
      J = Palin(dt_save+1)
      IF (rank==0) PRINT'(a,ES16.4)', " Initial Cost Functional = ", J

      !!! Perturbation for Kappa Test
      IF (rank==0) PRINT *, " *** Initial Condition Perturbation ***"
      wpert = w0

      IF (rank==0) PRINT *, " Solving Adjoint System..."
      ! Solve adjoint system
      CALL IMEX_bwd(gradL2_hat)

      ! Compute the H1 Sobolev gradient
      CALL SobolevGrad(gradL2_hat, grad_H1)
      ! L2 denominator for Kappa Test
      kap_denom = inner_product(grad_H1, wpert, Grad_type)
      IF (rank==0) PRINT'(a,ES16.4)', " Denominator of Kappa Test = ", kap_denom
      CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)

      ! Loop through values of epsilon for Kappa Test
      DO ii = -15, -1, 1
        ! Value of epsilon
        ep = 10.0_pr**(ii)
        ! Initial condition with perturbation
        wprime = w0 + ep*wpert

        ! Compute mean of the initial condition in physical space
        mean_val = mean_vort(wprime)
        IF (rank==0) PRINT'(a,ES16.4)', " Mean of perturbation = ", mean_val

        ! Solve perturbation system
        CALL IMEX_fwd(wprime, Palin)

        ! Compute value of cost functional
        J1 = Palin(dt_save+1)

        ! Kappa Test numerator
        kap_numer = (J1 - J)/ep
        IF (rank==0) THEN
          PRINT'(a,ES16.4)', " Epsilon                   = ", ep
          PRINT'(a,ES16.4)', " Initial Cost Functional   = ", J
          PRINT'(a,ES16.4)', " Perturbed Cost Functional = ", J1
          PRINT'(a,ES16.4)', " Kappa Test Numerator      = ", kap_numer
          PRINT'(a,ES20.8)', " Kappa Value               = ", kap_numer/kap_denom
          PRINT *, "============================================== "
        END IF
      END DO

    END SUBROUTINE KappaTest

END MODULE
